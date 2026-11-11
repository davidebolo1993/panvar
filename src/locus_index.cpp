#include "panvar/locus_index.hpp"

#include "panvar/anchors.hpp"
#include "panvar/background.hpp"
#include "panvar/background_graph.hpp"
#include "panvar/fasta_index.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/gz_reader.hpp"
#include "panvar/kmer.hpp"
#include "panvar/ref_path.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace panvar::recruit {

namespace {

namespace fs = std::filesystem;

constexpr std::size_t kK = 31, kS = 11;
constexpr std::uint64_t kMagic = 0x3130305850585650ULL;   // "PVPXP001"
constexpr std::uint64_t kVersion = 3;   // 2: decoy windows stop at the locus; 3: clusters split at it
const std::size_t kLengths[] = {100, 150, 250, 1000, 10000};
constexpr std::size_t kWindowStep = 25;

std::string hex64(std::uint64_t v) {
    char b[17];
    std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(v));
    return b;
}

std::string file_checksum(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("recruit: cannot read " + path);
    std::vector<char> buf(1 << 20);
    std::uint64_t h = 0x9ae16a3b2f90404fULL, n = 0;
    while (in) {
        in.read(buf.data(), static_cast<std::streamsize>(buf.size()));
        const std::size_t got = static_cast<std::size_t>(in.gcount());
        std::size_t i = 0;
        for (; i + 8 <= got; i += 8) { std::uint64_t w; std::memcpy(&w, buf.data() + i, 8); h = kmer::mix64(h ^ w); }
        std::uint64_t w = 0;
        std::memcpy(&w, buf.data() + i, got - i);
        h = kmer::mix64(h ^ w ^ (static_cast<std::uint64_t>(got) << 48));
        n += got;
    }
    return hex64(kmer::mix64(h ^ n));
}

bool is_gfa_file(const std::string& path) {
    GzLineReader in(path);
    std::string line;
    while (in.ok() && in.getline(line)) {
        if (line.empty()) continue;
        return line.size() > 1 && line[1] == '\t' && std::string("HSLPWCJ#").find(line[0]) != std::string::npos;
    }
    return false;
}

std::string locus_name_of(const std::string& gfa) {
    std::string n = fs::path(gfa).filename().string();
    for (const char* ext : {".gz", ".gfa", ".gfa1"})
        if (n.size() > std::strlen(ext) && n.compare(n.size() - std::strlen(ext), std::strlen(ext), ext) == 0) n.resize(n.size() - std::strlen(ext));
    return n;
}

// Which counts to read: the compacted dictionary beside the store when it is current, else the store.
std::unique_ptr<bg::Background> open_counts(const std::string& store, std::string& kind, const std::function<void(const std::string&)>& log) {
    std::string base = store;
    while (base.size() > 1 && base.back() == '/') base.pop_back();
    if (fs::is_regular_file(base + ".sshash")) {
        try {
            auto b = bg::Background::open(base + ".sshash");
            kind = "compacted";
            return b;
        } catch (const std::exception& e) {
            if (log) log(std::string("not using the compacted background: ") + e.what());
        }
    }
    kind = "store";
    return bg::Background::open(store);
}

struct LocusWork {
    LocusIndex idx;
    std::vector<std::uint64_t> counts;                 // background counts per k-mer
    std::vector<std::uint64_t> sub, sub_trim;          // own occurrences in the extended / trimmed spans
    std::map<std::string, std::vector<anchors::Span>> spans;   // per genome
};

// Every k-mer of every path: nodes it starts on, largest copy number within a path, paths holding it.
void collect_kmers(const Graph& g, LocusIndex& idx, std::vector<std::uint32_t>& present) {
    std::vector<std::string> ids;
    for (const auto& [id, n] : g.nodes) ids.push_back(id);
    const bool numeric = std::all_of(ids.begin(), ids.end(), [](const std::string& s) {
        return !s.empty() && s.size() < 19 && s.find_first_not_of("0123456789") == std::string::npos;
    });
    std::sort(ids.begin(), ids.end(), [&](const std::string& a, const std::string& b) {
        return numeric ? std::stoull(a) < std::stoull(b) : a < b;
    });
    std::unordered_map<std::string, std::uint32_t> node_index;
    for (std::size_t i = 0; i < ids.size(); ++i) node_index[ids[i]] = static_cast<std::uint32_t>(i);
    idx.nodes = ids;

    struct Agg { std::uint32_t mult = 0, paths = 0; std::vector<std::uint32_t> nodes; };
    std::unordered_map<std::uint64_t, Agg> agg;
    std::vector<std::pair<std::uint64_t, std::uint32_t>> occ;
    for (const auto& p : g.paths) {
        const std::string seq = spell_path_steps_sequence(g, p.steps);
        std::vector<std::uint64_t> step_end;
        std::uint64_t pos = 0;
        for (const auto& st : p.steps) { pos += g.nodes.at(st.node_id).sequence.size(); step_end.push_back(pos); }
        occ.clear();
        std::size_t si = 0;
        kmer::for_each_canonical_kmer(seq.data(), seq.size(), kK, [&](std::uint64_t c, std::size_t start, bool) {
            while (si < step_end.size() && step_end[si] <= start) ++si;
            occ.emplace_back(c, node_index.at(p.steps[si].node_id));
        });
        std::sort(occ.begin(), occ.end());
        for (std::size_t i = 0; i < occ.size();) {
            std::size_t e = i;
            Agg& a = agg[occ[i].first];
            while (e < occ.size() && occ[e].first == occ[i].first) {
                if (a.nodes.empty() || std::find(a.nodes.begin(), a.nodes.end(), occ[e].second) == a.nodes.end()) a.nodes.push_back(occ[e].second);
                ++e;
            }
            a.mult = std::max<std::uint32_t>(a.mult, static_cast<std::uint32_t>(e - i));
            ++a.paths;
            i = e;
        }
    }
    idx.codes.reserve(agg.size());
    for (const auto& [c, a] : agg) idx.codes.push_back(c);
    std::sort(idx.codes.begin(), idx.codes.end());
    std::map<std::vector<std::uint32_t>, std::uint32_t> sets;
    idx.nodeset_first.push_back(0);
    for (const std::uint64_t c : idx.codes) {
        Agg& a = agg[c];
        std::sort(a.nodes.begin(), a.nodes.end());
        auto [it, fresh] = sets.emplace(a.nodes, static_cast<std::uint32_t>(sets.size()));
        if (fresh) {
            idx.nodeset_nodes.insert(idx.nodeset_nodes.end(), a.nodes.begin(), a.nodes.end());
            idx.nodeset_first.push_back(static_cast<std::uint32_t>(idx.nodeset_nodes.size()));
        }
        idx.nodeset.push_back(it->second);
        idx.mult.push_back(static_cast<std::uint16_t>(std::min<std::uint32_t>(a.mult, 65535)));
        present.push_back(a.paths);
    }
}

// Window fractions per read length, from each path's k-mer labels.
void compute_ceilings(const Graph& g, LocusIndex& idx) {
    std::vector<std::vector<double>> p3(std::size(kLengths)), r2(std::size(kLengths)), s2(std::size(kLengths));
    for (const auto& p : g.paths) {
        const std::string seq = spell_path_steps_sequence(g, p.steps);
        if (seq.size() < kK) continue;
        const std::size_t nk = seq.size() - kK + 1;
        std::vector<std::int8_t> lab(nk, -1);
        std::vector<char> single(nk, 0);
        kmer::for_each_canonical_kmer(seq.data(), seq.size(), kK, [&](std::uint64_t c, std::size_t start, bool) {
            const std::int64_t i = idx.find(c);
            if (i < 0) return;
            lab[start] = static_cast<std::int8_t>(idx.cls[i]);
            single[start] = idx.mult[i] == 1;
        });
        std::vector<std::uint32_t> pc(nk + 1, 0);
        for (std::size_t i = 0; i < nk; ++i) pc[i + 1] = pc[i] + (lab[i] == kPrivate);
        // runs of private (and of private single-copy) k-mers: start, end (inclusive); a window's sites are
        // ceil(clipped length / k) summed over the runs it meets
        auto runs_of = [&](bool single_only) {
            std::vector<std::pair<std::size_t, std::size_t>> r;
            for (std::size_t i = 0; i < nk;) {
                const auto ok = [&](std::size_t j) { return lab[j] == kPrivate && (!single_only || single[j]); };
                if (!ok(i)) { ++i; continue; }
                std::size_t e = i;
                while (e + 1 < nk && ok(e + 1)) ++e;
                r.emplace_back(i, e);
                i = e + 1;
            }
            return r;
        };
        const auto rp = runs_of(false), rs = runs_of(true);
        auto prefix = [&](const std::vector<std::pair<std::size_t, std::size_t>>& r) {
            std::vector<std::uint64_t> p(r.size() + 1, 0);
            for (std::size_t i = 0; i < r.size(); ++i) p[i + 1] = p[i] + (r[i].second - r[i].first + kK) / kK;
            return p;
        };
        const auto pp = prefix(rp), pq = prefix(rs);
        auto sites = [&](const std::vector<std::pair<std::size_t, std::size_t>>& r, const std::vector<std::uint64_t>& p, std::size_t a, std::size_t e) {
            // runs ending at or after a, starting at or before e
            const std::size_t lo = static_cast<std::size_t>(std::lower_bound(r.begin(), r.end(), a, [](const auto& x, std::size_t v) { return x.second < v; }) - r.begin());
            const std::size_t hi = static_cast<std::size_t>(std::upper_bound(r.begin(), r.end(), e, [](std::size_t v, const auto& x) { return v < x.first; }) - r.begin());
            if (lo >= hi) return std::uint64_t{0};
            std::uint64_t n = p[hi] - p[lo];
            auto full = [&](std::size_t i) { return (r[i].second - r[i].first + kK) / kK; };
            auto clip = [&](std::size_t i) { const std::size_t s0 = std::max(r[i].first, a), e0 = std::min(r[i].second, e); return (e0 - s0 + kK) / kK; };
            n = n - full(lo) + clip(lo);
            if (hi - 1 != lo) n = n - full(hi - 1) + clip(hi - 1);
            return n;
        };
        for (std::size_t li = 0; li < std::size(kLengths); ++li) {
            const std::size_t L = kLengths[li];
            if (seq.size() < L) continue;
            std::size_t n = 0, a = 0, b = 0, c = 0;
            for (std::size_t w = 0; w + L <= seq.size(); w += kWindowStep) {
                const std::size_t e = w + L - kK;   // last k-mer start in the window
                ++n;
                a += pc[e + 1] - pc[w] >= 3;
                b += sites(rp, pp, w, e) >= 2;
                c += sites(rs, pq, w, e) >= 2;
            }
            if (n == 0) continue;
            p3[li].push_back(static_cast<double>(a) / n);
            r2[li].push_back(static_cast<double>(b) / n);
            s2[li].push_back(static_cast<double>(c) / n);
        }
    }
    auto minmed = [](std::vector<double> v, double& mn, double& med) {
        if (v.empty()) { mn = med = 0; return; }
        std::sort(v.begin(), v.end());
        mn = v.front();
        med = v[v.size() / 2];
    };
    for (std::size_t li = 0; li < std::size(kLengths); ++li) {
        if (p3[li].empty()) continue;
        Ceiling c;
        c.length = kLengths[li];
        minmed(p3[li], c.p3_min, c.p3_median);
        minmed(r2[li], c.sites2_min, c.sites2_median);
        minmed(s2[li], c.single2_min, c.single2_median);
        idx.ceilings.push_back(c);
    }
}

// Binary writing and reading of the .plx file.
struct Out {
    std::ofstream f;
    void u64(std::uint64_t v) { f.write(reinterpret_cast<const char*>(&v), 8); }
    void str(const std::string& s) { u64(s.size()); f.write(s.data(), static_cast<std::streamsize>(s.size())); }
    template <class T> void vec(const std::vector<T>& v) { u64(v.size()); f.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(T))); }
};
struct In {
    std::ifstream f;
    std::string path;
    std::uint64_t u64() { std::uint64_t v = 0; if (!f.read(reinterpret_cast<char*>(&v), 8)) throw std::runtime_error("recruit: " + path + " is truncated"); return v; }
    std::string str() { std::string s(u64(), '\0'); if (!s.empty() && !f.read(&s[0], static_cast<std::streamsize>(s.size()))) throw std::runtime_error("recruit: " + path + " is truncated"); return s; }
    template <class T> std::vector<T> vec() {
        std::vector<T> v(u64());
        if (!v.empty() && !f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(T)))) throw std::runtime_error("recruit: " + path + " is truncated");
        return v;
    }
};

}  // namespace

std::int64_t LocusIndex::find(std::uint64_t canonical) const {
    const auto it = std::lower_bound(codes.begin(), codes.end(), canonical);
    return it != codes.end() && *it == canonical ? it - codes.begin() : -1;
}

std::string LocusIndex::key() const {
    std::string d;
    for (const auto& g : decoy_genomes) d += g + ",";
    return "v" + std::to_string(kVersion) + "|" + graph_checksum + "|" + reference_path + "|" + store_id + "|" + std::to_string(content) + "|" +
           coordinate_reference + "|" + std::to_string(informative_max) + "|" + std::to_string(pad) + "|" + std::to_string(dust) + "|" + std::to_string(decoy_min_hits) + "|" + d;
}

std::string expected_key(const BuildOptions& o, const std::string& gfa) {
    const bg::Manifest m = bg::read_manifest(o.background);
    LocusIndex k;
    k.graph_checksum = file_checksum(gfa);
    const Graph g = parse_gfa(gfa, ParseGfaOptions{true, false});
    std::vector<std::string> names;
    for (const auto& p : g.paths) names.push_back(p.name);
    const auto hit = match_reference_name(names, o.reference_query);
    if (hit.size() != 1) return {};
    k.reference_path = names[hit[0]];
    k.store_id = m.store_id;
    k.content = bg::content_id(m);
    k.coordinate_reference = o.coordinate_reference.empty() && m.references.size() == 1 ? m.references[0].first : o.coordinate_reference;
    k.informative_max = o.informative_max;
    k.pad = o.pad;
    k.dust = o.dust;
    k.decoy_min_hits = o.decoy_min_hits;
    k.decoy_genomes = o.decoy_genomes.empty() ? std::vector<std::string>{k.coordinate_reference} : o.decoy_genomes;
    return k.key();
}

std::vector<LocusIndex> build(const BuildOptions& o) {
    auto log = [&](const std::string& s) { if (o.log) o.log(s); };
    const bg::Manifest m = bg::read_manifest(o.background);
    std::string cref = o.coordinate_reference;
    if (cref.empty()) {
        if (m.references.size() != 1)
            throw std::runtime_error("recruit: give --coordinate-reference (the background has " + std::to_string(m.references.size()) + ")");
        cref = m.references[0].first;
    }
    if (std::none_of(m.references.begin(), m.references.end(), [&](const auto& r) { return r.first == cref; }))
        throw std::runtime_error("recruit: " + cref + " is not a coordinate reference of the background");
    std::map<std::string, bg::Genome> genomes;
    for (const auto& s : m.segments) for (const auto& g : s.genomes) genomes[g.name] = g;
    const std::vector<std::string> decoys = o.decoy_genomes.empty() ? std::vector<std::string>{cref} : o.decoy_genomes;
    for (const auto& d : decoys)
        if (!genomes.count(d)) throw std::runtime_error("recruit: decoy genome " + d + " is not in the background");
    std::string kind;
    const auto counts = open_counts(o.background, kind, o.log);
    log("background: " + counts->kind() + ", coordinate reference " + cref);
    const anchors::ReferenceAnchors ref = anchors::load_reference((fs::path(o.background) / "anchors" / cref / "reference.anc").string());
    std::vector<std::unique_ptr<anchors::SegmentAnchors>> seg_anchors;
    std::vector<const anchors::SegmentAnchors*> sa;
    for (const auto& s : m.segments) {
        seg_anchors.push_back(std::make_unique<anchors::SegmentAnchors>((fs::path(o.background) / s.dir / ("anchors-" + cref + ".bin")).string()));
        sa.push_back(seg_anchors.back().get());
    }

    // 1-3: each locus's k-mers, their background counts, their copies in every genome
    std::vector<LocusWork> work(o.loci.size());
    std::vector<std::unique_ptr<Graph>> graphs;
    for (std::size_t li = 0; li < o.loci.size(); ++li) {
        LocusWork& w = work[li];
        LocusIndex& x = w.idx;
        x.graph_path = fs::absolute(o.loci[li]).string();
        x.name = locus_name_of(o.loci[li]);
        x.graph_checksum = file_checksum(o.loci[li]);
        graphs.push_back(std::make_unique<Graph>(parse_gfa(o.loci[li])));
        const Graph& g = *graphs.back();
        std::vector<std::string> names;
        for (const auto& p : g.paths) names.push_back(p.name);
        const auto hit = match_reference_name(names, o.reference_query);
        if (hit.size() != 1)
            throw std::runtime_error("recruit: " + o.loci[li] + ": -r " + o.reference_query + " names " + std::to_string(hit.size()) + " paths, not one");
        x.reference_path = names[hit[0]];
        const std::string ref_seq = spell_path_steps_sequence(g, g.paths[hit[0]].steps);
        const ParsedReferencePath rp = parse_reference_path_label(x.reference_path, ref_seq.size(), g.paths[hit[0]].type == 'W');
        if (!rp.has_interval) throw std::runtime_error("recruit: " + x.reference_path + " carries no interval (contig:start-end)");
        x.contig = rp.chrom;
        x.start = rp.region_start_1based - 1;
        x.end = x.start + ref_seq.size();
        x.paths = g.paths.size();
        x.background = fs::absolute(o.background).string();
        x.store_id = m.store_id;
        x.content = bg::content_id(m);
        x.background_kind = kind;
        x.coordinate_reference = cref;
        x.informative_max = o.informative_max;
        x.pad = o.pad;
        x.dust = o.dust;
        x.decoy_min_hits = o.decoy_min_hits;
        x.decoy_genomes = decoys;
        std::vector<std::uint32_t> present;
        collect_kmers(g, x, present);
        x.backbone.resize(x.codes.size());
        for (std::size_t i = 0; i < x.codes.size(); ++i) x.backbone[i] = present[i] == x.paths && x.mult[i] == 1;
        w.counts = counts->count_batch(x.codes);
        w.sub.assign(x.codes.size(), 0);
        w.sub_trim.assign(x.codes.size(), 0);
        for (const auto& sp : anchors::locate(ref, sa, x.contig, x.start, x.end)) w.spans[sp.genome].push_back(sp);
        anchors::Span self;   // the coordinate reference: its exact interval
        self.genome = cref; self.contig = x.contig; self.start = x.start; self.end = x.end; self.inner_start = x.start; self.inner_end = x.end;
        w.spans[cref] = {self};
        log(x.name + ": " + std::to_string(x.codes.size()) + " k-mers from " + std::to_string(x.paths) + " paths; " + x.contig + ":" +
            std::to_string(x.start) + "-" + std::to_string(x.end) + " located in " + std::to_string(w.spans.size()) + " genomes");
    }

    // subtraction: read every span from its genome (one file at a time, every locus at once)
    for (const auto& [gname, ge] : genomes) {
        std::vector<fasta::Interval> iv;
        std::vector<std::pair<std::size_t, std::size_t>> owner;   // (locus, span)
        for (std::size_t li = 0; li < work.size(); ++li) {
            const auto it = work[li].spans.find(gname);
            if (it == work[li].spans.end()) continue;
            for (std::size_t s = 0; s < it->second.size(); ++s) {
                iv.push_back({it->second[s].contig, it->second[s].start, it->second[s].end});
                owner.emplace_back(li, s);
            }
        }
        if (iv.empty()) continue;
        std::vector<std::string> seqs;
        if (is_gfa_file(ge.path)) {
            const bg::graph::Haplotypes h(ge.path);
            for (const auto& i : iv) seqs.push_back(h.fetch(i.contig, i.start, i.end));
        } else {
            if (!fasta::IndexedFasta::indexed(ge.path)) log("note: " + ge.path + " has no .fai (and .gzi): read in full to cut its spans");
            seqs = fasta::fetch_intervals(ge.path, iv);
        }
        for (std::size_t k = 0; k < iv.size(); ++k) {
            LocusWork& w = work[owner[k].first];
            const anchors::Span& sp = w.spans[gname][owner[k].second];
            std::uint64_t found = 0;
            kmer::for_each_canonical_kmer(seqs[k].data(), seqs[k].size(), kK, [&](std::uint64_t c, std::size_t at, bool) {
                const std::int64_t i = w.idx.find(c);
                if (i < 0) return;
                ++w.sub[i];
                ++found;
                const std::uint64_t pos = iv[k].start + at;
                if (pos >= sp.inner_start && pos + kK <= sp.inner_end) ++w.sub_trim[i];
            });
            w.idx.spans.push_back({gname, sp.contig, sp.start, sp.end, sp.forward, found});
        }
    }

    // classes
    for (LocusWork& w : work) {
        LocusIndex& x = w.idx;
        const std::size_t n = x.codes.size();
        x.cls.resize(n);
        x.outside.resize(n);
        auto class_of = [&](std::uint64_t out) { return out == 0 ? kPrivate : out <= x.informative_max ? kShared : kRepeat; };
        for (std::size_t i = 0; i < n; ++i) {
            if (w.sub[i] > w.counts[i]) { ++x.negative; w.sub[i] = w.counts[i]; }
            const std::uint64_t out = w.counts[i] - w.sub[i];
            const std::uint64_t out_trim = w.counts[i] - std::min(w.sub_trim[i], w.counts[i]);
            x.outside[i] = static_cast<std::uint32_t>(std::min<std::uint64_t>(out, 0xFFFFFFFFULL));
            std::uint8_t c = class_of(out);
            if (c != class_of(out_trim) && c == kPrivate) { c = kShared; ++x.span_demoted; }
            if (c == kPrivate && kmer::dust_score(x.codes[i], kK) >= x.dust) { c = kShared; ++x.dust_demoted; }
            x.cls[i] = c;
        }
        if (x.negative > 0)
            log("warning: " + x.name + ": " + std::to_string(x.negative) + " k-mers occur more often in the located copies than in the background (a FASTA other than the one counted?)");
    }

    // decoys: the homologous windows of the decoy genomes around their S k-mers, outside the locus
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> shared;   // S k-mer -> loci
    for (std::size_t li = 0; li < work.size(); ++li)
        for (std::size_t i = 0; i < work[li].idx.codes.size(); ++i)
            if (work[li].idx.cls[i] == kShared) shared[work[li].idx.codes[i]].push_back(static_cast<std::uint32_t>(li));
    std::vector<std::unordered_set<std::uint64_t>> candidates(work.size());
    for (const auto& dname : decoys) {
        const bg::Genome& ge = genomes.at(dname);
        if (is_gfa_file(ge.path)) throw std::runtime_error("recruit: decoy genome " + dname + " is a graph; decoys must be FASTA");
        GzLineReader in(ge.path);
        if (!in.ok()) throw std::runtime_error("recruit: cannot read " + ge.path);
        std::string line, name, seq;
        auto scan = [&] {
            if (seq.empty()) return;
            std::vector<std::vector<std::uint64_t>> hits(work.size());
            kmer::for_each_canonical_kmer(seq.data(), seq.size(), kK, [&](std::uint64_t c, std::size_t at, bool) {
                const auto it = shared.find(c);
                if (it == shared.end()) return;
                for (const std::uint32_t li : it->second) {
                    bool own = false;   // inside this locus's own copy in this genome
                    const auto sp = work[li].spans.find(dname);
                    if (sp != work[li].spans.end())
                        for (const auto& s : sp->second) own = own || (s.contig == name && at + kK > s.start && at < s.end);
                    if (!own) hits[li].push_back(at);
                }
            });
            for (std::size_t li = 0; li < work.size(); ++li) {
                if (hits[li].empty()) continue;
                // clusters of hits no more than a pad apart; a cluster with enough hits is a homologous stretch,
                // its window padded on both sides
                const std::uint64_t pad = work[li].idx.pad;
                std::uint64_t first = 0, last = 0, count = 0;
                auto flush = [&] {
                    if (count < work[li].idx.decoy_min_hits) return;
                    std::uint64_t ws = first > pad ? first - pad : 0, we = std::min<std::uint64_t>(seq.size(), last + kK + pad);
                    // the padding stops at the locus's own copy in this genome: that sequence is the locus
                    const auto own = work[li].spans.find(dname);
                    if (own != work[li].spans.end())
                        for (const auto& sp : own->second) {
                            if (sp.contig != name || sp.end <= ws || sp.start >= we) continue;
                            if (last + kK <= sp.start) we = std::min<std::uint64_t>(we, sp.start);
                            else if (first >= sp.end) ws = std::max<std::uint64_t>(ws, sp.end);
                        }
                    if (we <= ws) return;
                    ++work[li].idx.decoy_windows;
                    work[li].idx.decoy_bases += we - ws;
                    work[li].idx.decoy_window_list.push_back({dname, name, ws, we, count});
                    kmer::for_each_canonical_kmer(seq.data() + ws, we - ws, kK, [&](std::uint64_t c, std::size_t, bool) {
                        if (work[li].idx.find(c) < 0) candidates[li].insert(c);
                    });
                };
                // the locus's own copies in this genome split clusters: hits on its two sides are two clusters
                std::vector<std::pair<std::uint64_t, std::uint64_t>> own_here;
                const auto own_it = work[li].spans.find(dname);
                if (own_it != work[li].spans.end())
                    for (const auto& sp : own_it->second) if (sp.contig == name) own_here.emplace_back(sp.start, sp.end);
                auto across_own = [&](std::uint64_t a, std::uint64_t b) {   // a < b
                    for (const auto& [os, oe] : own_here) if (a < os && b >= oe) return true;
                    return false;
                };
                for (const std::uint64_t at : hits[li]) {
                    if (count > 0 && at <= last + pad && !across_own(last, at)) { last = at; ++count; continue; }
                    flush();
                    first = last = at;
                    count = 1;
                }
                flush();
            }
        };
        while (in.getline(line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                scan();
                const auto sp = line.find_first_of(" \t");
                name = line.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
                seq.clear();
                continue;
            }
            seq += line;
        }
        scan();
        log("decoys: scanned " + dname);
    }
    for (std::size_t li = 0; li < work.size(); ++li) {
        LocusIndex& x = work[li].idx;
        std::vector<std::uint64_t> cand(candidates[li].begin(), candidates[li].end());
        std::sort(cand.begin(), cand.end());
        const auto cc = counts->count_batch(cand);
        for (std::size_t i = 0; i < cand.size(); ++i)
            if (cc[i] <= x.informative_max) x.decoy.push_back(cand[i]);
        // screen: closed syncmers among P and S
        for (std::size_t i = 0; i < x.codes.size(); ++i)
            if (x.cls[i] != kRepeat && kmer::is_closed_syncmer(x.codes[i], kK, kS)) x.screen.push_back(x.codes[i]);
        compute_ceilings(*graphs[li], x);
        std::size_t np = 0, ns = 0, nr = 0;
        for (const auto c : x.cls) (c == kPrivate ? np : c == kShared ? ns : nr)++;
        log(x.name + ": P " + std::to_string(np) + ", S " + std::to_string(ns) + ", R " + std::to_string(nr) + "; decoys " +
            std::to_string(x.decoy.size()) + " from " + std::to_string(x.decoy_windows) + " windows (" + std::to_string(x.decoy_bases) +
            " bp), screen " + std::to_string(x.screen.size()));
    }
    std::vector<LocusIndex> out;
    for (auto& w : work) out.push_back(std::move(w.idx));
    return out;
}

void save(const LocusIndex& x, const std::string& path) {
    const std::string tmp = path + ".tmp";
    {
        Out o{std::ofstream(tmp, std::ios::binary)};
        if (!o.f) throw std::runtime_error("recruit: cannot write " + tmp);
        o.u64(0);
        o.u64(kVersion);
        o.str(x.key());
        for (const auto* s : {&x.name, &x.graph_path, &x.graph_checksum, &x.reference_path, &x.contig, &x.background, &x.store_id,
                              &x.background_kind, &x.coordinate_reference})
            o.str(*s);
        o.u64(x.start); o.u64(x.end); o.u64(x.content); o.u64(x.informative_max); o.u64(x.pad);
        std::uint64_t dbits;
        std::memcpy(&dbits, &x.dust, 8);
        o.u64(dbits);
        o.u64(x.decoy_genomes.size());
        for (const auto& d : x.decoy_genomes) o.str(d);
        o.u64(x.paths);
        o.vec(x.codes); o.vec(x.cls); o.vec(x.outside); o.vec(x.mult); o.vec(x.backbone); o.vec(x.nodeset);
        o.vec(x.nodeset_first); o.vec(x.nodeset_nodes);
        o.u64(x.nodes.size());
        for (const auto& n : x.nodes) o.str(n);
        o.vec(x.decoy); o.vec(x.screen);
        o.u64(x.spans.size());
        for (const auto& s : x.spans) { o.str(s.genome); o.str(s.contig); o.u64(s.start); o.u64(s.end); o.u64(s.forward); o.u64(s.subtracted); }
        o.vec(x.ceilings);
        o.u64(x.dust_demoted); o.u64(x.span_demoted); o.u64(x.negative);
        o.u64(x.decoy_min_hits); o.u64(x.decoy_windows); o.u64(x.decoy_bases);
        o.u64(x.decoy_window_list.size());
        for (const auto& w : x.decoy_window_list) { o.str(w.genome); o.str(w.contig); o.u64(w.start); o.u64(w.end); o.u64(w.hits); }
        o.f.seekp(0);
        o.u64(kMagic);
        o.f.flush();
        if (!o.f) throw std::runtime_error("recruit: cannot write " + tmp);
    }
    fs::rename(tmp, path);
}

LocusIndex load(const std::string& path) {
    In i{std::ifstream(path, std::ios::binary), path};
    if (!i.f) throw std::runtime_error("recruit: cannot read " + path);
    if (i.u64() != kMagic) throw std::runtime_error("recruit: " + path + " is not a locus index of this version (or was interrupted)");
    if (i.u64() != kVersion) throw std::runtime_error("recruit: " + path + ": unknown version");
    LocusIndex x;
    const std::string key = i.str();
    for (auto* s : {&x.name, &x.graph_path, &x.graph_checksum, &x.reference_path, &x.contig, &x.background, &x.store_id,
                    &x.background_kind, &x.coordinate_reference})
        *s = i.str();
    x.start = i.u64(); x.end = i.u64(); x.content = i.u64();
    x.informative_max = static_cast<std::uint32_t>(i.u64());
    x.pad = i.u64();
    const std::uint64_t dbits = i.u64();
    std::memcpy(&x.dust, &dbits, 8);
    x.decoy_genomes.resize(i.u64());
    for (auto& d : x.decoy_genomes) d = i.str();
    x.paths = i.u64();
    x.codes = i.vec<std::uint64_t>(); x.cls = i.vec<std::uint8_t>(); x.outside = i.vec<std::uint32_t>(); x.mult = i.vec<std::uint16_t>();
    x.backbone = i.vec<std::uint8_t>(); x.nodeset = i.vec<std::uint32_t>();
    x.nodeset_first = i.vec<std::uint32_t>(); x.nodeset_nodes = i.vec<std::uint32_t>();
    x.nodes.resize(i.u64());
    for (auto& n : x.nodes) n = i.str();
    x.decoy = i.vec<std::uint64_t>(); x.screen = i.vec<std::uint64_t>();
    x.spans.resize(i.u64());
    for (auto& s : x.spans) { s.genome = i.str(); s.contig = i.str(); s.start = i.u64(); s.end = i.u64(); s.forward = i.u64() != 0; s.subtracted = i.u64(); }
    x.ceilings = i.vec<Ceiling>();
    x.dust_demoted = i.u64(); x.span_demoted = i.u64(); x.negative = i.u64();
    x.decoy_min_hits = static_cast<std::uint32_t>(i.u64()); x.decoy_windows = i.u64(); x.decoy_bases = i.u64();
    x.decoy_window_list.resize(i.u64());
    for (auto& w : x.decoy_window_list) { w.genome = i.str(); w.contig = i.str(); w.start = i.u64(); w.end = i.u64(); w.hits = i.u64(); }
    const std::size_t n = x.codes.size();
    if (x.cls.size() != n || x.outside.size() != n || x.mult.size() != n || x.backbone.size() != n || x.nodeset.size() != n ||
        key != x.key())
        throw std::runtime_error("recruit: " + path + " is corrupt");
    return x;
}

void write_report(const LocusIndex& x, const std::string& path) {
    std::ofstream r(path);
    std::size_t np = 0, ns = 0, nr = 0, bb = 0, multi = 0;
    for (std::size_t i = 0; i < x.codes.size(); ++i) {
        (x.cls[i] == kPrivate ? np : x.cls[i] == kShared ? ns : nr)++;
        bb += x.backbone[i];
        multi += x.mult[i] > 1;
    }
    auto pct = [&](std::size_t v) { char b[16]; std::snprintf(b, sizeof b, "%.1f", 100.0 * static_cast<double>(v) / static_cast<double>(std::max<std::size_t>(1, x.codes.size()))); return std::string(b); };
    r << "section\tkey\tvalue\n";
    r << "locus\tname\t" << x.name << "\nlocus\tgraph\t" << x.graph_path << "\nlocus\treference_path\t" << x.reference_path
      << "\nlocus\tinterval\t" << x.contig << ':' << x.start << '-' << x.end << "\nlocus\tpaths\t" << x.paths
      << "\nlocus\tbackground\t" << x.background << " (" << x.background_kind << ", content " << std::hex << x.content << std::dec << ")"
      << "\nlocus\tkmers\t" << x.codes.size() << "\nlocus\tprivate\t" << np << " (" << pct(np) << "%)"
      << "\nlocus\tshared\t" << ns << " (" << pct(ns) << "%)" << "\nlocus\trepeat\t" << nr << " (" << pct(nr) << "%)"
      << "\nlocus\tbackbone\t" << bb << "\nlocus\tmulti_copy\t" << multi << "\nlocus\tnode_sets\t" << (x.nodeset_first.empty() ? 0 : x.nodeset_first.size() - 1)
      << "\nlocus\tdust_demoted\t" << x.dust_demoted << "\nlocus\tspan_demoted\t" << x.span_demoted << "\nlocus\tnegative_residuals\t" << x.negative
      << "\nlocus\tdecoys\t" << x.decoy.size() << "\nlocus\tdecoy_windows\t" << x.decoy_windows << " (" << x.decoy_bases << " bp)"
      << "\nlocus\tscreen\t" << x.screen.size() << '\n';
    r << "span\tgenome\tcontig\tstart\tend\tstrand\tsubtracted\n";
    for (const auto& s : x.spans)
        r << "span\t" << s.genome << '\t' << s.contig << '\t' << s.start << '\t' << s.end << '\t' << (s.forward ? '+' : '-') << '\t' << s.subtracted << '\n';
    r << "decoy_window\tgenome\tcontig\tstart\tend\thits\n";
    for (const auto& w : x.decoy_window_list)
        r << "decoy_window\t" << w.genome << '\t' << w.contig << '\t' << w.start << '\t' << w.end << '\t' << w.hits << '\n';
    r << "ceiling\tread_length\tp3_min\tp3_median\tsites2_min\tsites2_median\tsingle2_min\tsingle2_median\n";
    for (const auto& c : x.ceilings)
        r << "ceiling\t" << c.length << '\t' << c.p3_min << '\t' << c.p3_median << '\t' << c.sites2_min << '\t' << c.sites2_median << '\t'
          << c.single2_min << '\t' << c.single2_median << '\n';
    if (!r) throw std::runtime_error("recruit: cannot write " + path);
}

}  // namespace panvar::recruit
