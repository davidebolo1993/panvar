#include "panvar/locus_index.hpp"

#include "panvar/anchors.hpp"
#include "panvar/background.hpp"
#include "panvar/background_graph.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/fasta_index.hpp"
#include "panvar/flat_map.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/gz_reader.hpp"
#include "panvar/kmer.hpp"
#include "panvar/recruit_calibration.hpp"
#include "panvar/ref_path.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <tuple>
#include <unistd.h>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>

namespace panvar::recruit {

namespace {

namespace fs = std::filesystem;

constexpr std::size_t kK = 31, kS = 11;
constexpr std::uint64_t kMagic = 0x3130305850585650ULL;   // "PVPXP001"
constexpr std::uint64_t kVersion = 16;   // 2: decoy windows stop at the locus; 3: clusters split at it; 4: S per haplotype; 5: calibration; 6: its decoy pairs cross window edges; 7: locus pieces off the located copies; 8: and past estimated ends; 9: and everywhere for a locus with array pieces outside; 10: stretches reach contig ends; 11: pieces tied to the locus, short contigs whole, no decoys in pieces; 12: copies walked out from their anchors; 13: across repeats; 14: a checksum; 15: pieces only next to a copy or at a contig end; 16: the score's calibration, genomes held out
const std::size_t kLengths[] = {100, 150, 250, 1000, 10000};
constexpr std::size_t kWindowStep = 25;
// A stretch of a genome is a piece of the locus when, over every window of this many bases in it, this
// share of the k-mer starts are locus k-mers (a held-out haplotype's new alleles, about one per kb, leave
// about 97%) and it is tied to the locus (see the piece scan): a close paralog, or one whose differences
// cluster, can be as dense.
constexpr std::uint64_t kPieceWindow = 5000;
constexpr double kPieceDensity = 0.95;
// A stretch is next to a located copy, or to another piece, when at most this many bases lie between them:
// less than a window of other sequence, as where an array's density dips
constexpr std::uint64_t kPieceGap = kPieceWindow;
// So a genome indexed for reading by interval is read only where a piece can be: around each located copy and
// at each contig's ends, this many bases at first, widened while a piece comes within the margin of a part's
// edge (closer, a window or a neighbouring piece past the edge could still change it)
constexpr std::uint64_t kPiecePad = 8 * kPieceWindow;
constexpr std::uint64_t kPieceMargin = kPieceGap + 3 * kPieceWindow;
// The piece scan looks up every kPieceSample-th k-mer start first: a window of kPieceWindow starts holding a share
// kPieceDensity of the locus's has, at worst, all its misses on sampled starts, and still that many sampled hits
constexpr std::uint64_t kPieceSample = 8;
static_assert(kPieceWindow % kPieceSample == 0, "a window holds a whole number of sampled starts");
// A located copy is certain between its outermost inside anchors; past them its ends are estimated, and may
// overshoot into flank. From the anchored stretch the copy is walked out over windows of this many k-mer starts
// in which this share are the locus's: a copy's own sequence, new alleles and all, keeps it up, a flank holds
// almost none. Repeats count only to cross from one stretch of the locus to the next. Reference anchors cluster
// where a locus is like the reference, so at a locus whose haplotypes differ from it most of each copy lies past
// them.
constexpr std::uint64_t kWalkWindow = 1000;
constexpr double kWalkDensity = 0.5;
// A locus is looked for in every genome when this share of its k-mers are tandem (this many copies or more
// within a haplotype) and left with a few copies outside the located ones (at most this many times their own
// copy number): the mark of array pieces next to copies that look whole. A locus without an array has a few
// hundredths of a percent.
constexpr std::uint16_t kTandemCopies = 3;
constexpr std::uint64_t kTandemOutside = 10;
constexpr double kTandemShare = 0.005;

// A located copy walked out from its anchored stretch [inner_start, inner_end) through the stretches dense in the
// locus's k-mers (kWalkWindow, kWalkDensity), within [start, end). Repeats (above `repeat_above` background copies)
// do not carry the walk on their own: a stretch dense only with them joins the copy when locus sequence goes on past
// it (a repeat element inside the locus), not when it only runs on into the flank. `hits` are the locus's k-mer
// occurrences in the span (position, k-mer), sorted by position.
std::pair<std::uint64_t, std::uint64_t> walk_out(const std::vector<std::pair<std::uint64_t, std::uint32_t>>& hits,
                                                 const std::vector<std::uint64_t>& counts, std::uint64_t repeat_above,
                                                 std::uint64_t start, std::uint64_t end, std::uint64_t inner_start, std::uint64_t inner_end) {
    std::vector<std::uint64_t> p, all;   // positions of the locus's k-mers: repeats left out, and all of them
    for (const auto& [pos, i] : hits) {
        if (all.empty() || all.back() != pos) all.push_back(pos);
        if (counts[i] <= repeat_above && (p.empty() || p.back() != pos)) p.push_back(pos);
    }
    std::uint64_t lo = inner_start, hi = inner_end;
    if (end < start + kK || p.empty()) return {lo, hi};
    const std::uint64_t starts = end - start - kK + 1;
    if (starts < kWalkWindow) {   // no window fits: the whole span, when dense enough over all of it
        if (static_cast<double>(p.size()) >= kWalkDensity * static_cast<double>(starts)) return {start, end};
        return {lo, hi};
    }
    const std::size_t need = static_cast<std::size_t>(std::ceil(kWalkDensity * static_cast<double>(kWalkWindow)));
    auto dense_of = [&](const std::vector<std::uint64_t>& q) {   // the union of the windows holding enough of them
        std::vector<std::pair<std::uint64_t, std::uint64_t>> dense;
        for (std::size_t i = 0, j = 0; i < q.size(); ++i) {
            if (j < i) j = i;
            while (j < q.size() && q[j] < q[i] + kWalkWindow) ++j;
            if (j - i < need) continue;
            const std::uint64_t ws = q[i], we = q[j - 1] + kK;
            if (!dense.empty() && ws <= dense.back().second) dense.back().second = std::max(dense.back().second, we);
            else dense.emplace_back(ws, we);
        }
        return dense;
    };
    const auto dense = dense_of(p), bridges = dense_of(all);   // a stretch of `dense` lies within one of `bridges`
    auto touch = [](std::uint64_t a, std::uint64_t b, std::uint64_t c, std::uint64_t d) { return a <= d && b >= c; };
    for (bool grown = true; grown;) {   // stretches touching the copy, or reached from it across repeats, join it
        grown = false;
        for (const auto& [a, b] : dense) {
            if (a >= lo && b <= hi) continue;
            bool reach = touch(a, b, lo, hi);
            for (std::size_t k = 0; !reach && k < bridges.size(); ++k)
                reach = touch(bridges[k].first, bridges[k].second, lo, hi) && touch(bridges[k].first, bridges[k].second, a, b);
            if (reach) {
                lo = std::min(lo, a);
                hi = std::max(hi, b);
                grown = true;
            }
        }
    }
    return {std::max(lo, start), std::min(hi, end)};
}

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

// A locus is a graph (GFA) or its haplotypes (FASTA). A FASTA locus is read as a graph of one segment per record, each
// the only step of a path named as the record: the index needs only the paths' sequences, and a graph's paths spelled
// out give the same index but for the nodes each private k-mer starts on, which a FASTA locus does not have.
Graph read_locus(const std::string& path, bool sequences, bool& fasta) {
    fasta = !is_gfa_file(path);
    if (!fasta) return parse_gfa(path, ParseGfaOptions{true, sequences});
    GzLineReader in(path);
    if (!in.ok()) throw std::runtime_error("recruit: cannot read " + path);
    Graph g;
    std::unordered_set<std::string> names;
    std::string line;
    Node* node = nullptr;
    while (in.getline(line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') {
            const auto sp = line.find_first_of(" \t");
            const std::string name = line.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
            if (name.empty()) throw std::runtime_error("recruit: " + path + ": a FASTA record without a name");
            if (!names.insert(name).second) throw std::runtime_error("recruit: " + path + ": two records named " + name);
            const std::string id = std::to_string(g.paths.size() + 1);
            node = &g.nodes[id];
            node->id = id;
            g.paths.push_back(PathRecord{name, 'P', {PathStep{id, false}}});
            continue;
        }
        if (node == nullptr) throw std::runtime_error("recruit: " + path + " is neither a GFA graph nor a FASTA file");
        if (sequences) node->sequence += line;
    }
    if (g.paths.empty()) throw std::runtime_error("recruit: " + path + ": no locus graph or FASTA records");
    return g;
}

// Which counts to read: the compacted dictionary beside the store when it holds the manifest's genomes (`content`),
// else the store, which must hold them too: a commit between reading the manifest and the counts (another job's
// --add) would mix two sets of genomes in one index.
std::unique_ptr<bg::Background> open_counts(const std::string& store, std::uint64_t content, std::string& kind,
                                            const std::function<void(const std::string&)>& log) {
    std::string base = store;
    while (base.size() > 1 && base.back() == '/') base.pop_back();
    if (fs::is_regular_file(base + ".sshash")) {
        try {
            auto b = bg::Background::open(base + ".sshash");
            if (b->content() == content) {
                kind = "compacted";
                return b;
            }
            if (log) log("not using the compacted background: it was compacted from other genomes than the store holds now");
        } catch (const std::exception& e) {
            if (log) log(std::string("not using the compacted background: ") + e.what());
        }
    }
    kind = bg::is_package(store) ? "compacted" : "store";   // a package's counts are its compacted dictionary
    auto b = bg::Background::open(store);
    if (b->content() != content)
        throw std::runtime_error("recruit: the background " + store + " changed while it was read (a commit by another job): run again");
    return b;
}

struct LocusWork {
    LocusIndex idx;
    std::vector<std::uint64_t> counts;                 // background counts per k-mer
    std::vector<std::uint64_t> sub, sub_trim;          // own occurrences in the extended / trimmed spans
    std::uint64_t ref_kmers = 0;                       // the reference path's k-mers (its interval must hold them)
    std::string ref_seq;                               // the reference path's sequence (its interval must spell it)
    std::string error;                                 // why this locus failed after the shared scans (empty: fine)
    std::map<std::string, std::vector<anchors::Span>> spans;   // per genome: the located copies, then the pieces
                                                               // found off them (the locus, both, from there on)
    std::vector<std::uint64_t> decoy_counts;           // background counts of the decoy k-mers (parallel to idx.decoy)
};

// Every k-mer of every path: nodes it starts on (none for a FASTA locus), largest copy number within a path, paths
// holding it.
void collect_kmers(const Graph& g, LocusIndex& idx, std::vector<std::uint32_t>& present, bool with_nodes) {
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
    if (with_nodes) idx.nodes = ids;

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
                if (with_nodes && (a.nodes.empty() || std::find(a.nodes.begin(), a.nodes.end(), occ[e].second) == a.nodes.end()))
                    a.nodes.push_back(occ[e].second);
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
        // two sites in a window, counted as the scorer counts them (a site opens at a private k-mer and takes in
        // those starting less than k after it, whatever lies between): the window's last private k-mer starts k or
        // more after its first. Nearest private (and private single-copy) k-mer at or after / at or before each start.
        constexpr std::size_t kNone = static_cast<std::size_t>(-1);
        auto nearest = [&](bool single_only, std::vector<std::size_t>& next, std::vector<std::size_t>& prev) {
            const auto ok = [&](std::size_t j) { return lab[j] == kPrivate && (!single_only || single[j]); };
            next.assign(nk + 1, kNone);
            prev.assign(nk, kNone);
            for (std::size_t i = nk; i-- > 0;) next[i] = ok(i) ? i : next[i + 1];
            for (std::size_t i = 0; i < nk; ++i) prev[i] = ok(i) ? i : (i > 0 ? prev[i - 1] : kNone);
        };
        std::vector<std::size_t> next_p, prev_p, next_s, prev_s;
        nearest(false, next_p, prev_p);
        nearest(true, next_s, prev_s);
        auto two_sites = [&](const std::vector<std::size_t>& next, const std::vector<std::size_t>& prev, std::size_t a, std::size_t e) {
            const std::size_t f = next[a], l = prev[e];
            return f != kNone && f <= e && l != kNone && l >= f + kK;
        };
        for (std::size_t li = 0; li < std::size(kLengths); ++li) {
            const std::size_t L = kLengths[li];
            if (seq.size() < L) continue;
            std::size_t n = 0, a = 0, b = 0, c = 0;
            for (std::size_t w = 0; w + L <= seq.size(); w += kWindowStep) {
                const std::size_t e = w + L - kK;   // last k-mer start in the window
                ++n;
                a += pc[e + 1] - pc[w] >= 3;
                b += two_sites(next_p, prev_p, w, e);
                c += two_sites(next_s, prev_s, w, e);
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
// A running checksum over the bytes written or read after the magic word, stored last: a damaged index is
// refused when loaded rather than scoring with wrong labels.
std::uint64_t hash_bytes_into(std::uint64_t h, const void* p, std::size_t n) {
    const auto* b = static_cast<const unsigned char*>(p);
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) { std::uint64_t w; std::memcpy(&w, b + i, 8); h = kmer::mix64(h ^ w); }
    std::uint64_t w = 0;
    std::memcpy(&w, b + i, n - i);
    return kmer::mix64(h ^ w ^ (static_cast<std::uint64_t>(n) << 56));
}
constexpr std::uint64_t kChecksumSeed = 0x6a09e667f3bcc909ULL;
struct Out {
    std::ofstream f;
    std::uint64_t h = kChecksumSeed;
    void bytes(const void* p, std::size_t n) { f.write(static_cast<const char*>(p), static_cast<std::streamsize>(n)); h = hash_bytes_into(h, p, n); }
    void u64(std::uint64_t v) { bytes(&v, 8); }
    void raw64(std::uint64_t v) { f.write(reinterpret_cast<const char*>(&v), 8); }   // outside the checksum
    void f64(double v) { std::uint64_t b; std::memcpy(&b, &v, 8); u64(b); }
    void str(const std::string& s) { u64(s.size()); bytes(s.data(), s.size()); }
    template <class T> void vec(const std::vector<T>& v) { u64(v.size()); bytes(v.data(), v.size() * sizeof(T)); }
};
struct In {
    std::ifstream f;
    std::string path;
    std::uint64_t h = kChecksumSeed;
    void bytes(void* p, std::size_t n) {
        if (n > 0 && !f.read(static_cast<char*>(p), static_cast<std::streamsize>(n))) throw std::runtime_error("recruit: " + path + " is truncated");
        h = hash_bytes_into(h, p, n);
    }
    std::uint64_t u64() { std::uint64_t v = 0; bytes(&v, 8); return v; }
    double f64() { const std::uint64_t b = u64(); double v; std::memcpy(&v, &b, 8); return v; }
    std::uint64_t raw64() { std::uint64_t v = 0; if (!f.read(reinterpret_cast<char*>(&v), 8)) throw std::runtime_error("recruit: " + path + " is truncated"); return v; }
    // a length read from the file, bounded by what is left of it, so a damaged length is refused, not allocated
    std::uint64_t length(std::size_t unit) {
        const std::uint64_t n = u64();
        const auto at = f.tellg();
        f.seekg(0, std::ios::end);
        const auto end = f.tellg();
        f.seekg(at);
        if (unit > 0 && n > static_cast<std::uint64_t>(end - at) / unit) throw std::runtime_error("recruit: " + path + " is corrupt");
        return n;
    }
    std::string str() { std::string s(length(1), '\0'); bytes(s.data(), s.size()); return s; }   // empty ones hashed too, as written
    template <class T> std::vector<T> vec() {
        std::vector<T> v(length(sizeof(T)));
        bytes(v.data(), v.size() * sizeof(T));
        return v;
    }
};

}  // namespace

std::int64_t LocusIndex::find(std::uint64_t canonical) const {
    const auto it = std::lower_bound(codes.begin(), codes.end(), canonical);
    return it != codes.end() && *it == canonical ? it - codes.begin() : -1;
}

std::string LocusIndex::key() const {
    std::string d, cg;
    for (const auto& g : decoy_genomes) d += g + ",";
    for (const auto& g : calibration_genomes_given) cg += g + ",";
    return "v" + std::to_string(kVersion) + "|" + graph_checksum + "|" + reference_path + "|" + store_id + "|" + std::to_string(content) + "|" +
           coordinate_reference + "|" + std::to_string(informative_max) + "|" + std::to_string(pad) + "|" + std::to_string(dust) + "|" + std::to_string(decoy_min_hits) + "|" + d +
           (calibrated ? "|c1|" + (cg.empty() ? std::string("auto") : cg) : "|c0");
}

std::string locus_name(const std::string& path) {
    std::string n = fs::path(path).filename().string();
    auto strip = [&](const char* ext) {
        const std::size_t e = std::strlen(ext);
        if (n.size() > e && n.compare(n.size() - e, e, ext) == 0) { n.resize(n.size() - e); return true; }
        return false;
    };
    strip(".gz");
    for (const char* ext : {".gfa", ".gfa1", ".fa", ".fasta", ".fna"})
        if (strip(ext)) break;
    return n;
}

std::string expected_key(const BuildOptions& o, const std::string& gfa) {
    const bg::Manifest m = bg::read_manifest(o.background);
    LocusIndex k;
    k.graph_checksum = file_checksum(gfa);
    bool fasta = false;
    const Graph g = read_locus(gfa, false, fasta);
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
    k.calibrated = o.calibrate;
    k.calibration_genomes_given = o.calibration_genomes;
    return k.key();
}

// A genome the index reads again (to cut its spans, to scan for decoys) must be the file the background
// counted: a different one would subtract or place k-mers the counts never held. Its size is checked (its
// sequence checksum would mean reading every genome whole).
void check_unchanged(const bg::Genome& g) {
    std::error_code ec;
    const std::uint64_t now = fs::file_size(g.path, ec);
    if (ec) throw std::runtime_error("recruit: genome " + g.name + ": cannot read " + g.path + " (the background counted it there)");
    if (g.bytes != 0 && now != g.bytes)
        throw std::runtime_error("recruit: genome " + g.name + ": " + g.path + " is " + std::to_string(now) + " bytes, " +
                                 std::to_string(g.bytes) + " when the background counted it: the file was replaced");
}

// Calls f(i) for i in [0, n) on up to `threads` threads, each i on one thread; the first exception thrown
// is rethrown once all have stopped.
void parallel_for(std::size_t n, std::size_t threads, const std::function<void(std::size_t)>& f) {
    std::atomic<std::size_t> next{0};
    std::exception_ptr failed;
    std::mutex m;
    auto run = [&] {
        for (std::size_t i; (i = next++) < n;) {
            try {
                f(i);
            } catch (...) {
                std::lock_guard<std::mutex> lock(m);
                if (!failed) failed = std::current_exception();
                next = n;
            }
        }
    };
    std::vector<std::thread> pool;
    for (std::size_t t = 1; t < std::min(threads, n); ++t) pool.emplace_back(run);
    run();
    for (auto& t : pool) t.join();
    if (failed) std::rethrow_exception(failed);
}

// A calibration genome's off-target sequence for a locus: its homologous windows with context (cut where the
// locus's own copy begins, as the windows are) and kCalibrationGenomeFlank bases on each side of its copies,
// none of the locus's copies or pieces; read as merged blocks, so no base is counted twice.
struct GenomeRegion {
    struct Rec {
        std::size_t block = 0;
        std::uint64_t start = 0, end = 0, core_start = 0, core_end = 0;   // on the contig
        bool window = false;
    };
    std::vector<fasta::Interval> blocks;   // by contig and start, disjoint
    std::vector<Rec> recs;
    std::uint64_t windows = 0;
};

GenomeRegion region_of(const LocusWork& w, const std::string& genome) {
    const LocusIndex& x = w.idx;
    GenomeRegion r;
    const auto own_it = w.spans.find(genome);
    static const std::vector<anchors::Span> kNone;
    const std::vector<anchors::Span>& own = own_it == w.spans.end() ? kNone : own_it->second;
    std::vector<GenomeRegion::Rec> recs;
    std::vector<std::string> contig_of;
    for (const auto* list : {&x.decoy_window_list, &x.calibration_window_list})
        for (const DecoyWindow& win : *list) {
            if (win.genome != genome) continue;
            std::uint64_t s = win.start > kCalibrationContext ? win.start - kCalibrationContext : 0, e = win.end + kCalibrationContext;
            for (const auto& sp : own) {
                if (sp.contig != win.contig || sp.end <= s || sp.start >= e) continue;
                if (sp.end <= win.start) s = std::max<std::uint64_t>(s, sp.end);
                else if (sp.start >= win.end) e = std::min<std::uint64_t>(e, sp.start);
            }
            if (e <= s) continue;
            ++r.windows;
            recs.push_back({0, s, e, win.start, win.end, true});
            contig_of.push_back(win.contig);
        }
    for (const GenomeSpan& sp : x.spans) {
        if (sp.genome != genome || sp.piece) continue;
        const std::pair<std::uint64_t, std::uint64_t> sides[2] = {{sp.start > kCalibrationGenomeFlank ? sp.start - kCalibrationGenomeFlank : 0, sp.start},
                                                                  {sp.end, sp.end + kCalibrationGenomeFlank}};
        for (const auto& [a, b] : sides) {
            // the side less every copy and piece of the locus on the contig
            std::vector<std::pair<std::uint64_t, std::uint64_t>> parts{{a, b}};
            for (const auto& o : own) {
                if (o.contig != sp.contig) continue;
                std::vector<std::pair<std::uint64_t, std::uint64_t>> next;
                for (const auto& [ps, pe] : parts) {
                    if (o.end <= ps || o.start >= pe) { next.emplace_back(ps, pe); continue; }
                    if (o.start > ps) next.emplace_back(ps, o.start);
                    if (o.end < pe) next.emplace_back(o.end, pe);
                }
                parts.swap(next);
            }
            for (const auto& [ps, pe] : parts)
                if (pe > ps) {
                    recs.push_back({0, ps, pe, ps, pe, false});
                    contig_of.push_back(sp.contig);
                }
        }
    }
    std::vector<std::size_t> order(recs.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return std::tie(contig_of[a], recs[a].start, recs[a].end) < std::tie(contig_of[b], recs[b].start, recs[b].end);
    });
    for (const std::size_t i : order) {
        GenomeRegion::Rec rec = recs[i];
        if (r.blocks.empty() || r.blocks.back().contig != contig_of[i] || rec.start > r.blocks.back().end)
            r.blocks.push_back({contig_of[i], rec.start, rec.end});
        else
            r.blocks.back().end = std::max(r.blocks.back().end, rec.end);
        rec.block = r.blocks.size() - 1;
        r.recs.push_back(rec);
    }
    return r;
}

// What a calibration genome holds for a locus, read from its blocks: the off-target records, and the k-mers that
// change label while it is held out.
struct GenomeRead {
    CalibrationGenome stats;
    GenomeOffTarget held;
    std::vector<std::uint32_t> homolog;   // its distinct homolog k-mers: locus (not private) indices, then decoys after them
};

GenomeRead read_region(const LocusWork& w, const FlatMap<std::uint32_t>& code_index, const GenomeRegion& r, const std::vector<std::string>& seqs,
                       bool keep_records) {
    const LocusIndex& x = w.idx;
    GenomeRead out;
    out.stats.windows = r.windows;
    std::vector<std::uint32_t> in_locus(x.codes.size(), 0), in_decoy(x.decoy.size(), 0);
    for (const std::string& b : seqs) {
        out.stats.bases += b.size();
        kmer::for_each_canonical_kmer(b.data(), b.size(), kK, [&](std::uint64_t c, std::size_t, bool) {
            if (const std::uint32_t* i = code_index.find(c)) { ++in_locus[*i]; return; }
            const auto it = std::lower_bound(x.decoy.begin(), x.decoy.end(), c);
            if (it != x.decoy.end() && *it == c) ++in_decoy[static_cast<std::size_t>(it - x.decoy.begin())];
        });
    }
    for (std::size_t i = 0; i < in_locus.size(); ++i) {
        if (in_locus[i] == 0 || x.cls[i] == kPrivate) continue;
        out.homolog.push_back(static_cast<std::uint32_t>(i));
        // every occurrence outside the locus is in this genome: without it, the k-mer would be private
        if (x.outside[i] == in_locus[i]) out.held.relabel.push_back(static_cast<std::uint32_t>(i));
    }
    for (std::size_t j = 0; j < in_decoy.size(); ++j) {
        if (in_decoy[j] == 0) continue;
        out.homolog.push_back(static_cast<std::uint32_t>(x.codes.size() + j));
        if (w.decoy_counts[j] == in_decoy[j]) out.held.unlabel.push_back(x.decoy[j]);
    }
    out.stats.homolog_kmers = out.homolog.size();
    out.stats.relabelled = out.held.relabel.size();
    out.stats.unlabelled = out.held.unlabel.size();
    if (!keep_records) return out;
    for (const GenomeRegion::Rec& rec : r.recs) {
        const fasta::Interval& b = r.blocks[rec.block];
        const std::string& seq = seqs[rec.block];
        const std::uint64_t a = rec.start - b.start, e = std::min<std::uint64_t>(rec.end - b.start, seq.size());
        if (e <= a) continue;   // past the contig's end
        OffTarget t;
        t.decoy = rec.window;
        t.seq = seq.substr(a, e - a);
        t.core_start = std::min<std::uint64_t>(rec.core_start - rec.start, t.seq.size());
        t.core_end = std::min<std::uint64_t>(rec.core_end - rec.start, t.seq.size());
        out.held.off.push_back(std::move(t));
    }
    return out;
}

// The genomes each locus holds out in calibration, and their off-target sequence. Given, every locus takes the list.
// By default each locus picks among the decoy genomes other than the coordinate reference, up to kCalibrationGenomesMax:
// first the one whose windows and flanks hold the most distinct homolog k-mers (locus k-mers that are not private,
// and decoy k-mers), then each next the one adding the most k-mers absent from those chosen (ties: the earlier in the
// list), until none adds any. Genomes are read one at a time, every locus at once; with more candidates than
// can be taken, twice (to choose, then to keep the chosen ones' sequence).
void read_calibration_genomes(std::vector<LocusWork>& work, const std::map<std::string, bg::Genome>& genomes, const std::string& cref,
                              const std::vector<std::string>& decoys, const std::vector<std::string>& given,
                              const std::vector<FlatMap<std::uint32_t>>& code_index, std::size_t threads,
                              const std::function<void(const std::string&)>& log, std::vector<std::vector<GenomeOffTarget>>& out) {
    std::vector<std::string> candidates = given;
    if (given.empty())
        for (const auto& d : decoys)
            if (d != cref) candidates.push_back(d);
    std::size_t most = kCalibrationGenomesMax;
    if (const char* m = std::getenv("PANVAR_TEST_CALIBRATION_GENOMES_MAX"); m != nullptr) most = std::max<std::size_t>(1, std::strtoull(m, nullptr, 10));
    const bool choose = given.empty() && candidates.size() > most;
    // per locus, the genomes it takes (indices into candidates)
    std::vector<std::vector<std::size_t>> chosen(work.size());
    auto read_genome = [&](std::size_t ci, const std::vector<std::size_t>& loci, bool keep, std::vector<GenomeRead>& got) {
        const std::string& gname = candidates[ci];
        const bg::Genome& ge = genomes.at(gname);
        check_unchanged(ge);
        std::vector<GenomeRegion> regions(loci.size());
        std::vector<fasta::Interval> all;
        for (std::size_t k = 0; k < loci.size(); ++k) {
            regions[k] = region_of(work[loci[k]], gname);
            all.insert(all.end(), regions[k].blocks.begin(), regions[k].blocks.end());
        }
        std::vector<std::string> seqs = all.empty() ? std::vector<std::string>() : fasta::fetch_intervals(ge.path, all);
        std::vector<std::vector<std::string>> per(loci.size());
        for (std::size_t k = 0, at = 0; k < loci.size(); ++k)
            for (std::size_t b = 0; b < regions[k].blocks.size(); ++b) per[k].push_back(std::move(seqs[at++]));
        got.assign(loci.size(), {});
        parallel_for(loci.size(), threads, [&](std::size_t k) {
            got[k] = read_region(work[loci[k]], code_index[loci[k]], regions[k], per[k], keep);
            got[k].stats.genome = gname;
        });
    };
    std::vector<std::size_t> every(work.size());
    for (std::size_t li = 0; li < work.size(); ++li) every[li] = li;
    if (choose) {
        std::vector<std::vector<std::vector<std::uint32_t>>> homolog(work.size(), std::vector<std::vector<std::uint32_t>>(candidates.size()));
        for (std::size_t ci = 0; ci < candidates.size(); ++ci) {
            std::vector<GenomeRead> got;
            read_genome(ci, every, false, got);
            for (std::size_t li = 0; li < work.size(); ++li) homolog[li][ci] = std::move(got[li].homolog);
            log("calibration genomes: candidate read: " + candidates[ci]);
        }
        for (std::size_t li = 0; li < work.size(); ++li) {
            const LocusIndex& x = work[li].idx;
            std::vector<char> covered(x.codes.size() + x.decoy.size(), 0), taken(candidates.size(), 0);
            while (chosen[li].size() < most) {
                std::size_t best = candidates.size(), best_gain = 0;
                for (std::size_t ci = 0; ci < candidates.size(); ++ci) {
                    if (taken[ci]) continue;
                    std::size_t gain = 0;
                    for (const std::uint32_t id : homolog[li][ci]) gain += !covered[id];
                    if (gain > best_gain) { best = ci; best_gain = gain; }
                }
                if (best == candidates.size()) break;   // none adds a homolog k-mer
                taken[best] = 1;
                chosen[li].push_back(best);
                for (const std::uint32_t id : homolog[li][best]) covered[id] = 1;
            }
        }
    } else {
        for (std::size_t li = 0; li < work.size(); ++li)
            for (std::size_t ci = 0; ci < candidates.size(); ++ci) chosen[li].push_back(ci);
    }
    // the chosen genomes' sequence, each genome once for every locus taking it
    std::vector<std::vector<std::pair<std::size_t, GenomeRead>>> found(work.size());   // (rank in chosen, what was read)
    for (std::size_t ci = 0; ci < candidates.size(); ++ci) {
        std::vector<std::size_t> loci;
        for (std::size_t li = 0; li < work.size(); ++li)
            if (work[li].error.empty() && std::find(chosen[li].begin(), chosen[li].end(), ci) != chosen[li].end()) loci.push_back(li);
        if (loci.empty()) continue;
        std::vector<GenomeRead> got;
        read_genome(ci, loci, true, got);
        for (std::size_t k = 0; k < loci.size(); ++k) {
            const std::size_t rank = static_cast<std::size_t>(std::find(chosen[loci[k]].begin(), chosen[loci[k]].end(), ci) - chosen[loci[k]].begin());
            found[loci[k]].emplace_back(rank, std::move(got[k]));
        }
        log("calibration genomes: read: " + candidates[ci] + " (for " + std::to_string(loci.size()) + " loci)");
    }
    for (std::size_t li = 0; li < work.size(); ++li) {
        LocusIndex& x = work[li].idx;
        std::sort(found[li].begin(), found[li].end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        x.calibration_genomes.clear();
        for (auto& [rank, gr] : found[li]) {
            gr.held.genome = x.calibration_genomes.size();
            x.calibration_genomes.push_back(gr.stats);
            out[li].push_back(std::move(gr.held));
        }
        std::string names;
        for (const auto& c : x.calibration_genomes) names += (names.empty() ? "" : ",") + c.genome;
        if (work[li].error.empty())
            log("locus: " + x.name + "; calibration genomes: " +
                (names.empty() ? std::string("none (no genome to hold out besides the coordinate reference: off-target reads gain private sites by sequencing errors only)")
                               : names));
    }
}

std::vector<LocusIndex> build(const BuildOptions& o) {
    std::mutex log_mutex;
    auto log = [&](const std::string& s) {
        if (!o.log) return;
        std::lock_guard<std::mutex> lock(log_mutex);
        o.log(s);
    };
    const std::size_t threads = std::max<std::size_t>(1, o.threads);
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
    // background haplotypes: a FASTA genome is one, a graph one per sample#haplotype of its P and W lines
    // (PanSN path names: sample#haplotype#contig, or sample#contig for a haploid sample); a graph whose
    // P lines are not PanSN names cannot say how many haplotypes it holds and counts as one
    std::uint64_t haplotypes = 0;
    for (const auto& [gname, ge] : genomes) {
        if (!fs::is_regular_file(ge.path)) throw std::runtime_error("recruit: genome " + gname + ": no file at " + ge.path);
        if (!is_gfa_file(ge.path)) { ++haplotypes; continue; }
        std::set<std::string> haps;
        bool plain = false;
        GzLineReader in(ge.path);
        if (!in.ok()) throw std::runtime_error("recruit: cannot read " + ge.path);
        for (std::string line; in.getline(line);) {
            if (line.size() < 2 || line[1] != '\t' || (line[0] != 'P' && line[0] != 'W')) continue;
            const std::size_t a = line.find('\t', 2);
            if (line[0] == 'W') {
                const std::size_t b = line.find('\t', a + 1);
                if (a == std::string::npos || b == std::string::npos) continue;
                haps.insert(line.substr(2, a - 2) + "#" + line.substr(a + 1, b - a - 1));   // sample#haplotype
            } else {
                const std::string name = line.substr(2, a == std::string::npos ? std::string::npos : a - 2);
                const std::size_t h1 = name.find('#');
                if (h1 == std::string::npos) { plain = true; continue; }
                const std::size_t h2 = name.find('#', h1 + 1);
                haps.insert(h2 == std::string::npos ? name.substr(0, h1) + "#0" : name.substr(0, h2));
            }
        }
        if (plain && haps.empty()) log("note: genome: " + ge.path + "; its paths are not PanSN names (sample#haplotype#contig): counted as one haplotype");
        else if (plain) log("note: genome: " + ge.path + "; some paths are not PanSN names: only the PanSN ones are counted as haplotypes");
        haplotypes += std::max<std::size_t>(1, haps.size());
    }
    haplotypes = std::max<std::uint64_t>(1, haplotypes);
    const std::vector<std::string> decoys = o.decoy_genomes.empty() ? std::vector<std::string>{cref} : o.decoy_genomes;
    for (const auto& d : decoys)
        if (!genomes.count(d)) throw std::runtime_error("recruit: decoy genome " + d + " is not in the background");
    // calibration genomes given that are not decoy genomes: scanned for homologous windows too, which give no decoys
    std::vector<std::string> scanned = decoys;
    for (const auto& c : o.calibration_genomes) {
        if (!genomes.count(c)) throw std::runtime_error("recruit: calibration genome " + c + " is not in the background");
        if (std::count(o.calibration_genomes.begin(), o.calibration_genomes.end(), c) > 1)
            throw std::runtime_error("recruit: calibration genome " + c + " is given twice");
        if (std::find(scanned.begin(), scanned.end(), c) == scanned.end()) scanned.push_back(c);
    }
    std::string kind;
    const auto counts = open_counts(o.background, bg::content_id(m), kind, o.log);
    log("background: " + counts->kind() + "; coordinate reference: " + cref + "; haplotypes: " + std::to_string(haplotypes));
    const anchors::ReferenceAnchors ref = anchors::load_reference((fs::path(o.background) / "anchors" / cref / "reference.anc").string());
    std::vector<std::unique_ptr<anchors::SegmentAnchors>> seg_anchors;
    std::vector<const anchors::SegmentAnchors*> sa;
    for (const auto& s : m.segments) {
        seg_anchors.push_back(std::make_unique<anchors::SegmentAnchors>((fs::path(o.background) / s.dir / ("anchors-" + cref + ".bin")).string()));
        sa.push_back(seg_anchors.back().get());
    }

    // 1-3: each locus's k-mers, their background counts, their copies in every genome
    std::vector<LocusWork> work(o.loci.size());
    std::vector<std::unique_ptr<Graph>> graphs(o.loci.size());
    // progress bars when the run logs (not --quiet)
    auto bar_label = [&](const std::string& what) { return o.log ? cli::log_prefix("recruit") + what : std::string(); };
    cli::ProgressBar locate_bar(bar_label("locating loci"), o.loci.size());
    parallel_for(o.loci.size(), threads, [&](std::size_t li) {
        LocusWork& w = work[li];
        LocusIndex& x = w.idx;
        x.graph_path = fs::absolute(o.loci[li]).string();
        x.name = locus_name(o.loci[li]);
        x.graph_checksum = file_checksum(o.loci[li]);
        bool fasta = false;
        graphs[li] = std::make_unique<Graph>(read_locus(o.loci[li], true, fasta));
        const Graph& g = *graphs[li];
        // paths are spelled by joining their segments, so a link that overlaps would invent k-mers
        for (const auto& [id, n] : g.nodes)
            for (const auto* side : {&n.start, &n.end})
                for (const auto& nb : *side)
                    if (nb.overlap > 0)
                        throw std::runtime_error("recruit: " + o.loci[li] + ": a link from " + id + " to " + nb.node_id + " overlaps by " +
                                                 std::to_string(nb.overlap) + " bases; only 0M or * links are supported");
        std::vector<std::string> names;
        for (const auto& p : g.paths) names.push_back(p.name);
        const auto hit = match_reference_name(names, o.reference_query);
        if (hit.size() != 1)
            throw std::runtime_error("recruit: " + o.loci[li] + ": -r " + o.reference_query + " names " + std::to_string(hit.size()) + " paths, not one");
        x.reference_path = names[hit[0]];
        const std::string ref_seq = spell_path_steps_sequence(g, g.paths[hit[0]].steps);
        const ParsedReferencePath rp = parse_reference_path_label(x.reference_path, ref_seq.size(), g.paths[hit[0]].type == 'W');
        if (!rp.has_interval) throw std::runtime_error("recruit: " + x.reference_path + " carries no interval (contig:start-end)");
        if (rp.length_mismatch)
            throw std::runtime_error("recruit: " + x.reference_path + ": the interval's length is not the path's, " + std::to_string(ref_seq.size()) + " bases");
        kmer::for_each_canonical_kmer(ref_seq.data(), ref_seq.size(), kK, [&](std::uint64_t, std::size_t, bool) { ++w.ref_kmers; });
        w.ref_seq = ref_seq;
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
        x.haplotypes = haplotypes;
        x.pad = o.pad;
        x.dust = o.dust;
        x.decoy_min_hits = o.decoy_min_hits;
        x.decoy_genomes = decoys;
        x.calibration_genomes_given = o.calibration_genomes;
        x.calibrated = o.calibrate;
        std::vector<std::uint32_t> present;
        collect_kmers(g, x, present, !fasta);
        x.backbone.resize(x.codes.size());
        for (std::size_t i = 0; i < x.codes.size(); ++i) x.backbone[i] = present[i] == x.paths && x.mult[i] == 1;
        w.counts = counts->count_batch(x.codes);
        w.sub.assign(x.codes.size(), 0);
        w.sub_trim.assign(x.codes.size(), 0);
        for (const auto& sp : anchors::locate(ref, sa, x.contig, x.start, x.end)) w.spans[sp.genome].push_back(sp);
        anchors::Span self;   // the coordinate reference: its exact interval
        self.genome = cref; self.contig = x.contig; self.start = x.start; self.end = x.end; self.inner_start = x.start; self.inner_end = x.end;
        w.spans[cref] = {self};
        log("locus: " + x.name + "; k-mers: " + std::to_string(x.codes.size()) + "; paths: " + std::to_string(x.paths) + "; interval: " + x.contig + ":" +
            std::to_string(x.start) + "-" + std::to_string(x.end) + "; located in genomes: " + std::to_string(w.spans.size()));
        locate_bar.tick();
    });
    locate_bar.done();

    // subtraction: read every span from its genome (genomes in parallel, every locus at once; what each finds is
    // applied in genome order, so the index does not depend on the threads)
    std::vector<FlatMap<std::uint32_t>> code_index(work.size());   // k-mer -> its index in the locus's codes
    parallel_for(work.size(), threads, [&](std::size_t li) {
        const auto& codes = work[li].idx.codes;
        code_index[li].reserve(std::max<std::size_t>(codes.size(), 1));
        for (std::size_t i = 0; i < codes.size(); ++i) {
            bool existing = false;
            code_index[li].slot(codes[i], existing) = static_cast<std::uint32_t>(i);
        }
    });
    struct Subtracted {
        std::vector<std::tuple<std::uint32_t, std::uint32_t, bool>> occ;   // (locus, k-mer, inside the walked-out copy)
        std::vector<std::pair<std::size_t, GenomeSpan>> spans;
    };
    std::vector<std::string> gnames;
    for (const auto& g : genomes) gnames.push_back(g.first);
    std::vector<Subtracted> subtracted(gnames.size());
    std::vector<char> ready(gnames.size(), 0);
    std::size_t applied = 0;   // genomes before this one are applied: each is applied once it and all before it are done
    std::mutex apply_mutex;
    // occurrences are summed as soon as a genome is done (the order does not matter); spans keep the genomes'
    // order, so they wait for the genomes before
    auto apply_ready = [&] {
        while (applied < gnames.size() && ready[applied]) {
            Subtracted& r = subtracted[applied++];
            for (auto& [li, gs] : r.spans) work[li].idx.spans.push_back(std::move(gs));
            r = Subtracted{};
        }
    };
    cli::ProgressBar cut_bar(bar_label("reading the copies from genomes"), gnames.size());
    parallel_for(gnames.size(), threads, [&](std::size_t gi) {
        const std::string& gname = gnames[gi];
        const bg::Genome& ge = genomes.at(gname);
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
        auto done = [&] {
            std::lock_guard<std::mutex> lock(apply_mutex);
            Subtracted& r = subtracted[gi];
            for (const auto& [li, i, trimmed] : r.occ) {
                ++work[li].sub[i];
                if (trimmed) ++work[li].sub_trim[i];
            }
            decltype(r.occ)().swap(r.occ);
            ready[gi] = 1;
            apply_ready();
        };
        if (iv.empty()) { done(); return; }
        check_unchanged(ge);
        std::vector<std::vector<bg::graph::Haplotypes::Piece>> pieces(iv.size());   // a FASTA span is one piece
        if (is_gfa_file(ge.path)) {
            const bg::graph::Haplotypes h(ge.path);
            for (std::size_t k = 0; k < iv.size(); ++k) pieces[k] = h.fetch(iv[k].contig, iv[k].start, iv[k].end);
        } else {
            if (!fasta::IndexedFasta::indexed(ge.path)) {
                std::lock_guard<std::mutex> lock(apply_mutex);
                log("note: genome: " + ge.path + "; no .fai (and .gzi): read in full to cut its spans");
            }
            std::vector<std::string> seqs = fasta::fetch_intervals(ge.path, iv);
            for (std::size_t k = 0; k < iv.size(); ++k) pieces[k].push_back({iv[k].start, std::move(seqs[k])});
        }
        Subtracted& out = subtracted[gi];
        for (std::size_t k = 0; k < iv.size(); ++k) {
            const std::size_t li = owner[k].first;
            LocusWork& w = work[li];
            anchors::Span& sp = w.spans.find(gname)->second[owner[k].second];   // this genome's own: no other thread
            // the locus's k-mer occurrences in the span, by position
            std::vector<std::pair<std::uint64_t, std::uint32_t>> hits;
            for (const auto& piece : pieces[k])
                kmer::for_each_canonical_kmer(piece.seq.data(), piece.seq.size(), kK, [&](std::uint64_t c, std::size_t at, bool) {
                    const std::uint32_t* i = code_index[li].find(c);
                    if (i != nullptr) hits.emplace_back(piece.start + at, *i);
                });
            std::sort(hits.begin(), hits.end());
            const auto [lo, hi] = walk_out(hits, w.counts, static_cast<std::uint64_t>(w.idx.informative_max) * w.idx.haplotypes,
                                           sp.start, sp.end, sp.inner_start, sp.inner_end);
            for (const auto& [pos, i] : hits) out.occ.emplace_back(static_cast<std::uint32_t>(li), i, pos >= lo && pos + kK <= hi);
            // the reference's own interval must hold the reference path: a label off its sequence would subtract
            // nothing there and turn the locus's k-mers shared. Where the span covers the label's interval, the
            // interval must spell the path (another path's k-mers there would pass a count); else its k-mers.
            if (gname == cref && owner[k].second == 0 && sp.contig == w.idx.contig) {
                const std::string* covering = nullptr;
                std::uint64_t at = 0;
                for (const auto& piece : pieces[k])
                    if (piece.start <= w.idx.start && piece.start + piece.seq.size() >= w.idx.end) { covering = &piece.seq; at = w.idx.start - piece.start; }
                if (covering != nullptr) {
                    auto same = [&](const std::string& want) {
                        if (want.size() != w.idx.end - w.idx.start) return false;
                        for (std::size_t q = 0; q < want.size(); ++q) {
                            const int a = std::toupper(static_cast<unsigned char>((*covering)[at + q])), b = std::toupper(static_cast<unsigned char>(want[q]));
                            const bool acgt = std::strchr("ACGT", a) != nullptr && std::strchr("ACGT", b) != nullptr;   // N or IUPAC: either way
                            if (acgt && a != b) return false;
                        }
                        return true;
                    };
                    if (!same(w.ref_seq) && !same(reverse_complement(w.ref_seq)))
                        throw std::runtime_error("recruit: " + w.idx.reference_path + ": " + cref + " " + w.idx.contig + ":" + std::to_string(w.idx.start) + "-" +
                                                 std::to_string(w.idx.end) + " does not spell the path; the label does not match the reference");
                }
            }
            if (gname == cref && owner[k].second == 0 && hits.size() < w.ref_kmers)
                throw std::runtime_error("recruit: " + w.idx.reference_path + ": " + cref + " " + sp.contig + ":" + std::to_string(sp.start) + "-" +
                                         std::to_string(sp.end) + " holds " + std::to_string(hits.size()) + " of the path's " +
                                         std::to_string(w.ref_kmers) + " k-mers; the label does not match the reference");
            GenomeSpan gs{gname, sp.contig, sp.start, sp.end, sp.forward, hits.size()};
            gs.anchored_start = sp.inner_start;
            gs.anchored_end = sp.inner_end;
            gs.locus_start = lo;
            gs.locus_end = hi;
            out.spans.emplace_back(li, std::move(gs));
            sp.inner_start = lo;   // the piece scan's uncertain ends are those past the walked-out copy
            sp.inner_end = hi;
        }
        done();
        cut_bar.tick();
    });
    cut_bar.done();
    apply_ready();

    // pieces: a genome whose copy of a locus is broken (a copy reaching a contig end, more than one copy, or
    // none) may hold more of it out of the anchors' reach, such as a tandem array cut by a contig break; a copy
    // with an end estimated from the reference (no flank anchor where expected, as past an array longer than
    // the reference's) may go on past that end on its contig. A locus whose tandem k-mers are often left with a
    // few copies outside the located ones has such pieces next to copies that look whole (an array contig of
    // its own): every genome is scanned for it. A scanned genome, or contig, is searched for stretches dense in
    // the locus's k-mers and tied to it (on a located copy's contig, or reaching a contig end): their k-mers off
    // the located copies are subtracted, and those in a copy's uncertain ends are the locus, so no longer
    // uncertain.
    {
        std::vector<char> everywhere(work.size(), 0);
        for (std::size_t li = 0; li < work.size(); ++li) {
            const LocusWork& w = work[li];
            std::size_t tandem_out = 0;
            for (std::size_t i = 0; i < w.idx.codes.size(); ++i) {
                if (w.idx.mult[i] < kTandemCopies) continue;
                const std::uint64_t out = w.counts[i] - std::min(w.sub[i], w.counts[i]);
                if (out >= 1 && out <= kTandemOutside * w.idx.mult[i]) ++tandem_out;
            }
            if (static_cast<double>(tandem_out) >= kTandemShare * static_cast<double>(w.idx.codes.size())) {
                everywhere[li] = 1;
                log("locus: " + w.idx.name + "; tandem k-mers with a few copies outside: " + std::to_string(tandem_out) + "; every genome scanned for pieces");
            }
        }
        struct Need {
            std::vector<std::size_t> whole;                               // loci looked for in every contig
            std::map<std::string, std::vector<std::size_t>> contigs;     // contig -> loci looked for in it
            std::map<std::string, std::uint64_t> lengths;
        };
        std::vector<std::pair<std::string, Need>> scan_list;
        for (const auto& g : genomes) {
            const std::string& gname = g.first;
            if (gname == cref) continue;
            Need nd;
            for (std::size_t li = 0; li < work.size(); ++li) {
                const auto it = work[li].spans.find(gname);
                if (everywhere[li] || it == work[li].spans.end() || it->second.size() > 1 ||
                    std::any_of(it->second.begin(), it->second.end(), [](const anchors::Span& sp) { return sp.open; })) {
                    nd.whole.push_back(li);
                    continue;
                }
                for (const auto& sp : it->second)
                    if (sp.estimated) {
                        nd.contigs[sp.contig].push_back(li);
                        nd.lengths[sp.contig] = sp.contig_length;
                    }
            }
            if (!nd.whole.empty() || !nd.contigs.empty()) scan_list.emplace_back(gname, std::move(nd));
        }
        std::vector<char> wanted(work.size(), 0);
        for (const auto& e : scan_list) {
            for (const std::size_t li : e.second.whole) wanted[li] = 1;
            for (const auto& c : e.second.contigs) for (const std::size_t li : c.second) wanted[li] = 1;
        }
        FlatMap<std::uint32_t> lk;                       // locus k-mer -> an entry of lk_owners, the loci holding it
        std::vector<std::vector<std::uint32_t>> lk_owners;
        {
            std::size_t n = 0;
            for (std::size_t li = 0; li < work.size(); ++li) if (wanted[li]) n += work[li].idx.codes.size();
            lk.reserve(std::max<std::size_t>(n, 1));
            for (std::size_t li = 0; li < work.size(); ++li) {
                if (!wanted[li]) continue;
                for (const std::uint64_t c : work[li].idx.codes) {
                    bool existing = false;
                    std::uint32_t& id = lk.slot(c, existing);
                    if (!existing) { id = static_cast<std::uint32_t>(lk_owners.size()); lk_owners.emplace_back(); }
                    lk_owners[id].push_back(static_cast<std::uint32_t>(li));
                }
            }
        }
        const std::size_t need = static_cast<std::size_t>(std::ceil(kPieceDensity * static_cast<double>(kPieceWindow)));
        const std::size_t sampled_need = kPieceWindow / kPieceSample - (kPieceWindow - need);   // see kPieceSample
        // what a genome's scan finds, applied in genome order once all are scanned (genomes run in parallel)
        struct Found {
            std::vector<std::pair<std::uint32_t, std::uint32_t>> off, ends;   // (locus, k-mer): off the copies / in an uncertain end
            std::vector<std::pair<std::size_t, GenomeSpan>> pieces;
        };
        // `seq` is [offset, offset + seq.size()) of a contig `len` bases long. A scan of part of a contig sets grow[0]
        // (grow[1]) when a piece reaches within kPieceMargin of the part's start (end) inside the contig: the part
        // must be widened there and scanned again
        auto scan = [&](const std::string& gname, const std::string& name, const std::string& seq, std::uint64_t offset,
                        std::uint64_t len, const std::vector<char>& here, std::size_t inner, Found& f, char* grow) {
            if (seq.size() < kK) return;
            const bool whole_only = len < kPieceWindow + kK - 1;   // no window fits
            const bool at_start = offset == 0, at_end = offset + seq.size() == len;
            // each locus's k-mer starts on the contig, in order: ranges scanned in parallel, joined in order. On a
            // contig longer than a window, only every kPieceSample-th start is looked up first; a dense window
            // (need of kPieceWindow starts) keeps at least sampled_need of its sampled starts, so only around
            // windows that do is every start looked up: the same stretches, an eighth of the lookups
            const std::size_t starts = seq.size() - kK + 1;
            const std::size_t parts = std::min<std::size_t>(inner, std::max<std::size_t>(1, starts / 65536));
            const std::uint64_t step = whole_only ? 1 : kPieceSample;
            std::vector<std::vector<std::vector<std::uint64_t>>> part_hits(parts, std::vector<std::vector<std::uint64_t>>(work.size()));
            auto scan_part = [&](std::size_t p) {
                const std::size_t a = starts * p / parts, b = starts * (p + 1) / parts;
                if (b <= a) return;
                kmer::for_each_canonical_kmer(seq.data() + a, b - a + kK - 1, kK, [&](std::uint64_t c, std::size_t off, bool) {
                    if ((a + off) % step != 0) return;
                    const std::uint32_t* id = lk.find(c);
                    if (id == nullptr) return;
                    for (const std::uint32_t li : lk_owners[*id]) if (here[li]) part_hits[p][li].push_back(a + off);
                });
            };
            parallel_for(parts, parts, scan_part);
            for (std::size_t li = 0; li < work.size(); ++li) {
                if (!here[li]) continue;
                std::vector<std::uint64_t> hits;
                for (std::size_t p = 0; p < parts; ++p) hits.insert(hits.end(), part_hits[p][li].begin(), part_hits[p][li].end());
                if (!whole_only) {
                    // candidate windows from the sampled starts, widened by two windows each way (a stretch's
                    // reach to a contig end is decided within a window of it), then every start in them
                    std::vector<std::pair<std::uint64_t, std::uint64_t>> regions;
                    for (std::size_t i = 0, j = 0; i < hits.size(); ++i) {
                        if (j < i) j = i;
                        while (j < hits.size() && hits[j] < hits[i] + kPieceWindow) ++j;
                        if (j - i < sampled_need) continue;
                        const std::uint64_t ra = hits[i] > 2 * kPieceWindow ? hits[i] - 2 * kPieceWindow : 0;
                        const std::uint64_t rb = std::min<std::uint64_t>(starts, hits[i] + 3 * kPieceWindow);
                        if (!regions.empty() && ra <= regions.back().second) regions.back().second = std::max(regions.back().second, rb);
                        else regions.emplace_back(ra, rb);
                    }
                    hits.clear();
                    for (const auto& r : regions) {
                        const std::uint64_t ra = r.first, rb = r.second;
                        kmer::for_each_canonical_kmer(seq.data() + ra, rb - ra + kK - 1, kK, [&](std::uint64_t c, std::size_t off, bool) {
                            const std::uint32_t* id = lk.find(c);
                            if (id == nullptr) return;
                            for (const std::uint32_t x : lk_owners[*id]) if (x == li) { hits.push_back(ra + off); break; }
                        });
                    }
                }
                std::vector<std::pair<std::uint64_t, std::uint64_t>> dense;
                if (whole_only) {
                    // a contig shorter than a window is a piece whole when at least half of its k-mer starts are the
                    // locus's, the rule for a stretch's last bases before a contig end
                    if (hits.empty() || 2 * hits.size() < starts) continue;
                    dense.emplace_back(0, seq.size());
                } else {
                    if (hits.size() < need) continue;
                    // dense stretches: the union of the windows of kPieceWindow k-mer starts holding enough of them
                    for (std::size_t i = 0, j = 0; i < hits.size(); ++i) {
                        if (j < i) j = i;
                        while (j < hits.size() && hits[j] < hits[i] + kPieceWindow) ++j;
                        if (j - i < need) continue;
                        const std::uint64_t ws = hits[i], we = hits[j - 1] + kK;
                        if (!dense.empty() && ws <= dense.back().second) dense.back().second = std::max(dense.back().second, we);
                        else dense.emplace_back(ws, we);
                    }
                    // a stretch less than a window from a contig end reaches that end when at least half of the k-mer
                    // starts in between are the locus's: no window fits there, and assemblies hold more errors near
                    // their contig ends
                    if (!dense.empty()) {
                        auto& first = dense.front();
                        if (at_start && first.first > 0 && first.first < kPieceWindow) {
                            const std::size_t n = static_cast<std::size_t>(std::lower_bound(hits.begin(), hits.end(), first.first) - hits.begin());
                            if (2 * n >= first.first) first.first = 0;
                        }
                        auto& last = dense.back();
                        if (at_end && last.second < seq.size() && seq.size() - last.second < kPieceWindow) {
                            const std::uint64_t from = last.second - kK + 1;   // the first start past the stretch's last k-mer
                            const std::size_t n = static_cast<std::size_t>(hits.end() - std::upper_bound(hits.begin(), hits.end(), last.second - kK));
                            if (2 * n >= starts - from) last.second = seq.size();
                        }
                    }
                }
                const LocusWork& w = work[li];
                std::vector<const anchors::Span*> own;   // the located copies on this contig
                const auto it = w.spans.find(gname);
                if (it != w.spans.end()) for (const auto& sp : it->second) if (sp.contig == name) own.push_back(&sp);
                // only a stretch tied to the locus is a piece: reaching a contig end (an array broken across
                // contigs, a stub, a contig of its own), or next to a located copy, directly or through other
                // such stretches (an array running past the copy's end). One apart from every copy with other
                // sequence on both sides may be a paralog the coordinate reference lacks, even on a copy's
                // contig, which can hold a whole chromosome: it stays outside (shared, never falsely private).
                for (auto& d : dense) { d.first += offset; d.second += offset; }   // contig coordinates from here
                std::vector<char> tied(dense.size(), 0);
                for (std::size_t si = 0; si < dense.size(); ++si) {
                    const std::uint64_t a = dense[si].first, b = dense[si].second;
                    tied[si] = a == 0 || b == len || std::any_of(own.begin(), own.end(), [&](const anchors::Span* sp) {
                        return a <= sp->end + kPieceGap && sp->start <= b + kPieceGap;
                    });
                }
                // through neighbours: the stretches are sorted and disjoint
                for (std::size_t si = 1; si < dense.size(); ++si)
                    if (tied[si - 1] && dense[si].first <= dense[si - 1].second + kPieceGap) tied[si] = 1;
                for (std::size_t si = dense.size(); si-- > 1;)
                    if (tied[si] && dense[si].first <= dense[si - 1].second + kPieceGap) tied[si - 1] = 1;
                for (std::size_t si = 0; si < dense.size(); ++si) {
                    if (!tied[si] || grow == nullptr) continue;
                    if (!at_start && dense[si].first < offset + kPieceMargin) grow[0] = 1;
                    if (!at_end && dense[si].second + kPieceMargin > offset + seq.size()) grow[1] = 1;
                }
                for (std::size_t si = 0; si < dense.size(); ++si) {
                    if (!tied[si]) continue;
                    const std::uint64_t a = dense[si].first, b = dense[si].second;
                    // the stretch off the copies, as pieces overlapping a neighbouring copy by k - 1 bases: a k-mer
                    // lies wholly in exactly one copy or piece
                    std::vector<std::pair<std::uint64_t, std::uint64_t>> off{{a, b}};
                    for (const anchors::Span* sp : own) {
                        std::vector<std::pair<std::uint64_t, std::uint64_t>> next;
                        for (const auto& x : off) {
                            if (sp->end <= x.first || sp->start >= x.second) { next.push_back(x); continue; }
                            if (x.first < sp->start) next.emplace_back(x.first, std::min(x.second, sp->start + kK - 1));
                            if (sp->end < x.second) next.emplace_back(std::max(x.first, sp->end - (kK - 1)), x.second);
                        }
                        off = std::move(next);
                    }
                    std::vector<std::uint64_t> found(off.size(), 0);
                    kmer::for_each_canonical_kmer(seq.data() + (a - offset), b - a, kK, [&](std::uint64_t c, std::size_t at, bool) {
                        const std::int64_t i = w.idx.find(c);
                        if (i < 0) return;
                        const std::uint64_t pos = a + at;
                        const auto key = std::make_pair(static_cast<std::uint32_t>(li), static_cast<std::uint32_t>(i));
                        for (const anchors::Span* sp : own)
                            if (pos >= sp->start && pos + kK <= sp->end) {   // counted with its copy
                                if (pos < sp->inner_start || pos + kK > sp->inner_end) f.ends.push_back(key);
                                return;
                            }
                        f.off.push_back(key);
                        for (std::size_t k = 0; k < off.size(); ++k)
                            if (pos >= off[k].first && pos + kK <= off[k].second) { ++found[k]; break; }
                    });
                    for (std::size_t k = 0; k < off.size(); ++k)
                        if (found[k] > 0) f.pieces.emplace_back(li, GenomeSpan{gname, name, off[k].first, off[k].second, true, found[k], true});
                }
            }
        };
        // genomes in parallel, each read by one thread; a lone genome splits its contigs among the threads
        std::vector<Found> found(scan_list.size());
        const std::size_t outer = std::max<std::size_t>(1, std::min(threads, scan_list.size()));
        const std::size_t inner = std::max<std::size_t>(1, threads / outer);
        std::mutex log_mutex;
        for (const auto& e : scan_list) check_unchanged(genomes.at(e.first));
        cli::ProgressBar scan_bar(bar_label("pieces: scanning genomes"), scan_list.size());
        parallel_for(scan_list.size(), outer, [&](std::size_t gi) {
            const std::string& gname = scan_list[gi].first;
            const Need& nd = scan_list[gi].second;
            const bg::Genome& ge = genomes.at(gname);
            if (is_gfa_file(ge.path)) {
                std::lock_guard<std::mutex> lock(log_mutex);
                log("note: genome: " + gname + "; a graph: not scanned for pieces of a locus");
                return;
            }
            auto here_for = [&](const std::string& contig) {
                std::vector<char> h(work.size(), 0);
                for (const std::size_t li : nd.whole) h[li] = 1;
                const auto it = nd.contigs.find(contig);
                if (it != nd.contigs.end()) for (const std::size_t li : it->second) h[li] = 1;
                return h;
            };
            if (fasta::IndexedFasta::indexed(ge.path)) {
                fasta::IndexedFasta fa(ge.path);
                const std::vector<std::string> names = [&] {
                    std::vector<std::string> v;
                    if (!nd.whole.empty()) return fa.contigs();
                    for (const auto& c : nd.contigs) v.push_back(c.first);
                    return v;
                }();
                for (const std::string& name : names) {
                    const std::uint64_t len = fa.length(name);
                    const std::vector<char> here = here_for(name);
                    std::vector<std::pair<std::uint64_t, std::uint64_t>> parts{{0, std::min(len, kPiecePad)},
                                                                               {len > kPiecePad ? len - kPiecePad : 0, len}};
                    for (std::size_t li = 0; li < work.size(); ++li) {
                        const auto it = here[li] ? work[li].spans.find(gname) : work[li].spans.end();
                        if (it == work[li].spans.end()) continue;
                        for (const auto& sp : it->second)
                            if (sp.contig == name) parts.emplace_back(sp.start > kPiecePad ? sp.start - kPiecePad : 0, std::min(len, sp.end + kPiecePad));
                    }
                    while (true) {
                        std::sort(parts.begin(), parts.end());
                        std::vector<std::pair<std::uint64_t, std::uint64_t>> merged;
                        for (const auto& x : parts) {
                            if (!merged.empty() && x.first <= merged.back().second + 2 * kPieceMargin) merged.back().second = std::max(merged.back().second, x.second);
                            else merged.push_back(x);
                        }
                        parts = std::move(merged);
                        Found here_found;
                        bool grew = false;
                        for (auto& x : parts) {
                            char grow[2] = {0, 0};
                            scan(gname, name, fa.fetch(name, x.first, x.second), x.first, len, here, inner, here_found, grow);
                            const std::uint64_t w = std::max(x.second - x.first, kPiecePad);
                            if (grow[0]) x.first = x.first > w ? x.first - w : 0;
                            if (grow[1]) x.second = std::min(len, x.second + w);
                            grew = grew || grow[0] || grow[1];
                        }
                        if (grew) continue;
                        // in the order a scan of the whole contig gives: by locus, then along the contig
                        std::stable_sort(here_found.pieces.begin(), here_found.pieces.end(), [](const auto& x, const auto& y) {
                            return x.first != y.first ? x.first < y.first : x.second.start < y.second.start;
                        });
                        Found& f = found[gi];
                        f.off.insert(f.off.end(), here_found.off.begin(), here_found.off.end());
                        f.ends.insert(f.ends.end(), here_found.ends.begin(), here_found.ends.end());
                        for (auto& pc : here_found.pieces) f.pieces.push_back(std::move(pc));
                        break;
                    }
                }
            } else if (!nd.whole.empty()) {
                GzLineReader in(ge.path);
                if (!in.ok()) throw std::runtime_error("recruit: cannot read " + ge.path);
                std::string line, name, seq;
                while (in.getline(line)) {
                    if (line.empty()) continue;
                    if (line[0] == '>') {
                        if (!name.empty()) scan(gname, name, seq, 0, seq.size(), here_for(name), inner, found[gi], nullptr);
                        const auto sp = line.find_first_of(" \t");
                        name = line.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
                        seq.clear();
                        continue;
                    }
                    seq += line;
                }
                if (!name.empty()) scan(gname, name, seq, 0, seq.size(), here_for(name), inner, found[gi], nullptr);
            } else {
                std::vector<fasta::Interval> iv;
                for (const auto& c : nd.contigs) iv.push_back({c.first, 0, nd.lengths.at(c.first)});
                const std::vector<std::string> seqs = fasta::fetch_intervals(ge.path, iv);
                for (std::size_t k = 0; k < iv.size(); ++k)
                    scan(gname, iv[k].contig, seqs[k], 0, seqs[k].size(), here_for(iv[k].contig), inner, found[gi], nullptr);
            }
            std::lock_guard<std::mutex> lock(log_mutex);
            // the occurrences counted as each genome is done, not held for every genome until the last
            for (const auto& [li, i] : found[gi].off) { ++work[li].sub[i]; ++work[li].sub_trim[i]; }
            for (const auto& [li, i] : found[gi].ends) ++work[li].sub_trim[i];
            decltype(found[gi].off)().swap(found[gi].off);
            decltype(found[gi].ends)().swap(found[gi].ends);
            log(nd.whole.empty() ? "pieces: contigs scanned: " + gname + " (" + std::to_string(nd.contigs.size()) + ")"
                                 : "pieces: genome scanned: " + gname);
            scan_bar.tick();
        });
        scan_bar.done();
        std::vector<std::size_t> scanned_for(work.size(), 0);
        std::vector<std::set<std::string>> piece_genomes(work.size());
        for (std::size_t gi = 0; gi < scan_list.size(); ++gi) {
            const Need& nd = scan_list[gi].second;
            std::vector<char> counted(work.size(), 0);
            for (const std::size_t li : nd.whole) counted[li] = 1;
            for (const auto& c : nd.contigs) for (const std::size_t li : c.second) counted[li] = 1;
            for (std::size_t li = 0; li < work.size(); ++li) scanned_for[li] += counted[li];
            for (auto& [li, sp] : found[gi].pieces) {
                piece_genomes[li].insert(sp.genome);
                // a piece is the locus wherever a copy is: no decoy window or off-target read is taken from it
                anchors::Span as;
                as.genome = sp.genome;
                as.contig = sp.contig;
                as.start = as.inner_start = sp.start;
                as.end = as.inner_end = sp.end;
                as.forward = sp.forward;
                work[li].spans[sp.genome].push_back(std::move(as));
                work[li].idx.spans.push_back(std::move(sp));
            }
        }
        for (std::size_t li = 0; li < work.size(); ++li) {
            if (scanned_for[li] == 0) continue;
            std::uint64_t n = 0, bp = 0;
            for (const auto& sp : work[li].idx.spans) if (sp.piece) { ++n; bp += sp.end - sp.start; }
            log("locus: " + work[li].idx.name + "; genomes scanned for pieces: " + std::to_string(scanned_for[li]) + "; pieces: " + std::to_string(n) +
                " (" + std::to_string(bp) + " bp) in " + std::to_string(piece_genomes[li].size()) + " genomes");
        }
    }

    // classes
    for (LocusWork& w : work) {
        LocusIndex& x = w.idx;
        const std::size_t n = x.codes.size();
        x.cls.resize(n);
        x.outside.resize(n);
        const std::uint64_t shared_max = static_cast<std::uint64_t>(x.informative_max) * x.haplotypes;
        auto class_of = [&](std::uint64_t out) { return out == 0 ? kPrivate : out <= shared_max ? kShared : kRepeat; };
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
        if (const char* f = std::getenv("PANVAR_TEST_FAIL_LOCUS"); f != nullptr && x.name == f) x.negative = std::max<std::uint64_t>(x.negative, 1);
        if (x.negative > 0) {   // this locus fails; the others go on, their indexes are kept
            w.error = "recruit: " + x.name + ": " + std::to_string(x.negative) +
                      " k-mers occur more often in the located copies than the background counted: a genome file differs from the one counted";
            log(w.error);
        }
    }

    // decoys: the homologous windows of the decoy genomes around their S k-mers, outside the locus
    // S k-mer -> an entry of `owners`, the loci holding it
    std::vector<std::vector<std::uint32_t>> owners;
    FlatMap<std::uint32_t> shared;
    {
        std::size_t n = 0;
        for (const auto& w : work) n += static_cast<std::size_t>(std::count(w.idx.cls.begin(), w.idx.cls.end(), kShared));
        shared.reserve(n);
        for (std::size_t li = 0; li < work.size(); ++li)
            for (std::size_t i = 0; i < work[li].idx.codes.size(); ++i) {
                if (work[li].idx.cls[i] != kShared) continue;
                bool existing = false;
                std::uint32_t& id = shared.slot(work[li].idx.codes[i], existing);
                if (!existing) {
                    id = static_cast<std::uint32_t>(owners.size());
                    owners.emplace_back();
                }
                owners[id].push_back(static_cast<std::uint32_t>(li));
            }
    }
    const std::size_t scan_threads = threads;
    std::vector<std::unordered_set<std::uint64_t>> candidates(work.size());
    cli::ProgressBar decoy_bar(bar_label("decoys: scanning genomes"), scanned.size());
    for (const auto& dname : scanned) {
        const bg::Genome& ge = genomes.at(dname);
        const bool is_decoy = std::find(decoys.begin(), decoys.end(), dname) != decoys.end();
        if (is_gfa_file(ge.path))
            throw std::runtime_error("recruit: " + std::string(is_decoy ? "decoy" : "calibration") + " genome " + dname + " is a graph; it must be FASTA");
        check_unchanged(ge);
        GzLineReader in(ge.path);
        if (!in.ok()) throw std::runtime_error("recruit: cannot read " + ge.path);
        std::string line, name, seq;
        auto scan = [&] {
            if (seq.empty()) return;
            // hits of S k-mers outside each locus's own copy, in position order: the contig's k-mer starts
            // are split into ranges scanned in parallel and joined in order
            std::vector<std::vector<std::uint64_t>> hits(work.size());
            const std::size_t starts = seq.size() >= kK ? seq.size() - kK + 1 : 0;
            const std::size_t parts = std::min<std::size_t>(scan_threads, std::max<std::size_t>(1, starts / 65536));
            std::vector<std::vector<std::vector<std::uint64_t>>> part_hits(parts, std::vector<std::vector<std::uint64_t>>(work.size()));
            auto scan_part = [&](std::size_t p) {
                const std::size_t a = starts * p / parts, b = starts * (p + 1) / parts;
                if (b <= a) return;
                kmer::for_each_canonical_kmer(seq.data() + a, b - a + kK - 1, kK, [&](std::uint64_t c, std::size_t off, bool) {
                    const std::uint32_t* id = shared.find(c);
                    if (id == nullptr) return;
                    const std::size_t at = a + off;
                    for (const std::uint32_t li : owners[*id]) {
                        bool own = false;   // inside this locus's own copy in this genome
                        const auto sp = work[li].spans.find(dname);
                        if (sp != work[li].spans.end())
                            for (const auto& s : sp->second) own = own || (s.contig == name && at + kK > s.start && at < s.end);
                        if (!own) part_hits[p][li].push_back(at);
                    }
                });
            };
            if (starts > 0) {
                parallel_for(parts, parts, scan_part);
                for (std::size_t li = 0; li < work.size(); ++li)
                    for (std::size_t p = 0; p < parts; ++p)
                        hits[li].insert(hits[li].end(), part_hits[p][li].begin(), part_hits[p][li].end());
            }
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
                    if (!is_decoy) {   // a calibration genome's homologous window: read in calibration, no decoys
                        work[li].idx.calibration_window_list.push_back({dname, name, ws, we, count});
                        return;
                    }
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
        log(std::string(is_decoy ? "decoys" : "calibration windows") + ": genome scanned: " + dname);
        decoy_bar.tick();
    }
    decoy_bar.done();
    parallel_for(work.size(), threads, [&](std::size_t li) {
        LocusIndex& x = work[li].idx;
        std::vector<std::uint64_t> cand(candidates[li].begin(), candidates[li].end());
        std::sort(cand.begin(), cand.end());
        const auto cc = counts->count_batch(cand);
        for (std::size_t i = 0; i < cand.size(); ++i)
            if (cc[i] <= static_cast<std::uint64_t>(x.informative_max) * x.haplotypes) {
                x.decoy.push_back(cand[i]);
                work[li].decoy_counts.push_back(cc[i]);
            }
        // screen: closed syncmers among P and S
        for (std::size_t i = 0; i < x.codes.size(); ++i)
            if (x.cls[i] != kRepeat && kmer::is_closed_syncmer(x.codes[i], kK, kS)) x.screen.push_back(x.codes[i]);
        compute_ceilings(*graphs[li], x);
        std::size_t np = 0, ns = 0, nr = 0;
        for (const auto c : x.cls) (c == kPrivate ? np : c == kShared ? ns : nr)++;
        log("locus: " + x.name + "; private: " + std::to_string(np) + "; shared: " + std::to_string(ns) + "; repeat: " + std::to_string(nr) +
            "; decoys: " + std::to_string(x.decoy.size()) + " from " + std::to_string(x.decoy_windows) + " windows (" + std::to_string(x.decoy_bases) +
            " bp); screen: " + std::to_string(x.screen.size()));
    });
    if (o.calibrate) {
        // off-target sequence: the coordinate reference beside each locus, and the decoy windows with context
        // (cut where the locus's own copy begins, as the windows are)
        std::vector<std::vector<OffTarget>> off(work.size());
        std::map<std::string, std::vector<fasta::Interval>> want;
        std::map<std::string, std::vector<std::pair<std::size_t, std::size_t>>> slot;   // genome -> (locus, entry)
        auto add = [&](std::size_t li, const std::string& genome, const std::string& contig, std::uint64_t s, std::uint64_t e, bool decoy,
                       std::uint64_t core_s, std::uint64_t core_e) {
            if (e <= s) return;
            OffTarget t;
            t.decoy = decoy;
            t.core_start = core_s - s;
            t.core_end = core_e - s;
            off[li].push_back(std::move(t));
            want[genome].push_back({contig, s, e});
            slot[genome].emplace_back(li, off[li].size() - 1);
        };
        for (std::size_t li = 0; li < work.size(); ++li) {
            const LocusIndex& x = work[li].idx;
            if (!is_gfa_file(genomes.at(cref).path)) {
                const std::uint64_t a = x.start > kCalibrationFlank ? x.start - kCalibrationFlank : 0;
                add(li, cref, x.contig, a, x.start, false, a, x.start);
                add(li, cref, x.contig, x.end, x.end + kCalibrationFlank, false, x.end, x.end + kCalibrationFlank);
            }
            for (const DecoyWindow& w : x.decoy_window_list) {
                std::uint64_t s = w.start > kCalibrationContext ? w.start - kCalibrationContext : 0, e = w.end + kCalibrationContext;
                const auto own = work[li].spans.find(w.genome);
                if (own != work[li].spans.end())
                    for (const auto& sp : own->second) {
                        if (sp.contig != w.contig || sp.end <= s || sp.start >= e) continue;
                        if (sp.end <= w.start) s = std::max<std::uint64_t>(s, sp.end);
                        else if (sp.start >= w.end) e = std::min<std::uint64_t>(e, sp.start);
                    }
                add(li, w.genome, w.contig, s, e, true, w.start, w.end);
            }
        }
        for (const auto& [gname, iv] : want) {
            check_unchanged(genomes.at(gname));
            std::vector<std::string> seqs = fasta::fetch_intervals(genomes.at(gname).path, iv);
            for (std::size_t k = 0; k < iv.size(); ++k) {
                OffTarget& t = off[slot[gname][k].first][slot[gname][k].second];
                t.seq = std::move(seqs[k]);
                if (!t.decoy) t.core_end = t.seq.size();   // a flank clamped at its contig's end
                t.core_end = std::min<std::uint64_t>(t.core_end, t.seq.size());
                t.core_start = std::min(t.core_start, t.core_end);
            }
        }
        // calibration genomes: each held out of the background in turn, its homologous windows and flanks read as a
        // new sample's off-target sequence (see read_calibration_genomes)
        std::vector<std::vector<GenomeOffTarget>> held_genomes(work.size());
        read_calibration_genomes(work, genomes, cref, decoys, o.calibration_genomes, code_index, threads, log, held_genomes);
        cli::ProgressBar cal_bar(bar_label("calibrating loci"), work.size());
        for (std::size_t li = 0; li < work.size(); ++li) {
            struct Tick { cli::ProgressBar& bar; ~Tick() { bar.tick(); } } tick{cal_bar};   // as each locus ends
            LocusIndex& x = work[li].idx;
            if (!work[li].error.empty()) continue;
            try {
                calibrate(x, *graphs[li], off[li], held_genomes[li], threads);
            } catch (const std::exception& e) {
                work[li].error = "recruit: " + x.name + ": calibration: " + e.what();
                log(work[li].error);
                continue;
            }
            if (std::none_of(graphs[li]->paths.begin(), graphs[li]->paths.end(), [](const PathRecord& p) { return p.name.find('#') != std::string::npos; }))
                log("note: locus: " + x.name + "; no path name is PanSN (sample#haplotype#contig): calibration holds out one path at a time, "
                    "so a sample's other haplotype stays in and held-out recall is optimistic");
            for (const Calibration& c : x.calibration) {
                const std::uint64_t scr = c.flank_screened + c.decoy_screened;
                char b[512];
                std::snprintf(b, sizeof b, "locus: %s; calibration: %s reads; held-out recall at 2 sites: %.4f (%llu fragments, %llu samples); off-target certain at 1 site: %llu, at 2: %llu; screened off-target: %llu",
                              x.name.c_str(), c.read_class.c_str(), c.locus ? static_cast<double>(c.locus_certain[2]) / static_cast<double>(c.locus) : 0.0,
                              static_cast<unsigned long long>(c.locus), static_cast<unsigned long long>(c.held_out),
                              static_cast<unsigned long long>(c.flank_certain[1] + c.decoy_certain[1]),
                              static_cast<unsigned long long>(c.flank_certain[2] + c.decoy_certain[2]), static_cast<unsigned long long>(scr));
                log(b);
            }
        }
    }
    std::vector<LocusIndex> out;
    for (auto& w : work) {
        if (!w.error.empty() && o.errors == nullptr) throw std::runtime_error(w.error);
        out.push_back(std::move(w.idx));
    }
    if (o.errors != nullptr) {
        o.errors->clear();
        for (const auto& w : work) o.errors->push_back(w.error);
    }
    return out;
}

void save(const LocusIndex& x, const std::string& path) {
    const std::string tmp = path + ".tmp." + std::to_string(::getpid());   // two jobs writing one index do not share it
    {
        Out o{std::ofstream(tmp, std::ios::binary)};
        if (!o.f) throw std::runtime_error("recruit: cannot write " + tmp);
        o.raw64(0);
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
        for (const auto& s : x.spans) { o.str(s.genome); o.str(s.contig); o.u64(s.start); o.u64(s.end); o.u64(s.forward); o.u64(s.subtracted); o.u64(s.piece);
                                     o.u64(s.anchored_start); o.u64(s.anchored_end); o.u64(s.locus_start); o.u64(s.locus_end); }
        o.vec(x.ceilings);
        o.u64(x.dust_demoted); o.u64(x.span_demoted); o.u64(x.negative);
        o.u64(x.decoy_min_hits); o.u64(x.decoy_windows); o.u64(x.decoy_bases); o.u64(x.haplotypes);
        o.u64(x.decoy_window_list.size());
        for (const auto& w : x.decoy_window_list) { o.str(w.genome); o.str(w.contig); o.u64(w.start); o.u64(w.end); o.u64(w.hits); }
        o.u64(x.calibrated);
        o.u64(x.calibration.size());
        for (const auto& c : x.calibration) {
            o.str(c.read_class);
            std::uint64_t ebits;
            std::memcpy(&ebits, &c.error_rate, 8);
            for (const std::uint64_t v : {ebits, c.read_length, c.insert, c.held_out, c.locus, c.flank, c.decoy, c.flank_screened, c.decoy_screened}) o.u64(v);
            for (const auto* a : {&c.locus_certain, &c.flank_certain, &c.decoy_certain})
                for (const std::uint64_t v : *a) o.u64(v);
            for (const double r : {c.r_in, c.r_out, c.q_in, c.q_out}) o.f64(r);
            o.u64(c.genome_fragments); o.u64(c.genome_screened);
            for (const ScoreHistogram* h : {&c.locus_score, &c.off_score, &c.genome_score}) { o.vec(h->bin); o.vec(h->count); }
            o.u64(c.genomes.size());
            for (const auto& g : c.genomes) { o.str(g.genome); o.u64(g.fragments); o.u64(g.screened); o.u64(g.with_private); }
        }
        o.u64(x.calibration_genomes_given.size());
        for (const auto& g : x.calibration_genomes_given) o.str(g);
        o.u64(x.calibration_window_list.size());
        for (const auto& w : x.calibration_window_list) { o.str(w.genome); o.str(w.contig); o.u64(w.start); o.u64(w.end); o.u64(w.hits); }
        o.u64(x.calibration_genomes.size());
        for (const auto& g : x.calibration_genomes) {
            o.str(g.genome);
            for (const std::uint64_t v : {g.windows, g.bases, g.homolog_kmers, g.relabelled, g.unlabelled}) o.u64(v);
        }
        o.raw64(o.h);
        o.f.seekp(0);
        o.raw64(kMagic);
        o.f.flush();
        if (!o.f) throw std::runtime_error("recruit: cannot write " + tmp);
    }
    fs::rename(tmp, path);
}

LocusIndex load(const std::string& path) {
    In i{std::ifstream(path, std::ios::binary), path};
    if (!i.f) throw std::runtime_error("recruit: cannot read " + path);
    if (i.raw64() != kMagic) throw std::runtime_error("recruit: " + path + " is not a locus index (or was interrupted)");
    if (const std::uint64_t v = i.u64(); v != kVersion)
        throw std::runtime_error("recruit: " + path + " is a locus index of format " + std::to_string(v) + "; this panvar reads format " +
                                 std::to_string(kVersion) + (v < kVersion ? ", which calibrates the score: rebuild it" : ": use the panvar that built it, or rebuild it"));
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
    x.decoy_genomes.resize(i.length(8));
    for (auto& d : x.decoy_genomes) d = i.str();
    x.paths = i.u64();
    x.codes = i.vec<std::uint64_t>(); x.cls = i.vec<std::uint8_t>(); x.outside = i.vec<std::uint32_t>(); x.mult = i.vec<std::uint16_t>();
    x.backbone = i.vec<std::uint8_t>(); x.nodeset = i.vec<std::uint32_t>();
    x.nodeset_first = i.vec<std::uint32_t>(); x.nodeset_nodes = i.vec<std::uint32_t>();
    x.nodes.resize(i.length(8));
    for (auto& n : x.nodes) n = i.str();
    x.decoy = i.vec<std::uint64_t>(); x.screen = i.vec<std::uint64_t>();
    x.spans.resize(i.length(8));
    for (auto& s : x.spans) { s.genome = i.str(); s.contig = i.str(); s.start = i.u64(); s.end = i.u64(); s.forward = i.u64() != 0; s.subtracted = i.u64(); s.piece = i.u64() != 0;
                                      s.anchored_start = i.u64(); s.anchored_end = i.u64(); s.locus_start = i.u64(); s.locus_end = i.u64(); }
    x.ceilings = i.vec<Ceiling>();
    x.dust_demoted = i.u64(); x.span_demoted = i.u64(); x.negative = i.u64();
    x.decoy_min_hits = static_cast<std::uint32_t>(i.u64()); x.decoy_windows = i.u64(); x.decoy_bases = i.u64(); x.haplotypes = i.u64();
    x.decoy_window_list.resize(i.length(8));
    for (auto& w : x.decoy_window_list) { w.genome = i.str(); w.contig = i.str(); w.start = i.u64(); w.end = i.u64(); w.hits = i.u64(); }
    x.calibrated = i.u64() != 0;
    x.calibration.resize(i.length(8));
    for (auto& c : x.calibration) {
        c.read_class = i.str();
        const std::uint64_t ebits = i.u64();
        std::memcpy(&c.error_rate, &ebits, 8);
        for (std::uint64_t* v : {&c.read_length, &c.insert, &c.held_out, &c.locus, &c.flank, &c.decoy, &c.flank_screened, &c.decoy_screened}) *v = i.u64();
        for (auto* a : {&c.locus_certain, &c.flank_certain, &c.decoy_certain})
            for (std::uint64_t& v : *a) v = i.u64();
        for (double* r : {&c.r_in, &c.r_out, &c.q_in, &c.q_out}) *r = i.f64();
        c.genome_fragments = i.u64(); c.genome_screened = i.u64();
        for (ScoreHistogram* h : {&c.locus_score, &c.off_score, &c.genome_score}) {
            h->bin = i.vec<std::int32_t>();
            h->count = i.vec<std::uint64_t>();
            if (h->bin.size() != h->count.size() || !std::is_sorted(h->bin.begin(), h->bin.end()))
                throw std::runtime_error("recruit: " + path + " is corrupt");
        }
        c.genomes.resize(i.length(8));
        for (auto& g : c.genomes) { g.genome = i.str(); g.fragments = i.u64(); g.screened = i.u64(); g.with_private = i.u64(); }
    }
    x.calibration_genomes_given.resize(i.length(8));
    for (auto& g : x.calibration_genomes_given) g = i.str();
    x.calibration_window_list.resize(i.length(8));
    for (auto& w : x.calibration_window_list) { w.genome = i.str(); w.contig = i.str(); w.start = i.u64(); w.end = i.u64(); w.hits = i.u64(); }
    x.calibration_genomes.resize(i.length(8));
    for (auto& g : x.calibration_genomes) {
        g.genome = i.str();
        for (std::uint64_t* v : {&g.windows, &g.bases, &g.homolog_kmers, &g.relabelled, &g.unlabelled}) *v = i.u64();
    }
    const std::uint64_t sum = i.h;
    if (i.raw64() != sum) throw std::runtime_error("recruit: " + path + " is corrupt (its checksum does not match)");
    const std::size_t n = x.codes.size();
    bool ok = x.cls.size() == n && x.outside.size() == n && x.mult.size() == n && x.backbone.size() == n && x.nodeset.size() == n &&
              key == x.key() && !x.nodeset_first.empty() && x.nodeset_first.back() == x.nodeset_nodes.size() &&
              std::is_sorted(x.codes.begin(), x.codes.end()) && std::is_sorted(x.nodeset_first.begin(), x.nodeset_first.end());
    for (std::size_t k = 0; ok && k < n; ++k) ok = x.nodeset[k] + 1 < x.nodeset_first.size() || x.cls[k] != kPrivate;
    for (std::size_t k = 0; ok && k < x.nodeset_nodes.size(); ++k) ok = x.nodeset_nodes[k] < x.nodes.size();
    if (!ok) throw std::runtime_error("recruit: " + path + " is corrupt");
    return x;
}

void write_report(const LocusIndex& x, const std::string& path) {
    const std::string tmp = path + ".tmp." + std::to_string(::getpid());
    std::ofstream r(tmp);
    std::size_t np = 0, ns = 0, nr = 0, bb = 0, multi = 0;
    for (std::size_t i = 0; i < x.codes.size(); ++i) {
        (x.cls[i] == kPrivate ? np : x.cls[i] == kShared ? ns : nr)++;
        bb += x.backbone[i];
        multi += x.mult[i] > 1;
    }
    std::uint64_t pieces = 0, piece_bp = 0;
    for (const auto& sp : x.spans) if (sp.piece) { ++pieces; piece_bp += sp.end - sp.start; }
    auto pct = [&](std::size_t v) { char b[16]; std::snprintf(b, sizeof b, "%.1f", 100.0 * static_cast<double>(v) / static_cast<double>(std::max<std::size_t>(1, x.codes.size()))); return std::string(b); };
    r << "section\tkey\tvalue\n";
    r << "locus\tname\t" << x.name << "\nlocus\tgraph\t" << x.graph_path << "\nlocus\treference_path\t" << x.reference_path
      << "\nlocus\tinterval\t" << x.contig << ':' << x.start << '-' << x.end << "\nlocus\tpaths\t" << x.paths
      << "\nlocus\tbackground\t" << x.background << " (" << x.background_kind << ", content " << std::hex << x.content << std::dec << ")"
      << "\nlocus\thaplotypes\t" << x.haplotypes << " in the background (shared: up to " << x.informative_max * x.haplotypes << " copies outside)"
      << "\nlocus\tkmers\t" << x.codes.size() << "\nlocus\tprivate\t" << np << " (" << pct(np) << "%)"
      << "\nlocus\tshared\t" << ns << " (" << pct(ns) << "%)" << "\nlocus\trepeat\t" << nr << " (" << pct(nr) << "%)"
      << "\nlocus\tbackbone\t" << bb << "\nlocus\tmulti_copy\t" << multi << "\nlocus\tnode_sets\t" << (x.nodeset_first.empty() ? 0 : x.nodeset_first.size() - 1)
      << "\nlocus\tdust_demoted\t" << x.dust_demoted << "\nlocus\tspan_demoted\t" << x.span_demoted << "\nlocus\tnegative_residuals\t" << x.negative
      << "\nlocus\tpieces\t" << pieces << " (" << piece_bp << " bp)"
      << "\nlocus\tdecoys\t" << x.decoy.size() << "\nlocus\tdecoy_windows\t" << x.decoy_windows << " (" << x.decoy_bases << " bp)"
      << "\nlocus\tscreen\t" << x.screen.size();
    std::string held;
    for (const auto& g : x.calibration_genomes) held += (held.empty() ? "" : ",") + g.genome;
    r << "\nlocus\tcalibration_genomes\t"
      << (!x.calibrated ? std::string("none (no calibration)")
          : !held.empty() ? held + (x.calibration_genomes_given.empty() ? " (chosen by the homolog k-mers each adds)" : " (given)")
                          : std::string("none (off-target reads gain private sites by sequencing errors only)"))
      << '\n';
    r << "span\tgenome\tcontig\tstart\tend\tstrand\tsubtracted\tanchored_start\tanchored_end\tlocus_start\tlocus_end\n";
    for (const auto& s : x.spans)
        if (!s.piece)
            r << "span\t" << s.genome << '\t' << s.contig << '\t' << s.start << '\t' << s.end << '\t' << (s.forward ? '+' : '-') << '\t' << s.subtracted
              << '\t' << s.anchored_start << '\t' << s.anchored_end << '\t' << s.locus_start << '\t' << s.locus_end << '\n';
    r << "piece\tgenome\tcontig\tstart\tend\tsubtracted\n";
    for (const auto& s : x.spans)
        if (s.piece) r << "piece\t" << s.genome << '\t' << s.contig << '\t' << s.start << '\t' << s.end << '\t' << s.subtracted << '\n';
    r << "decoy_window\tgenome\tcontig\tstart\tend\thits\n";
    for (const auto& w : x.decoy_window_list)
        r << "decoy_window\t" << w.genome << '\t' << w.contig << '\t' << w.start << '\t' << w.end << '\t' << w.hits << '\n';
    r << "ceiling\tread_length\tp3_min\tp3_median\tsites2_min\tsites2_median\tsingle2_min\tsingle2_median\n";
    for (const auto& c : x.ceilings)
        r << "ceiling\t" << c.length << '\t' << c.p3_min << '\t' << c.p3_median << '\t' << c.sites2_min << '\t' << c.sites2_median << '\t'
          << c.single2_min << '\t' << c.single2_median << '\n';
    // per read class and private-site cut: held-out locus fragments certain (recall), off-target ones certain
    // (false certain) of those screened, and that share's 95% upper bound
    r << "calibration\tclass\terror_rate\tread_length\tinsert\tcut\tlocus\tlocus_certain\trecall\tflank_screened\tflank_certain\tdecoy_screened\tdecoy_certain\tfalse_certain_bound\n";
    for (const auto& c : x.calibration)
        for (std::uint32_t m = 1; m <= kMaxCut; ++m) {
            const std::uint64_t fc = c.flank_certain[m] + c.decoy_certain[m];
            char b[64];
            std::snprintf(b, sizeof b, "%.4f\t%.3g", c.locus ? static_cast<double>(c.locus_certain[m]) / static_cast<double>(c.locus) : 0.0,
                          upper_bound(fc, c.flank_screened + c.decoy_screened));
            const std::string rb(b);
            r << "calibration\t" << c.read_class << '\t' << c.error_rate << '\t' << c.read_length << '\t' << c.insert << '\t' << m << '\t' << c.locus << '\t'
              << c.locus_certain[m] << '\t' << rb.substr(0, rb.find('\t')) << '\t' << c.flank_screened << '\t' << c.flank_certain[m] << '\t'
              << c.decoy_screened << '\t' << c.decoy_certain[m] << '\t' << rb.substr(rb.find('\t') + 1) << '\n';
        }
    // the score: the genomes held out, the windows read from those that are not decoy genomes, the site rates per
    // read class, each genome's off-target fragments, and the thresholds a few rates give
    r << "calibration_genome\tgenome\twindows\tbases\thomolog_kmers\trelabelled_private\tunlabelled_decoys\n";
    for (const auto& g : x.calibration_genomes)
        r << "calibration_genome\t" << g.genome << '\t' << g.windows << '\t' << g.bases << '\t' << g.homolog_kmers << '\t' << g.relabelled << '\t'
          << g.unlabelled << '\n';
    r << "calibration_window\tgenome\tcontig\tstart\tend\thits\n";
    for (const auto& w : x.calibration_window_list)
        r << "calibration_window\t" << w.genome << '\t' << w.contig << '\t' << w.start << '\t' << w.end << '\t' << w.hits << '\n';
    r << "score_rates\tclass\tr_in\tr_out\tq_in\tq_out\toff_screened\tgenome_screened\n";
    for (const auto& c : x.calibration) {
        char b[256];
        std::snprintf(b, sizeof b, "%.6g\t%.6g\t%.6g\t%.6g", c.r_in, c.r_out, c.q_in, c.q_out);
        r << "score_rates\t" << c.read_class << '\t' << b << '\t' << c.flank_screened + c.decoy_screened << '\t' << c.genome_screened << '\n';
    }
    r << "score_genome\tclass\tgenome\tfragments\tscreened\twith_private\n";
    for (const auto& c : x.calibration)
        for (const auto& g : c.genomes)
            r << "score_genome\t" << c.read_class << '\t' << g.genome << '\t' << g.fragments << '\t' << g.screened << '\t' << g.with_private << '\n';
    r << "score_tier\tclass\trate\tthreshold\tlocus_reached\trecall\toff_reached\toff_screened\tbound\n";
    for (const auto& c : x.calibration) {
        if (!has_score(c)) continue;
        for (const double rate : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5}) {
            const Tier t = score_tier(c, rate);
            char th[32], b[256];
            if (t.bin == kEveryBin) std::snprintf(th, sizeof th, "any");
            else std::snprintf(th, sizeof th, "%.4f", t.bin * kScoreBin);
            std::snprintf(b, sizeof b, "%g\t%s\t%llu\t%.4f\t%llu\t%llu\t%.3g", rate, th, static_cast<unsigned long long>(t.locus),
                          t.locus_total ? static_cast<double>(t.locus) / static_cast<double>(t.locus_total) : 0.0, static_cast<unsigned long long>(t.off),
                          static_cast<unsigned long long>(t.off_screened), t.bound());
            r << "score_tier\t" << c.read_class << '\t' << b << '\n';
        }
    }
    r.close();
    if (!r) throw std::runtime_error("recruit: cannot write " + tmp);
    fs::rename(tmp, path);
}

}  // namespace panvar::recruit
