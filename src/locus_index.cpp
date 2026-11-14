#include "panvar/locus_index.hpp"

#include "panvar/anchors.hpp"
#include "panvar/background.hpp"
#include "panvar/background_graph.hpp"
#include "panvar/fasta_index.hpp"
#include "panvar/flat_map.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/gz_reader.hpp"
#include "panvar/kmer.hpp"
#include "panvar/recruit_calibration.hpp"
#include "panvar/ref_path.hpp"

#include <algorithm>
#include <cmath>
#include <atomic>
#include <cstdio>
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
constexpr std::uint64_t kVersion = 13;   // 2: decoy windows stop at the locus; 3: clusters split at it; 4: S per haplotype; 5: calibration; 6: its decoy pairs cross window edges; 7: locus pieces off the located copies; 8: and past estimated ends; 9: and everywhere for a locus with array pieces outside; 10: stretches reach contig ends; 11: pieces tied to the locus, short contigs whole, no decoys in pieces; 12: transcript k-mers; 13: those counted in the transcriptome and the run's other loci
const std::size_t kLengths[] = {100, 150, 250, 1000, 10000};
constexpr std::size_t kWindowStep = 25;
// A stretch of a genome is a piece of the locus when, over every window of this many bases in it, this
// share of the k-mer starts are locus k-mers (a held-out haplotype's new alleles, about one per kb, leave
// about 97%) and it is tied to the locus (see the piece scan): a close paralog, or one whose differences
// cluster, can be as dense.
constexpr std::uint64_t kPieceWindow = 5000;
constexpr double kPieceDensity = 0.95;
// A locus is looked for in every genome when this share of its k-mers are tandem (this many copies or more
// within a haplotype) and left with a few copies outside the located ones (at most this many times their own
// copy number): the mark of array pieces next to copies that look whole. A locus without an array has a few
// hundredths of a percent.
constexpr std::uint16_t kTandemCopies = 3;
constexpr std::uint64_t kTandemOutside = 10;
constexpr double kTandemShare = 0.005;

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

// The genes a locus's transcripts are of: their ids (Ensembl-style ones without the version) and their names.
struct OwnGenes {
    std::set<std::string> ids, names;
};

struct LocusWork {
    LocusIndex idx;
    OwnGenes genes;                                    // the genes its transcripts are of
    std::vector<std::uint64_t> tx_count, tx_sub, cross;   // per k-mer (transcript-only ones): in the transcriptome,
                                                          // in its own genes' records there, in other loci
    std::vector<std::uint64_t> counts;                 // background counts per k-mer
    std::vector<std::uint64_t> sub, sub_trim;          // own occurrences in the extended / trimmed spans
    std::map<std::string, std::vector<anchors::Span>> spans;   // per genome: the located copies, then the pieces
                                                               // found off them (the locus, both, from there on)
};

// The table beside a transcripts FASTA: <prefix>.hst.tsv (or .hst.tsv.gz) for <prefix>.hst.fa[.gz]; empty
// when there is none.
std::string transcripts_table(const std::string& fa) {
    std::string p = fa;
    if (p.size() > 3 && p.compare(p.size() - 3, 3, ".gz") == 0) p.resize(p.size() - 3);
    const std::string tail = ".hst.fa";
    if (p.size() <= tail.size() || p.compare(p.size() - tail.size(), tail.size(), tail) != 0) return {};
    const std::string base = p.substr(0, p.size() - tail.size()) + ".hst.tsv";
    for (const std::string& c : {base, base + ".gz"})
        if (fs::is_regular_file(c)) return c;
    return {};
}

// The transcripts' checksum: the FASTA's, and its table's where there is one, as nodes are placed through it.
std::string transcripts_checksum_of(const std::string& fa) {
    const std::string t = transcripts_table(fa);
    return file_checksum(fa) + (t.empty() ? std::string() : "+" + file_checksum(t));
}

// A gene id as it is compared: an Ensembl-style id (ENS, species letters, G, digits) loses its version
// (ENSG00000000001.2), as annotations and FASTA headers differ in it. Any other id or name is kept whole: a
// suffix of a symbol or of a clone-based name tells genes apart.
std::string gene_id_key(const std::string& id) {
    const std::size_t dot = id.rfind('.');
    if (dot == std::string::npos || dot + 1 >= id.size() || id.compare(0, 3, "ENS") != 0 ||
        id.find_first_not_of("0123456789", dot + 1) != std::string::npos)
        return id;
    std::size_t i = 3;
    while (i < dot && std::isupper(static_cast<unsigned char>(id[i]))) ++i;
    if (i == 3 || id[i - 1] != 'G' || i == dot || id.find_first_not_of("0123456789", i) != dot) return id;
    return id.substr(0, dot);
}

bool ensembl_gene_id(const std::string& s) { return gene_id_key(s + ".1") == s; }

// Whether a transcript FASTA header is of one of a locus's own genes. The header's gene id decides where it
// has one and the locus knows its genes' ids: a record of another gene id is another gene, whatever its
// symbol. The symbol decides only without an id on either side. Read from the header: of a first word split
// at |, the second field as the gene id and the sixth (GENCODE) or third as the symbol; of the words after it,
// key=value or key:value with keys gene_id, gene (an id when Ensembl-style, else a symbol), gene_symbol and
// gene_name, up to a free-text description.
bool header_is_own(const std::string& header, const OwnGenes& own) {
    std::vector<std::string> ids, symbols;
    std::size_t a = 1, b = std::min(header.find_first_of(" \t", a), header.size());
    {
        std::vector<std::string> f;
        for (std::size_t x = a; x <= b;) {
            const std::size_t y = std::min(header.find('|', x), b);
            f.push_back(header.substr(x, y - x));
            x = y + 1;
        }
        if (f.size() >= 2 && !f[1].empty()) ids.push_back(f[1]);
        if (f.size() >= 6 && !f[5].empty()) symbols.push_back(f[5]);
        else if (f.size() == 3 && !f[2].empty()) symbols.push_back(f[2]);
    }
    while (b < header.size()) {
        a = header.find_first_not_of(" \t", b);
        if (a == std::string::npos) break;
        b = std::min(header.find_first_of(" \t", a), header.size());
        const std::size_t sep = header.find_first_of("=:", a);
        if (sep == std::string::npos || sep >= b) continue;
        const std::string key = header.substr(a, sep - a), value = header.substr(sep + 1, b - sep - 1);
        if (key == "description") break;
        if (value.empty()) continue;
        if (key == "gene_id" || (key == "gene" && ensembl_gene_id(gene_id_key(value)))) ids.push_back(value);
        else if (key == "gene" || key == "gene_symbol" || key == "gene_name") symbols.push_back(value);
    }
    if (!ids.empty() && !own.ids.empty())
        return std::any_of(ids.begin(), ids.end(), [&](const std::string& x) { return own.ids.count(gene_id_key(x)) != 0; });
    return std::any_of(symbols.begin(), symbols.end(), [&](const std::string& x) { return own.names.count(x) != 0; });
}

// A transcript's sequence, and the graph node of each of its bases (-1 where unknown).
using TranscriptFn = std::function<void(const std::string& seq, const std::function<std::int32_t(std::size_t)>& node_at)>;

// The transcripts of a panvar transcripts FASTA (plain or gzip), one at a time. Every written transcript is
// taken but those whose splice site differs from the reference's (the haplotype may not splice there, so
// their junctions are guesses) unless `splice_changed`, or only status=ok with `ok_only`; a record without a
// status is taken. `by_status` counts the records taken and all, per status. Their bases are placed on the
// graph's nodes through the table beside the FASTA (the haplotype's path and each exon's first base and
// length on it), where it exists; `short_rows` counts its rows cut short. A sequence equal to one already
// passed on is counted as taken but not passed on again: most transcripts are the same on many haplotypes.
void for_each_transcript(const std::string& fa, const Graph& g, const std::vector<std::string>& node_names, bool ok_only,
                         bool splice_changed, std::uint64_t& used, std::uint64_t& skipped,
                         std::map<std::string, std::array<std::uint64_t, 2>>& by_status, OwnGenes& genes,
                         std::uint64_t& short_rows, const TranscriptFn& fn) {
    // the exons of each record on its path, and its gene's id, from the table
    struct Place { std::string path, gene_id, gene_name; bool forward = true; std::vector<std::uint64_t> starts, lengths; };
    std::unordered_map<std::string, Place> place;
    const std::string table = transcripts_table(fa);
    if (!table.empty()) {
        GzLineReader in(table);
        if (!in.ok()) throw std::runtime_error("recruit: cannot read " + table);
        std::string line;
        in.getline(line);
        std::map<std::string, std::size_t> col;
        {
            std::stringstream ss(line);
            std::size_t k = 0;
            for (std::string f; std::getline(ss, f, '\t');) col[f] = k++;
        }
        for (const char* need : {"transcript_id", "gene_name", "haplotype", "copy", "strand", "orientation", "exon_starts", "exon_lengths"})
            if (!col.count(need)) throw std::runtime_error("recruit: " + table + " has no " + need + " column");
        const bool has_source = col.count("source") && col.count("record");
        while (in.getline(line)) {
            std::vector<std::string> f;
            std::stringstream ss(line);
            for (std::string x; std::getline(ss, x, '\t');) f.push_back(x);
            if (f.size() < col.size()) { ++short_rows; continue; }   // a cut row: its transcript is placed on no node
            if (f[col["exon_starts"]] == ".") continue;
            std::string key = f[col["transcript_id"]] + "|" + f[col["gene_name"]] + "|" + f[col["haplotype"]] + "|" + f[col["copy"]];
            if (has_source && f[col["source"]] == "annotated") key += "|annotated:" + f[col["record"]];
            Place p;
            p.path = f[col["haplotype"]];
            if (col.count("gene_id")) p.gene_id = f[col["gene_id"]];
            p.gene_name = f[col["gene_name"]];
            p.forward = (f[col["strand"]] == "+") == (f[col["orientation"]] != "-");
            for (const auto& [name, v] : {std::make_pair(std::string("exon_starts"), &p.starts), std::make_pair(std::string("exon_lengths"), &p.lengths)}) {
                std::stringstream xs(f[col[name]]);
                for (std::string x; std::getline(xs, x, ',');) v->push_back(std::stoull(x));
            }
            place.emplace(key, std::move(p));
        }
    }
    std::unordered_map<std::string, std::size_t> path_of, node_index;
    for (std::size_t p = 0; p < g.paths.size(); ++p) path_of[g.paths[p].name] = p;
    for (std::size_t i = 0; i < node_names.size(); ++i) node_index[node_names[i]] = i;
    std::map<std::size_t, std::vector<std::uint64_t>> step_end;   // per path used, each step's end
    std::unordered_set<std::uint64_t> seen;                        // sequences passed on, by a 64-bit hash
    GzLineReader in(fa);
    if (!in.ok()) throw std::runtime_error("recruit: cannot read " + fa);
    std::string line, name, status, seq;
    auto flush = [&] {
        if (name.empty()) return;
        static const std::set<std::string> written = {"ok", "boundary_moved", "clipped_at_edge", "intron_lost"};
        const bool take = status.empty() || (ok_only ? status == "ok" : written.count(status) != 0 || (splice_changed && status == "splice_site_changed"));
        auto& tally = by_status[status.empty() ? "no_status" : status];
        ++tally[1];
        if (!take || seq.size() < kK) { ++skipped; return; }
        ++used;
        ++tally[0];
        {   // the record's gene: its id and name in the table; without one, the name in the record's own name
            // (transcript|gene|haplotype|copy)
            const auto usable = [](const std::string& v) { return !v.empty() && v != "."; };
            const auto at = place.find(name);
            const std::size_t a = name.find('|'), b = a == std::string::npos ? a : name.find('|', a + 1);
            if (at != place.end()) {
                if (usable(at->second.gene_id)) genes.ids.insert(gene_id_key(at->second.gene_id));
                if (usable(at->second.gene_name)) genes.names.insert(at->second.gene_name);
            } else if (b != std::string::npos && usable(name.substr(a + 1, b - a - 1))) {
                genes.names.insert(name.substr(a + 1, b - a - 1));
            }
        }
        for (char& c : seq) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        // a 64-bit hash: two different transcripts colliding would drop one, odds ~1e-19 per pair here
        const std::uint64_t h = std::hash<std::string>{}(seq) * 0x9e3779b97f4a7c15ULL ^ seq.size();
        if (!seen.insert(h).second) return;
        const auto pl = place.find(name);
        const auto pi = pl == place.end() ? path_of.end() : path_of.find(pl->second.path);
        if (pi == path_of.end()) {
            fn(seq, [](std::size_t) { return std::int32_t(-1); });
            return;
        }
        const PathRecord& path = g.paths[pi->second];
        auto& ends = step_end[pi->second];
        if (ends.empty()) {
            std::uint64_t pos = 0;
            for (const auto& st : path.steps) { pos += g.nodes.at(st.node_id).sequence.size(); ends.push_back(pos); }
        }
        const Place& p = pl->second;
        std::vector<std::uint64_t> cum(1, 0);   // the transcript offset each exon starts at, as written
        for (const std::uint64_t l : p.lengths) cum.push_back(cum.back() + l);
        fn(seq, [&](std::size_t t) -> std::int32_t {
            if (t >= cum.back()) return -1;
            const std::size_t e = static_cast<std::size_t>(std::upper_bound(cum.begin(), cum.end(), t) - cum.begin()) - 1;
            const std::uint64_t o = t - cum[e];
            const std::uint64_t pos = p.forward ? p.starts[e] + o : p.starts[e] + p.lengths[e] - 1 - o;
            const std::size_t si = static_cast<std::size_t>(std::upper_bound(ends.begin(), ends.end(), pos) - ends.begin());
            return si < path.steps.size() ? static_cast<std::int32_t>(node_index.at(path.steps[si].node_id)) : -1;
        });
    };
    while (in.getline(line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            flush();
            seq.clear();
            status.clear();
            const std::size_t sp = line.find_first_of(" \t");
            name = line.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
            const std::size_t st = line.find(" status=");
            if (st != std::string::npos) status = line.substr(st + 8, line.find(' ', st + 8) - st - 8);
            continue;
        }
        seq += line;
    }
    flush();
}

// Every k-mer of every path: nodes it starts on, largest copy number within a path, paths holding it. With
// transcripts, their k-mers no path holds join as transcript-only: one copy, in no path, on the node of
// their first base where it is known.
void collect_kmers(const Graph& g, LocusIndex& idx, std::vector<std::uint32_t>& present,
                   const std::function<void(const TranscriptFn&)>& each_transcript = {}) {
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
    std::unordered_set<std::uint64_t> from_transcripts;
    const bool with_transcripts = static_cast<bool>(each_transcript);
    if (with_transcripts)
        each_transcript([&](const std::string& seq, const std::function<std::int32_t(std::size_t)>& node_at) {
            kmer::for_each_canonical_kmer(seq.data(), seq.size(), kK, [&](std::uint64_t c, std::size_t start, bool) {
                if (agg.count(c) && !from_transcripts.count(c)) return;   // a genomic k-mer
                Agg& a = agg[c];
                if (from_transcripts.insert(c).second) a.mult = 1;
                if (a.nodes.empty()) {
                    const std::int32_t n = node_at(start);
                    if (n >= 0) a.nodes.push_back(static_cast<std::uint32_t>(n));
                }
            });
        });
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
        if (with_transcripts) idx.transcript_only.push_back(from_transcripts.count(c) != 0);
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
           coordinate_reference + "|" + std::to_string(informative_max) + "|" + std::to_string(pad) + "|" + std::to_string(dust) + "|" + std::to_string(decoy_min_hits) + "|" + d + (calibrated ? "|c1" : "|c0") +
           (transcripts_checksum.empty() ? "" : "|t" + transcripts_checksum + (transcripts_ok_only ? "o" : transcripts_splice_changed ? "s" : "a") + "|x" + hex64(transcriptome));
}

std::string expected_key(const BuildOptions& o, const std::string& gfa, const std::string& transcripts) {
    const bg::Manifest m = bg::read_manifest(o.background);
    LocusIndex k;
    k.graph_checksum = file_checksum(gfa);
    if (!transcripts.empty()) {
        k.transcripts_checksum = transcripts_checksum_of(transcripts);
        k.transcripts_ok_only = o.transcripts_ok_only;
        k.transcripts_splice_changed = o.transcripts_splice_changed;
        k.transcriptome = bg::transcriptome_id(m);
    }
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
    k.calibrated = o.calibrate;
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

bool transcript_only(const LocusIndex& x, std::size_t i) { return !x.transcript_only.empty() && x.transcript_only[i]; }

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
    std::string kind;
    const auto counts = open_counts(o.background, kind, o.log);
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
    parallel_for(o.loci.size(), threads, [&](std::size_t li) {
        LocusWork& w = work[li];
        LocusIndex& x = w.idx;
        x.graph_path = fs::absolute(o.loci[li]).string();
        x.name = locus_name_of(o.loci[li]);
        x.graph_checksum = file_checksum(o.loci[li]);
        graphs[li] = std::make_unique<Graph>(parse_gfa(o.loci[li]));
        const Graph& g = *graphs[li];
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
        x.haplotypes = haplotypes;
        x.pad = o.pad;
        x.dust = o.dust;
        x.decoy_min_hits = o.decoy_min_hits;
        x.decoy_genomes = decoys;
        x.calibrated = o.calibrate;
        std::vector<std::uint32_t> present;
        const std::string tx = li < o.transcripts.size() ? o.transcripts[li] : std::string();
        if (tx.empty()) {
            collect_kmers(g, x, present);
        } else {
            std::vector<std::string> node_names;   // as collect_kmers orders them
            for (const auto& [id, n] : g.nodes) node_names.push_back(id);
            const bool numeric = std::all_of(node_names.begin(), node_names.end(), [](const std::string& v) {
                return !v.empty() && v.size() < 19 && v.find_first_not_of("0123456789") == std::string::npos;
            });
            std::sort(node_names.begin(), node_names.end(), [&](const std::string& a, const std::string& b) {
                return numeric ? std::stoull(a) < std::stoull(b) : a < b;
            });
            x.transcripts = fs::absolute(tx).string();
            x.transcripts_checksum = transcripts_checksum_of(tx);
            x.transcripts_ok_only = o.transcripts_ok_only;
            x.transcripts_splice_changed = o.transcripts_splice_changed;
            x.transcriptome = bg::transcriptome_id(m);
            std::map<std::string, std::array<std::uint64_t, 2>> by_status;
            std::uint64_t short_rows = 0;
            collect_kmers(g, x, present, [&](const TranscriptFn& fn) {
                for_each_transcript(tx, g, node_names, o.transcripts_ok_only, o.transcripts_splice_changed, x.transcripts_used,
                                    x.transcripts_skipped, by_status, w.genes, short_rows, fn);
            });
            if (short_rows > 0)
                log("warning: locus: " + x.name + "; " + transcripts_table(tx) + ": " + std::to_string(short_rows) +
                    " row(s) with fewer fields than the header, skipped: their transcripts' k-mers are placed on no node");
            for (const auto& [status, n] : by_status)
                x.transcripts_by_status += (x.transcripts_by_status.empty() ? "" : ", ") + status + " " + std::to_string(n[0]) + "/" + std::to_string(n[1]);
            log("locus: " + x.name + "; transcripts taken of all, by status: " + x.transcripts_by_status);
            const std::size_t n = static_cast<std::size_t>(std::count(x.transcript_only.begin(), x.transcript_only.end(), 1));
            log("locus: " + x.name + "; transcripts: " + std::to_string(x.transcripts_used) + " taken, " + std::to_string(x.transcripts_skipped) +
                " left out; k-mers only in transcripts: " + std::to_string(n));
        }
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
    });

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
        check_unchanged(ge);
        std::vector<std::vector<bg::graph::Haplotypes::Piece>> pieces(iv.size());   // a FASTA span is one piece
        if (is_gfa_file(ge.path)) {
            const bg::graph::Haplotypes h(ge.path);
            for (std::size_t k = 0; k < iv.size(); ++k) pieces[k] = h.fetch(iv[k].contig, iv[k].start, iv[k].end);
        } else {
            if (!fasta::IndexedFasta::indexed(ge.path)) log("note: genome: " + ge.path + "; no .fai (and .gzi): read in full to cut its spans");
            std::vector<std::string> seqs = fasta::fetch_intervals(ge.path, iv);
            for (std::size_t k = 0; k < iv.size(); ++k) pieces[k].push_back({iv[k].start, std::move(seqs[k])});
        }
        for (std::size_t k = 0; k < iv.size(); ++k) {
            LocusWork& w = work[owner[k].first];
            const anchors::Span& sp = w.spans[gname][owner[k].second];
            std::uint64_t found = 0;
            for (const auto& piece : pieces[k])
                kmer::for_each_canonical_kmer(piece.seq.data(), piece.seq.size(), kK, [&](std::uint64_t c, std::size_t at, bool) {
                    const std::int64_t i = w.idx.find(c);
                    if (i < 0) return;
                    ++w.sub[i];
                    ++found;
                    const std::uint64_t pos = piece.start + at;
                    if (pos >= sp.inner_start && pos + kK <= sp.inner_end) ++w.sub_trim[i];
                });
            w.idx.spans.push_back({gname, sp.contig, sp.start, sp.end, sp.forward, found});
        }
    }

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
            // the share is of the haplotypes' k-mers: a transcript's k-mers must not move the trigger
            const std::size_t genomic = w.idx.codes.size() -
                static_cast<std::size_t>(std::count(w.idx.transcript_only.begin(), w.idx.transcript_only.end(), 1));
            if (static_cast<double>(tandem_out) >= kTandemShare * static_cast<double>(genomic)) {
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
                for (std::size_t i = 0; i < work[li].idx.codes.size(); ++i) {
                    if (transcript_only(work[li].idx, i)) continue;   // pieces are found from the genomic walks
                    const std::uint64_t c = work[li].idx.codes[i];
                    bool existing = false;
                    std::uint32_t& id = lk.slot(c, existing);
                    if (!existing) { id = static_cast<std::uint32_t>(lk_owners.size()); lk_owners.emplace_back(); }
                    lk_owners[id].push_back(static_cast<std::uint32_t>(li));
                }
            }
        }
        const std::size_t need = static_cast<std::size_t>(std::ceil(kPieceDensity * static_cast<double>(kPieceWindow)));
        // what a genome's scan finds, applied in genome order once all are scanned (genomes run in parallel)
        struct Found {
            std::vector<std::pair<std::uint32_t, std::uint32_t>> off, ends;   // (locus, k-mer): off the copies / in an uncertain end
            std::vector<std::pair<std::size_t, GenomeSpan>> pieces;
        };
        auto scan = [&](const std::string& gname, const std::string& name, const std::string& seq, const std::vector<char>& here,
                        std::size_t inner, Found& f) {
            if (seq.size() < kK) return;
            const bool whole_only = seq.size() < kPieceWindow + kK - 1;   // no window fits
            // each locus's k-mer starts on the contig, in order: ranges scanned in parallel, joined in order
            const std::size_t starts = seq.size() - kK + 1;
            const std::size_t parts = std::min<std::size_t>(inner, std::max<std::size_t>(1, starts / 65536));
            std::vector<std::vector<std::vector<std::uint64_t>>> part_hits(parts, std::vector<std::vector<std::uint64_t>>(work.size()));
            auto scan_part = [&](std::size_t p) {
                const std::size_t a = starts * p / parts, b = starts * (p + 1) / parts;
                if (b <= a) return;
                kmer::for_each_canonical_kmer(seq.data() + a, b - a + kK - 1, kK, [&](std::uint64_t c, std::size_t off, bool) {
                    const std::uint32_t* id = lk.find(c);
                    if (id == nullptr) return;
                    for (const std::uint32_t li : lk_owners[*id]) if (here[li]) part_hits[p][li].push_back(a + off);
                });
            };
            std::vector<std::thread> pool;
            for (std::size_t p = 1; p < parts; ++p) pool.emplace_back(scan_part, p);
            scan_part(0);
            for (auto& t : pool) t.join();
            for (std::size_t li = 0; li < work.size(); ++li) {
                if (!here[li]) continue;
                std::vector<std::uint64_t> hits;
                for (std::size_t p = 0; p < parts; ++p) hits.insert(hits.end(), part_hits[p][li].begin(), part_hits[p][li].end());
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
                        if (first.first > 0 && first.first < kPieceWindow) {
                            const std::size_t n = static_cast<std::size_t>(std::lower_bound(hits.begin(), hits.end(), first.first) - hits.begin());
                            if (2 * n >= first.first) first.first = 0;
                        }
                        auto& last = dense.back();
                        if (last.second < seq.size() && seq.size() - last.second < kPieceWindow) {
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
                for (const auto& stretch : dense) {
                    const std::uint64_t a = stretch.first, b = stretch.second;
                    // only a stretch tied to the locus is a piece: on a contig of a located copy, or reaching a
                    // contig end (an array broken across contigs, a stub, a contig of its own). One with other
                    // sequence on both sides, on a contig without a copy, may be a paralog the coordinate
                    // reference lacks: it stays outside (shared, never falsely private).
                    if (own.empty() && a > 0 && b < seq.size()) continue;
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
                    kmer::for_each_canonical_kmer(seq.data() + a, b - a, kK, [&](std::uint64_t c, std::size_t at, bool) {
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
            if (!nd.whole.empty()) {
                GzLineReader in(ge.path);
                if (!in.ok()) throw std::runtime_error("recruit: cannot read " + ge.path);
                std::string line, name, seq;
                while (in.getline(line)) {
                    if (line.empty()) continue;
                    if (line[0] == '>') {
                        if (!name.empty()) scan(gname, name, seq, here_for(name), inner, found[gi]);
                        const auto sp = line.find_first_of(" \t");
                        name = line.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
                        seq.clear();
                        continue;
                    }
                    seq += line;
                }
                if (!name.empty()) scan(gname, name, seq, here_for(name), inner, found[gi]);
            } else {
                std::vector<fasta::Interval> iv;
                for (const auto& c : nd.contigs) iv.push_back({c.first, 0, nd.lengths.at(c.first)});
                const std::vector<std::string> seqs = fasta::fetch_intervals(ge.path, iv);
                for (std::size_t k = 0; k < iv.size(); ++k) scan(gname, iv[k].contig, seqs[k], here_for(iv[k].contig), inner, found[gi]);
            }
            std::lock_guard<std::mutex> lock(log_mutex);
            log(nd.whole.empty() ? "pieces: contigs scanned: " + gname + " (" + std::to_string(nd.contigs.size()) + ")"
                                 : "pieces: genome scanned: " + gname);
        });
        std::vector<std::size_t> scanned_for(work.size(), 0);
        std::vector<std::set<std::string>> piece_genomes(work.size());
        for (std::size_t gi = 0; gi < scan_list.size(); ++gi) {
            const Need& nd = scan_list[gi].second;
            std::vector<char> counted(work.size(), 0);
            for (const std::size_t li : nd.whole) counted[li] = 1;
            for (const auto& c : nd.contigs) for (const std::size_t li : c.second) counted[li] = 1;
            for (std::size_t li = 0; li < work.size(); ++li) scanned_for[li] += counted[li];
            for (const auto& [li, i] : found[gi].off) { ++work[li].sub[i]; ++work[li].sub_trim[i]; }
            for (const auto& [li, i] : found[gi].ends) ++work[li].sub_trim[i];
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

    // transcript-only k-mers elsewhere: in the background's transcriptomes, less the records of the locus's
    // own genes (so a paralog's junction stays counted), and in the other loci of this run (their haplotypes
    // and transcripts: a paralog's alleles the reference transcriptome lacks). Loci overlapping on the
    // reference are one place, not each other's elsewhere.
    {
        std::vector<std::size_t> with;
        for (std::size_t li = 0; li < work.size(); ++li)
            if (!work[li].idx.transcript_only.empty()) with.push_back(li);
        for (const std::size_t li : with) {
            LocusWork& w = work[li];
            w.tx_count.assign(w.idx.codes.size(), 0);
            w.tx_sub.assign(w.idx.codes.size(), 0);
            w.cross.assign(w.idx.codes.size(), 0);
        }
        if (!with.empty())
            for (const auto& ts : m.transcriptomes) {
                const bg::Segment seg((fs::path(o.background) / ts.dir).string());
                for (const std::size_t li : with) {
                    LocusWork& w = work[li];
                    for (std::size_t i = 0; i < w.idx.codes.size(); ++i)
                        if (w.idx.transcript_only[i]) w.tx_count[i] += seg.count(w.idx.codes[i]);
                }
                for (const auto& ge : ts.genomes) {
                    check_unchanged(ge);
                    for (const std::size_t li : with)
                        work[li].idx.transcriptome_names += (work[li].idx.transcriptome_names.empty() ? "" : ",") + ge.name;
                    GzLineReader in(ge.path);
                    if (!in.ok()) throw std::runtime_error("recruit: cannot read " + ge.path);
                    std::string line, seq;
                    std::vector<std::size_t> own;   // the loci whose genes the record is of
                    auto flush = [&] {
                        for (const std::size_t li : own) {
                            LocusWork& w = work[li];
                            ++w.idx.transcriptome_records;
                            kmer::for_each_canonical_kmer(seq.data(), seq.size(), kK, [&](std::uint64_t c, std::size_t, bool) {
                                const std::int64_t i = w.idx.find(c);
                                if (i >= 0 && w.idx.transcript_only[static_cast<std::size_t>(i)]) ++w.tx_sub[static_cast<std::size_t>(i)];
                            });
                        }
                        seq.clear();
                    };
                    while (in.getline(line)) {
                        if (line.empty()) continue;
                        if (line[0] == '>') {
                            flush();
                            own.clear();
                            for (const std::size_t li : with)
                                if (header_is_own(line, work[li].genes)) own.push_back(li);
                        } else if (!own.empty()) {
                            for (const char c : line) seq += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                        }
                    }
                    flush();
                }
            }
        for (const std::size_t li : with) {
            LocusWork& w = work[li];
            LocusIndex& x = w.idx;
            for (std::size_t lj = 0; lj < work.size(); ++lj) {
                const LocusIndex& y = work[lj].idx;
                if (lj == li || (x.contig == y.contig && x.start < y.end && y.start < x.end)) continue;
                std::uint64_t hits = 0;
                for (std::size_t i = 0; i < x.codes.size(); ++i)
                    if (x.transcript_only[i] && y.find(x.codes[i]) >= 0) { ++w.cross[i]; ++hits; }
                if (hits > 0) x.other_loci += (x.other_loci.empty() ? "" : ",") + y.name;
            }
            for (std::size_t i = 0; i < x.codes.size(); ++i) {
                if (!x.transcript_only[i]) continue;
                x.tx_in_transcriptome += w.tx_count[i] > w.tx_sub[i];
                x.tx_in_other_loci += w.cross[i] > 0;
            }
            log("locus: " + x.name + "; transcriptome: " + (x.transcriptome_names.empty() ? std::string("none in the background") : x.transcriptome_names +
                "; records of the locus's genes subtracted: " + std::to_string(x.transcriptome_records)) +
                "; transcript-only k-mers in other transcripts: " + std::to_string(x.tx_in_transcriptome) + "; in other loci of the run: " +
                std::to_string(x.tx_in_other_loci));
            if (!x.transcriptome_names.empty() && x.transcriptome_records == 0)
                log("warning: locus: " + x.name + "; no record of the transcriptome names a gene of its transcripts: the locus's own transcripts are "
                    "counted as elsewhere (the FASTA headers must carry the genes' ids, or their names where they carry no id)");
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
            // a transcript-only k-mer: nothing subtracted from the genomes' count (none holds a transcript), and
            // what other transcripts and loci hold added
            std::uint64_t elsewhere = 0;
            if (transcript_only(x, i)) {
                w.sub[i] = w.sub_trim[i] = 0;
                elsewhere = w.tx_count[i] - std::min(w.tx_sub[i], w.tx_count[i]) + w.cross[i];
            }
            if (w.sub[i] > w.counts[i]) { ++x.negative; w.sub[i] = w.counts[i]; }
            const std::uint64_t out = w.counts[i] - w.sub[i] + elsewhere;
            const std::uint64_t out_trim = w.counts[i] - std::min(w.sub_trim[i], w.counts[i]) + elsewhere;
            x.outside[i] = static_cast<std::uint32_t>(std::min<std::uint64_t>(out, 0xFFFFFFFFULL));
            std::uint8_t c = class_of(out);
            if (c != class_of(out_trim) && c == kPrivate) { c = kShared; ++x.span_demoted; }
            if (c == kPrivate && kmer::dust_score(x.codes[i], kK) >= x.dust) { c = kShared; ++x.dust_demoted; }
            x.cls[i] = c;
        }
        if (x.negative > 0)
            throw std::runtime_error("recruit: " + x.name + ": " + std::to_string(x.negative) +
                                     " k-mers occur more often in the located copies than the background counted: a genome file differs from the one counted");
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
                if (work[li].idx.cls[i] != kShared || transcript_only(work[li].idx, i)) continue;   // decoys from the genomic walks
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
    for (const auto& dname : decoys) {
        const bg::Genome& ge = genomes.at(dname);
        if (is_gfa_file(ge.path)) throw std::runtime_error("recruit: decoy genome " + dname + " is a graph; decoys must be FASTA");
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
                std::vector<std::thread> pool;
                for (std::size_t p = 1; p < parts; ++p) pool.emplace_back(scan_part, p);
                scan_part(0);
                for (auto& t : pool) t.join();
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
        log("decoys: genome scanned: " + dname);
    }
    parallel_for(work.size(), threads, [&](std::size_t li) {
        LocusIndex& x = work[li].idx;
        std::vector<std::uint64_t> cand(candidates[li].begin(), candidates[li].end());
        std::sort(cand.begin(), cand.end());
        const auto cc = counts->count_batch(cand);
        for (std::size_t i = 0; i < cand.size(); ++i)
            if (cc[i] <= static_cast<std::uint64_t>(x.informative_max) * x.haplotypes) x.decoy.push_back(cand[i]);
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
        for (std::size_t li = 0; li < work.size(); ++li) {
            LocusIndex& x = work[li].idx;
            calibrate(x, *graphs[li], off[li], threads);
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
        for (const auto& s : x.spans) { o.str(s.genome); o.str(s.contig); o.u64(s.start); o.u64(s.end); o.u64(s.forward); o.u64(s.subtracted); o.u64(s.piece); }
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
        }
        if (!x.transcripts_checksum.empty()) {   // only with transcripts, so an index without is as before
            o.str(x.transcripts); o.str(x.transcripts_checksum);
            o.u64(x.transcripts_ok_only); o.u64(x.transcripts_used); o.u64(x.transcripts_skipped);
            o.vec(x.transcript_only);
            o.u64(x.transcriptome); o.u64(x.transcriptome_records); o.u64(x.tx_in_transcriptome); o.u64(x.tx_in_other_loci);
            o.str(x.transcriptome_names); o.str(x.other_loci);
            o.u64(x.transcripts_splice_changed); o.str(x.transcripts_by_status);
        }
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
    for (auto& s : x.spans) { s.genome = i.str(); s.contig = i.str(); s.start = i.u64(); s.end = i.u64(); s.forward = i.u64() != 0; s.subtracted = i.u64(); s.piece = i.u64() != 0; }
    x.ceilings = i.vec<Ceiling>();
    x.dust_demoted = i.u64(); x.span_demoted = i.u64(); x.negative = i.u64();
    x.decoy_min_hits = static_cast<std::uint32_t>(i.u64()); x.decoy_windows = i.u64(); x.decoy_bases = i.u64(); x.haplotypes = i.u64();
    x.decoy_window_list.resize(i.u64());
    for (auto& w : x.decoy_window_list) { w.genome = i.str(); w.contig = i.str(); w.start = i.u64(); w.end = i.u64(); w.hits = i.u64(); }
    x.calibrated = i.u64() != 0;
    x.calibration.resize(i.u64());
    for (auto& c : x.calibration) {
        c.read_class = i.str();
        const std::uint64_t ebits = i.u64();
        std::memcpy(&c.error_rate, &ebits, 8);
        for (std::uint64_t* v : {&c.read_length, &c.insert, &c.held_out, &c.locus, &c.flank, &c.decoy, &c.flank_screened, &c.decoy_screened}) *v = i.u64();
        for (auto* a : {&c.locus_certain, &c.flank_certain, &c.decoy_certain})
            for (std::uint64_t& v : *a) v = i.u64();
    }
    if (i.f.peek() != std::char_traits<char>::eof()) {
        x.transcripts = i.str(); x.transcripts_checksum = i.str();
        x.transcripts_ok_only = i.u64() != 0; x.transcripts_used = i.u64(); x.transcripts_skipped = i.u64();
        x.transcript_only = i.vec<std::uint8_t>();
        x.transcriptome = i.u64(); x.transcriptome_records = i.u64(); x.tx_in_transcriptome = i.u64(); x.tx_in_other_loci = i.u64();
        x.transcriptome_names = i.str(); x.other_loci = i.str();
        x.transcripts_splice_changed = i.u64() != 0; x.transcripts_by_status = i.str();
    }
    const std::size_t n = x.codes.size();
    if (x.cls.size() != n || x.outside.size() != n || x.mult.size() != n || x.backbone.size() != n || x.nodeset.size() != n ||
        (!x.transcript_only.empty() && x.transcript_only.size() != n) || key != x.key())
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
      << "\nlocus\tscreen\t" << x.screen.size() << '\n';
    if (!x.transcripts_checksum.empty()) {
        std::size_t t[3] = {0, 0, 0};
        for (std::size_t i = 0; i < x.codes.size(); ++i) if (x.transcript_only[i]) ++t[x.cls[i]];
        r << "locus\ttranscripts\t" << x.transcripts << " (" << x.transcripts_used << " taken, " << x.transcripts_skipped << " left out"
          << (x.transcripts_ok_only ? "; status ok only" : x.transcripts_splice_changed ? "; with changed splice sites" : "") << ")"
          << "\nlocus\ttranscripts_by_status\t" << x.transcripts_by_status << " (taken/all)"
          << "\nlocus\ttranscript_only\t" << t[0] + t[1] + t[2] << " (private " << t[0] << ", shared " << t[1] << ", repeat " << t[2] << ")"
          << "\nlocus\ttranscriptome\t" << (x.transcriptome_names.empty() ? std::string("none in the background")
                                             : x.transcriptome_names + " (records of the locus's genes subtracted: " + std::to_string(x.transcriptome_records) + ")")
          << "\nlocus\ttranscript_only_elsewhere\tin other transcripts: " << x.tx_in_transcriptome << "; in other loci of the run"
          << (x.other_loci.empty() ? "" : " (" + x.other_loci + ")") << ": " << x.tx_in_other_loci << '\n';
    }
    r << "span\tgenome\tcontig\tstart\tend\tstrand\tsubtracted\n";
    for (const auto& s : x.spans)
        if (!s.piece) r << "span\t" << s.genome << '\t' << s.contig << '\t' << s.start << '\t' << s.end << '\t' << (s.forward ? '+' : '-') << '\t' << s.subtracted << '\n';
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
    if (!r) throw std::runtime_error("recruit: cannot write " + path);
}

}  // namespace panvar::recruit
