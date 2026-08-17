#include "panvar/benchmark_command.hpp"

#include "panvar/align.hpp"
#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/output.hpp"
#include "panvar/parallel.hpp"
#include "panvar/ref_path.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panvar {
namespace {

void print_benchmark_help() {
    std::cout
        << "Usage:\n"
        << "  panvar benchmark -i <graph.gfa> (-b <prefix> | -c <bubbles.csv>) \\\n"
        << "                   --variant-nodes <call.variant_nodes.tsv> -r <name> -o <prefix> [options]\n\n"
        << "Round-trip scores of the caller's own output, over the same bubbles, the same truth and the\n"
        << "same QV scale, differing only in what the reconstruction is allowed to use.\n\n"
        << "  graph     Reference walk with the haplotype's OWN steps substituted at every block that\n"
        << "            shares a node with ANY call at that bubble. Deliberately optimistic: sharing a\n"
        << "            node with some call is not the same as matching one, so a called block can\n"
        << "            authorise copying a neighbouring uncalled allele. Read it as a discovery upper\n"
        << "            bound -- can the graph hold this haplotype, and did the caller flag the\n"
        << "            divergent blocks -- never as a genotyping score. No genotype is read.\n"
        << "  called    The same substitution restricted to blocks a SPECIFIC record is attributed to\n"
        << "            AND that reach --min-sv-bp; nothing else of the haplotype's walk is copied. It\n"
        << "            still substitutes the haplotype's TRUE block, not the record's REF/ALT effect,\n"
        << "            so it is the ceiling those records would reach if each reproduced its block\n"
        << "            exactly -- not what they do reproduce. That is `genotype`.\n"
        << "  genotype  (--vcf) Reference sequence with only the edits this haplotype's GT names,\n"
        << "            applied from the VCF alone. A missed carrier keeps reference and a spurious one\n"
        << "            edits sequence that was already correct, so both error directions cost bases.\n"
        << "            This is what a consumer reconstructing a sample from the VCF actually gets.\n\n"
        << "All are aligned to the haplotype's true walk (edlib NW) and scored\n"
        << "QV = -10*log10(max(0.5, delta)/S), combined per haplotype length-weighted (sum delta / sum S)\n"
        << "and binned into the cosigt bands (<17, 17-23, 23-33, >33). The do-nothing baseline (plain\n"
        << "reference) is scored alongside and is the denominator of gap_closed.\n\n"
        << "Separately, a TRUTH EVENT LEDGER classifies what was there to be found, before any call is\n"
        << "consulted: each maximal divergent block between the reference walk and the haplotype's walk\n"
        << "is one event, sized in bp, and classified called / missed / below_threshold. Event size is\n"
        << "a property of the two walks, so it does not move with the alignment.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Passed graph the calls were made on (required)\n"
        << "  -b, --bubble-prefix-in <prefix>  panphorte output prefix (auto-uses <prefix>.bubbles.csv)\n"
        << "  -c, --bubbles-csv-in <path>      panphorte bubbles CSV (required if no prefix)\n"
        << "      --variant-nodes <path>       call's <prefix>.variant_nodes.tsv (required)\n"
        << "      --vcf <path>                 call's <prefix>.region.vcf. Adds the genotype-aware\n"
        << "                                   reconstruction above; without it only the graph and\n"
        << "                                   called reconstructions are scored\n"
        << "  -r, --reference-path <name>      Reference path name/substring, the diff baseline (required)\n"
        << "  -o, --out-prefix <prefix>        Output prefix for the QV tables (required)\n"
        << "      --called-bubbles-only        Score only bubbles carrying >=1 call. Ascertainment-biased\n"
        << "                                   by construction -- a bubble with real variation and no call\n"
        << "                                   is exactly the miss the ledger exists to count, and this\n"
        << "                                   drops it. Default scores every reference-traversed bubble.\n"
        << "      --dup-model <cn|cnbp>        How a copy-number record is reconstructed: lay down CN\n"
        << "                                   copies of RU_LEN in place of REF_CN (cn), or apply\n"
        << "                                   the per-sample CNBP bp delta (cnbp, default). Both tile an\n"
        << "                                   inferred reference span, so DUP reconstruction is a\n"
        << "                                   heuristic on length, not the inserted sequence.\n"
        << "      --min-sv-bp <N>              Event-size threshold for the truth ledger; set it to the\n"
        << "                                   value `call` ran with (default 50)\n"
        << "      --no-truth-events            Do not write the per-event ledger table\n"
        << "      --threads <N>                Worker threads over haplotypes (0 = auto)\n"
        << "  -q, --quiet                      Disable progress logs\n"
        << "  -h, --help                       Show this help\n";
}

// cosigt/locityper QV bands.
const char* qv_band(double qv) {
    if (qv < 17.0) return "<17";
    if (qv < 23.0) return "17-23";
    if (qv < 33.0) return "23-33";
    return ">33";
}

// QV = -10*log10(max(0.5, delta)/S). delta<=0.5 (a perfect/near-perfect match) is floored so the score
// stays finite and grows with region length, matching cosigt's convention.
double qv_from(std::size_t delta, std::size_t aln_len) {
    if (aln_len == 0) return 60.0;  // nothing to align; treat as perfect
    const double num = std::max(0.5, static_cast<double>(delta));
    return -10.0 * std::log10(num / static_cast<double>(aln_len));
}

// Ceiling QV for a perfect (delta<=0.5) reconstruction of total aligned length S: QV_max = 10*log10(2S).
// A short region can never reach the high cosigt bands, so we also report qv / qv_max -- length-fair,
// 1.0 for a perfect reconstruction of any size -- and bin that ratio into quintiles.
double qv_max_from(std::size_t aln_len) {
    if (aln_len == 0) return 60.0;
    return 10.0 * std::log10(2.0 * static_cast<double>(aln_len));
}

const char* qv_ratio_quintile(double ratio) {
    if (ratio < 0.2) return "0.0-0.2";
    if (ratio < 0.4) return "0.2-0.4";
    if (ratio < 0.6) return "0.4-0.6";
    if (ratio < 0.8) return "0.6-0.8";
    return "0.8-1.0";
}

// One emitted call, kept whole. Folding a bubble's calls into a single node set cannot distinguish
// "this truth event matched THAT record" from "this truth event touches something some record also
// touches", and the second is not evidence that anything was called correctly.
struct CalledRecord {
    std::string variant_id;
    std::string svtype;
    std::unordered_set<std::string> nodes;
};

struct CalledVariants {
    std::unordered_set<std::size_t> bubble_ids;                                // bubbles with >=1 call
    std::unordered_map<std::size_t, std::vector<CalledRecord>> records;        // bubble_id -> records
    std::unordered_map<std::size_t, std::unordered_set<std::string>> nodes;    // bubble_id -> union
    std::map<std::size_t, std::set<std::string>> svtypes;                      // bubble_id -> svtypes
};

// `bubble_members` is every node each bubble owns (interior plus both boundaries). A call's nodes are
// the whole basis on which a truth event is attributed to it, so a node that is not in the graph, or
// belongs to a different site, silently turns a called event into a missed one -- the input error and
// the caller failure produce the identical output. Nothing here is skipped or repaired: an input that
// cannot be trusted is refused.
CalledVariants load_variant_nodes(
        const std::string& path,
        const std::unordered_map<std::size_t, std::unordered_set<std::string>>& bubble_members) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Failed to open variant nodes: " + path);
    auto split_tab = [](const std::string& s) {
        std::vector<std::string> f;
        std::string cur;
        for (char c : s) { if (c == '\t') { f.push_back(cur); cur.clear(); } else cur += c; }
        f.push_back(cur);
        return f;
    };
    CalledVariants cv;
    std::string line;
    std::size_t ncol = 0, lineno = 0;
    std::unordered_set<std::string> seen_ids;
    static const char* kHead[] = {"variant_id", "bubble_id", "svtype", "node_ids"};
    while (std::getline(in, line)) {
        ++lineno;
        if (line.empty()) continue;
        const std::vector<std::string> f = split_tab(line);
        if (ncol == 0) {
            // A missing header is not cosmetic: the first data row would be eaten as one, dropping a
            // real call and reporting its event as missed.
            if (f.size() < 4)
                throw std::runtime_error(path + ": expected a header with at least 4 tab-separated "
                                         "columns, found " + std::to_string(f.size()));
            for (int i = 0; i < 4; ++i)
                if (f[static_cast<std::size_t>(i)] != kHead[i])
                    throw std::runtime_error(path + ": column " + std::to_string(i + 1) + " of the header is '" +
                                             f[static_cast<std::size_t>(i)] + "', expected '" + kHead[i] +
                                             "'. This is not a variant_nodes.tsv written by `call`");
            ncol = f.size();
            continue;
        }
        if (f.size() != ncol)
            throw std::runtime_error(path + ":" + std::to_string(lineno) + ": has " +
                                     std::to_string(f.size()) + " fields, the header has " +
                                     std::to_string(ncol));
        if (f[0].empty())
            throw std::runtime_error(path + ":" + std::to_string(lineno) + ": empty variant_id");
        if (!seen_ids.insert(f[0]).second)
            throw std::runtime_error(path + ":" + std::to_string(lineno) + ": duplicate variant_id '" +
                                     f[0] + "'; which record a truth event is attributed to would be undefined");
        std::size_t bid = 0;
        try {
            std::size_t used = 0;
            bid = static_cast<std::size_t>(std::stoull(f[1], &used));
            if (used != f[1].size()) throw std::invalid_argument("trailing");
        } catch (const std::exception&) {
            throw std::runtime_error(path + ":" + std::to_string(lineno) + ": bubble_id '" + f[1] +
                                     "' is not a number");
        }
        const auto mit = bubble_members.find(bid);
        if (mit == bubble_members.end())
            throw std::runtime_error(path + ":" + std::to_string(lineno) + ": bubble id " +
                                     std::to_string(bid) + " is not in the bubbles CSV; these outputs "
                                     "are from different runs");
        CalledRecord rec;
        rec.variant_id = f[0];
        rec.svtype = f[2];
        auto& un = cv.nodes[bid];
        std::string tok;
        auto take = [&]() {
            if (tok.empty()) return;
            if (!mit->second.count(tok))
                throw std::runtime_error(path + ":" + std::to_string(lineno) + ": record '" + f[0] +
                                         "' names node '" + tok + "', which bubble " + std::to_string(bid) +
                                         " does not contain. A stale node turns a called event into a "
                                         "missed one with nothing to show for it");
            rec.nodes.insert(tok);
            un.insert(tok);
            tok.clear();
        };
        for (char c : f[3]) { if (c == ',') take(); else tok += c; }
        take();
        if (rec.nodes.empty())
            throw std::runtime_error(path + ":" + std::to_string(lineno) + ": record '" + f[0] +
                                     "' names no nodes");
        cv.bubble_ids.insert(bid);
        cv.svtypes[bid].insert(f[2]);
        cv.records[bid].push_back(std::move(rec));
    }
    if (ncol == 0) throw std::runtime_error(path + ": file is empty");
    // Deterministic record order, so which record a truth event is attributed to cannot depend on the
    // order rows happen to sit in the file.
    for (auto& [bid, recs] : cv.records)
        std::sort(recs.begin(), recs.end(),
                  [](const CalledRecord& a, const CalledRecord& b) { return a.variant_id < b.variant_id; });
    return cv;
}

std::vector<std::string> split_on(const std::string& s, char sep) {
    std::vector<std::string> f;
    std::string cur;
    for (char c : s) { if (c == sep) { f.push_back(cur); cur.clear(); } else cur += c; }
    f.push_back(cur);
    return f;
}

// One VCF record, reduced to what a reconstruction needs. Genotypes are haploid because each VCF
// sample column IS one haplotype path -- the join to the graph is by name, with no phasing to infer.
struct VcfRecord {
    std::string id;
    std::size_t bubble_id = 0;
    bool has_bubble_id = false;
    std::size_t pos = 0;               // 1-based, genomic
    std::size_t end = 0;
    long long svlen = 0;
    std::string svtype;
    std::string insseq;
    std::string ref_allele;                 // REF column, when the record carries explicit sequence
    std::vector<std::string> alts;          // ALT alleles, explicit sequence only (empty if symbolic)
    std::size_t ref_cn = 0;            // copy-number records: copies the reference carries
    std::size_t ru_len = 0;            // copy-number records: length of one repeat unit
    // Wide enough for a bubble with hundreds of distinct alleles: the allele VCF indexes every one of
    // them, so an 8-bit slot silently wraps past 127 and drops the edit.
    std::vector<std::int32_t> gt;      // per sample: -1 missing, 0 reference, >=1 alt
    std::vector<long long> cn;         // per sample: copies this haplotype carries (-1 absent)
    std::vector<long long> cnbp;       // per sample: bp this haplotype gains/loses
};

struct VcfData {
    std::vector<std::string> samples;
    std::vector<VcfRecord> records;
};

VcfData load_vcf(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Failed to open VCF: " + path);
    VcfData vd;
    std::string line;
    bool saw_header = false;
    std::size_t lineno = 0;
    std::unordered_set<std::string> seen_ids;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line.rfind("##", 0) == 0) continue;
        if (line[0] == '#') {
            const std::vector<std::string> h = split_on(line, '\t');
            for (std::size_t i = 9; i < h.size(); ++i) vd.samples.push_back(h[i]);
            // A repeated sample column makes "which haplotype is this" depend on column order, and the
            // join below would silently keep only one of them.
            std::unordered_set<std::string> seen;
            for (const std::string& s : vd.samples)
                if (!seen.insert(s).second)
                    throw std::runtime_error("VCF has a duplicate sample column '" + s +
                                             "', so which haplotype a genotype belongs to is undefined: " + path);
            saw_header = true;
            continue;
        }
        ++lineno;
        const std::vector<std::string> f = split_on(line, '\t');
        if (!saw_header)
            throw std::runtime_error(path + ": a data line appears before the #CHROM header");
        // A short or long line means the columns are not where they are read from, so every genotype
        // after it belongs to the wrong haplotype.
        if (f.size() != 9 + vd.samples.size())
            throw std::runtime_error(path + ":" + std::to_string(lineno) + ": has " +
                                     std::to_string(f.size()) + " fields, the #CHROM header declares " +
                                     std::to_string(9 + vd.samples.size()) + " (9 fixed + " +
                                     std::to_string(vd.samples.size()) + " samples)");
        VcfRecord r;
        r.id = f[2];
        if (r.id != "." && !seen_ids.insert(r.id).second)
            throw std::runtime_error(path + ":" + std::to_string(lineno) + ": duplicate record ID '" +
                                     r.id + "'");
        r.pos = static_cast<std::size_t>(std::stoull(f[1]));
        // A record whose ALT is literal sequence is applied as a plain REF -> ALT substitution, which
        // is what the allele VCF emits: one record per bubble, every distinct allele spelled out.
        if (!f[4].empty() && f[4] != "." && f[4].find('<') == std::string::npos) {
            r.ref_allele = f[3];
            r.alts = split_on(f[4], ',');
        }
        for (const std::string& kv : split_on(f[7], ';')) {
            const std::size_t eq = kv.find('=');
            if (eq == std::string::npos) continue;
            const std::string k = kv.substr(0, eq), v = kv.substr(eq + 1);
            if (k == "BUBBLE_ID") { r.bubble_id = static_cast<std::size_t>(std::stoull(v)); r.has_bubble_id = true; }
            else if (k == "END") r.end = static_cast<std::size_t>(std::stoull(v));
            else if (k == "SVLEN") r.svlen = std::stoll(v);
            else if (k == "SVTYPE") r.svtype = v;
            else if (k == "INSSEQ") r.insseq = v;
            else if (k == "REF_CN") r.ref_cn = static_cast<std::size_t>(std::stoull(v));
            // RU_LEN only appears on CN_METHOD=REP now; a MODULE_BP record reports the same
            // quantity as CN_UNIT_BP, renamed because it is the SHARED per-copy content rather
            // than a literal repeat unit. --dup-model cn wants whichever the record carries.
            else if (k == "RU_LEN" || k == "CN_UNIT_BP")
                r.ru_len = static_cast<std::size_t>(std::stoull(v));
        }
        // Without BUBBLE_ID a record cannot be placed, and every such record would silently pile up
        // against bubble 0.
        if (!r.has_bubble_id)
            throw std::runtime_error(path + ":" + std::to_string(lineno) + ": record '" + r.id +
                                     "' carries no BUBBLE_ID, so it cannot be attached to a site");
        const std::vector<std::string> fmt = split_on(f[8], ':');
        std::size_t gt_i = fmt.size(), cnbp_i = fmt.size(), cn_i = fmt.size();
        for (std::size_t i = 0; i < fmt.size(); ++i) {
            if (fmt[i] == "GT") gt_i = i;
            else if (fmt[i] == "CNBP") cnbp_i = i;
            else if (fmt[i] == "CN") cn_i = i;
        }
        r.gt.assign(vd.samples.size(), -1);
        r.cn.assign(vd.samples.size(), -1);
        r.cnbp.assign(vd.samples.size(), 0);
        for (std::size_t s = 0; s < vd.samples.size(); ++s) {
            const std::vector<std::string> v = split_on(f[s + 9], ':');
            if (gt_i < v.size() && !v[gt_i].empty() && v[gt_i] != ".") {
                // Every sample column IS one haplotype, so a genotype is one allele index. A diploid
                // field would be read by stoi as its first allele -- "0/1" as 0 -- so a heterozygous
                // carrier would silently score as reference.
                const std::string& g = v[gt_i];
                if (g.find('/') != std::string::npos || g.find('|') != std::string::npos)
                    throw std::runtime_error(path + ":" + std::to_string(lineno) + ": sample '" +
                                             vd.samples[s] + "' has the diploid genotype '" + g +
                                             "'. benchmark joins one VCF column to one haplotype path, "
                                             "so genotypes must be haploid allele indices");
                std::size_t used = 0;
                int gi = 0;
                try { gi = std::stoi(g, &used); } catch (const std::exception&) { used = 0; }
                if (used != g.size() || gi < 0)
                    throw std::runtime_error(path + ":" + std::to_string(lineno) + ": sample '" +
                                             vd.samples[s] + "' has the genotype '" + g +
                                             "', which is not an allele index");
                r.gt[s] = static_cast<std::int32_t>(gi);
            }
            if (cnbp_i < v.size() && !v[cnbp_i].empty() && v[cnbp_i] != ".")
                r.cnbp[s] = std::stoll(v[cnbp_i]);
            if (cn_i < v.size() && !v[cn_i].empty() && v[cn_i] != ".")
                r.cn[s] = std::stoll(v[cn_i]);
        }
        vd.records.push_back(std::move(r));
    }
    if (!saw_header) throw std::runtime_error("VCF has no #CHROM header line: " + path);
    return vd;
}

// How many records of each kind could actually be laid down, so a genotype QV is never read without
// knowing what it was computed from.
struct GtStats {
    std::size_t applied = 0;
    std::size_t unplaceable = 0;   // POS outside the bubble's reference span (a merged or shifted bubble)
    std::size_t clamped = 0;       // edit ran past the span end and was truncated
    std::size_t unhandled = 0;     // svtype with no reconstruction rule
    std::size_t ref_mismatch = 0;  // record REF did not match the reference at POS (placement is wrong)
    std::size_t heuristic = 0;     // laid down by tiling an inferred span (DUP): length, not sequence
};

// Apply the edits `hap` is genotyped as carrying to the bubble's reference sequence. Nothing from the
// haplotype's own walk enters here -- that is the whole point of the second metric. Edits are laid down
// right to left so earlier ones do not shift the coordinates of later ones.
std::string apply_genotype(const std::string& ref_bubble,
                           std::size_t bubble_start_0,          // 0-based bp offset of ref_bubble in the reference path
                           std::size_t region_start_1based,
                           const std::vector<const VcfRecord*>& recs,
                           std::size_t hap,
                           bool cn_model,
                           GtStats& stats) {
    std::vector<const VcfRecord*> carried;
    for (const VcfRecord* r : recs)
        if (hap < r->gt.size() && r->gt[hap] >= 1) carried.push_back(r);
    std::sort(carried.begin(), carried.end(),
              [](const VcfRecord* a, const VcfRecord* b) { return a->pos > b->pos; });

    std::string out = ref_bubble;
    for (const VcfRecord* r : carried) {
        // POS is the anchor base; the event occupies POS+1 onwards. Everything is expressed relative to
        // the bubble's own reference sequence so the two metrics share a denominator.
        if (r->pos < region_start_1based) { ++stats.unplaceable; continue; }
        const std::size_t g = r->pos - region_start_1based;   // 0-based in the reference path
        if (g < bubble_start_0 || g >= bubble_start_0 + out.size()) { ++stats.unplaceable; continue; }
        const std::size_t at = g - bubble_start_0 + 1;        // first edited base, 0-based in `out`
        if (at > out.size()) { ++stats.unplaceable; continue; }

        if (!r->alts.empty()) {
            // Explicit-sequence allele: replace the record's REF span with the allele the GT names.
            // `at` is one past the anchor, and REF starts AT the anchor, so back up one base.
            const std::size_t start = at - 1;
            const std::int32_t gi = r->gt[hap];
            if (gi < 1 || static_cast<std::size_t>(gi) > r->alts.size()) { ++stats.unhandled; continue; }
            std::size_t n = r->ref_allele.size();
            if (start + n > out.size()) { n = out.size() - start; ++stats.clamped; }
            // The record's REF must be the reference sequence at POS. If it is not, the record is
            // being laid down in the wrong place and every score downstream is fiction.
            if (out.compare(start, n, r->ref_allele, 0, n) != 0) ++stats.ref_mismatch;
            out.erase(start, n);
            out.insert(start, r->alts[static_cast<std::size_t>(gi) - 1]);
            ++stats.applied;
        } else if (r->svtype == "DEL") {
            std::size_t n = static_cast<std::size_t>(r->svlen < 0 ? -r->svlen : r->svlen);
            if (at + n > out.size()) { n = out.size() - at; ++stats.clamped; }
            out.erase(at, n);
            ++stats.applied;
        } else if (r->svtype == "INS") {
            if (r->insseq.empty()) { ++stats.unhandled; continue; }
            out.insert(at, r->insseq);
            ++stats.applied;
        } else if (r->svtype == "INV") {
            // Symbolic inversion: the span POS+1..END comes back reverse-complemented. Unlike DUP this
            // is exact -- an inversion carries no new sequence, so the reference supplies all of it.
            std::size_t n = r->end > r->pos ? r->end - r->pos
                                            : static_cast<std::size_t>(r->svlen < 0 ? -r->svlen : r->svlen);
            if (n == 0) { ++stats.unhandled; continue; }
            if (at + n > out.size()) { n = out.size() - at; ++stats.clamped; }
            if (n == 0) { ++stats.unhandled; continue; }
            const std::string seg = out.substr(at, n);
            out.replace(at, n, reverse_complement(seg));
            ++stats.applied;
        } else if (r->svtype == "DUP") {
            // A copy-number record. The reference carries REF_CN copies of a RU_LEN unit here and this
            // haplotype carries CN of them, so lay down CN copies in place of REF_CN -- the array's
            // LENGTH is what this record is about, and it is what drives the score. Both branches tile
            // an inferred reference span, so what is reconstructed is the right LENGTH of approximately
            // right sequence: counted as `heuristic`, never as exact.
            const std::size_t avail = out.size() - at;
            const long long cn = hap < r->cn.size() ? r->cn[hap] : -1;
            if (cn_model && cn >= 0 && r->ru_len > 0 && r->ref_cn > 0) {
                const std::size_t unit_len = std::min(r->ru_len, avail);
                if (unit_len == 0) { ++stats.unhandled; continue; }
                const std::string unit = out.substr(at, unit_len);
                std::size_t drop = r->ref_cn * r->ru_len;
                if (drop > avail) { drop = avail; ++stats.clamped; }
                out.erase(at, drop);
                const std::size_t want = static_cast<std::size_t>(cn) * r->ru_len;
                std::string add;
                add.reserve(want);
                while (add.size() < want) add += unit.substr(0, std::min(unit.size(), want - add.size()));
                out.insert(at, add);
                ++stats.applied;
                ++stats.heuristic;
                continue;
            }
            // No CN or no unit length: fall back to the record's own bp delta over the annotated span.
            const long long d = hap < r->cnbp.size() ? r->cnbp[hap] : 0;
            if (d == 0) { ++stats.unhandled; continue; }
            const std::size_t unit = std::min(r->end > r->pos ? r->end - r->pos : 0, avail);
            if (unit == 0) { ++stats.unhandled; continue; }
            if (d < 0) {
                std::size_t n = static_cast<std::size_t>(-d);
                if (n > avail) { n = avail; ++stats.clamped; }
                out.erase(at, n);
            } else {
                std::string add;
                add.reserve(static_cast<std::size_t>(d));
                while (add.size() < static_cast<std::size_t>(d))
                    add += out.substr(at, std::min(unit, static_cast<std::size_t>(d) - add.size()));
                out.insert(at, add);
            }
            ++stats.applied;
            ++stats.heuristic;
        } else {
            ++stats.unhandled;
        }
    }
    return out;
}

// ---- One decomposition of the two walks, shared by the ledger and by both reconstructions ----
//
// LCS anchors on step tokens: the maximal runs of steps between consecutive shared anchors are the
// places the two walks disagree. Coordinate-free, so bubble orientation is irrelevant.
std::vector<std::pair<std::size_t, std::size_t>> walk_anchors(const std::vector<PathStep>& ref_steps,
                                                              const std::vector<PathStep>& hap_steps,
                                                              bool& coarse) {
    coarse = false;
    const std::vector<std::uint64_t> rt = build_walk_tokens(ref_steps);
    const std::vector<std::uint64_t> ht = build_walk_tokens(hap_steps);
    const std::size_t n = rt.size(), m = ht.size();
    std::vector<std::pair<std::size_t, std::size_t>> anchors;
    // Guard the O(nm) table on pathological walks (a huge DUP). Falling back means the whole bubble
    // becomes ONE block, which is a coarser decomposition, not a wrong one -- but it has to be counted,
    // because a single block cannot separate a called event from an uncalled one beside it.
    if (n == 0 || m == 0 || static_cast<std::uint64_t>(n) * m > 8000000ULL) {
        if (n != 0 && m != 0) coarse = true;
        return anchors;
    }
    std::vector<std::vector<std::uint32_t>> dp(n + 1, std::vector<std::uint32_t>(m + 1, 0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < m; ++j)
            dp[i + 1][j + 1] = (rt[i] == ht[j]) ? dp[i][j] + 1 : std::max(dp[i][j + 1], dp[i + 1][j]);
    std::size_t i = n, j = m;
    while (i > 0 && j > 0) {
        if (rt[i - 1] == ht[j - 1]) { anchors.emplace_back(i - 1, j - 1); --i; --j; }
        else if (dp[i - 1][j] >= dp[i][j - 1]) --i;
        else --j;
    }
    std::reverse(anchors.begin(), anchors.end());
    return anchors;
}

// One maximal divergent block: the reference takes ref_steps[ri0,ri1) where the haplotype takes
// hap_steps[hi0,hi1). This is the unit `call` works in -- a difference between two walks -- so it is
// the unit a truth ledger has to be built on. Its SIZE is a property of the two walks and of nothing
// else: it does not move with the alignment, the threshold or the calls.
struct TruthBlock {
    std::size_t ri0 = 0, ri1 = 0, hi0 = 0, hi1 = 0;
    std::size_t ref_bp = 0, hap_bp = 0;
    std::size_t size_bp = 0;   // max(ref_bp, hap_bp): call gates a linked DEL/INS pair on the larger arm
    int rec = -1;              // index into the bubble's records; -1 when no emitted call covers it
};

std::vector<TruthBlock> divergent_blocks(const Graph& graph,
                                         const std::vector<PathStep>& ref_steps,
                                         const std::vector<PathStep>& hap_steps,
                                         const std::vector<std::pair<std::size_t, std::size_t>>& anchors,
                                         const std::vector<CalledRecord>& recs) {
    auto node_bp = [&](const std::vector<PathStep>& steps, std::size_t lo, std::size_t hi) {
        std::size_t bp = 0;
        for (std::size_t k = lo; k < hi; ++k) {
            const auto it = graph.nodes.find(steps[k].node_id);
            if (it != graph.nodes.end()) bp += it->second.sequence.size();
        }
        return bp;
    };
    std::vector<TruthBlock> out;
    auto add = [&](std::size_t ri0, std::size_t ri1, std::size_t hi0, std::size_t hi1) {
        if (ri1 == ri0 && hi1 == hi0) return;
        TruthBlock b;
        b.ri0 = ri0; b.ri1 = ri1; b.hi0 = hi0; b.hi1 = hi1;
        b.ref_bp = node_bp(ref_steps, ri0, ri1);
        b.hap_bp = node_bp(hap_steps, hi0, hi1);
        b.size_bp = std::max(b.ref_bp, b.hap_bp);
        // Attribute to the FIRST record (by variant id) whose nodes this block actually contains. A
        // block that merely sits in the same bubble as a record is not attributed to it.
        for (std::size_t ri = 0; ri < recs.size() && b.rec < 0; ++ri) {
            for (std::size_t k = ri0; k < ri1 && b.rec < 0; ++k)
                if (recs[ri].nodes.count(ref_steps[k].node_id)) b.rec = static_cast<int>(ri);
            for (std::size_t k = hi0; k < hi1 && b.rec < 0; ++k)
                if (recs[ri].nodes.count(hap_steps[k].node_id)) b.rec = static_cast<int>(ri);
        }
        out.push_back(b);
    };
    std::size_t pi = 0, pj = 0;
    for (const auto& a : anchors) {
        add(pi, a.first, pj, a.second);
        pi = a.first + 1;
        pj = a.second + 1;
    }
    add(pi, ref_steps.size(), pj, hap_steps.size());
    return out;
}

// Walk the reference, substituting the haplotype's own steps at every block `take` accepts. `applied`
// reports whether anything was substituted.
std::vector<PathStep> reconstruct(const std::vector<PathStep>& ref_steps,
                                  const std::vector<PathStep>& hap_steps,
                                  const std::vector<std::pair<std::size_t, std::size_t>>& anchors,
                                  const std::vector<TruthBlock>& blocks,
                                  const std::function<bool(const TruthBlock&)>& take,
                                  bool& applied) {
    applied = false;
    std::unordered_map<std::size_t, const TruthBlock*> by_ri0;
    for (const TruthBlock& b : blocks) by_ri0[b.ri0 * 1000003ULL + b.hi0] = &b;
    std::vector<PathStep> out;
    std::size_t pi = 0, pj = 0;
    auto emit_block = [&](std::size_t ri1, std::size_t hi1) {
        const auto it = by_ri0.find(pi * 1000003ULL + pj);
        const bool substitute = it != by_ri0.end() && take(*it->second);
        if (substitute) {
            for (std::size_t k = pj; k < hi1; ++k) out.push_back(hap_steps[k]);
            if (ri1 != pi || hi1 != pj) applied = true;
        } else {
            for (std::size_t k = pi; k < ri1; ++k) out.push_back(ref_steps[k]);
        }
    };
    for (const auto& a : anchors) {
        emit_block(a.first, a.second);
        out.push_back(ref_steps[a.first]);  // shared anchor step
        pi = a.first + 1;
        pj = a.second + 1;
    }
    emit_block(ref_steps.size(), hap_steps.size());
    return out;
}

// The truth ledger for one (haplotype, bubble): what was there to be found, and what became of it.
struct Ledger {
    std::size_t events = 0;
    std::size_t called = 0, missed = 0, below = 0;
    std::size_t called_bp = 0, missed_bp = 0, below_bp = 0;
};

struct BubbleObs {
    std::size_t bubble_id = 0;
    bool carrier = false;
    std::size_t delta = 0;
    std::size_t aln_len = 0;
    std::size_t run_lt_bp = 0;   // residual in alignment runs < min-sv-bp
    std::size_t run_ge_bp = 0;   // residual in alignment runs >= min-sv-bp
    Ledger led;
    // The called-only reconstruction: the reference with exactly the retained calls implanted.
    std::size_t called_delta = 0;
    std::size_t called_aln_len = 0;
    // Genotype-aware reconstruction over the same bubble against the same truth, plus the do-nothing
    // baseline (plain reference) it has to beat. A VCF that closes none of the gap between the baseline
    // and the graph bound is worth nothing here, and one that scores WORSE than the baseline is applying
    // edits in the wrong place -- which is the invariant that makes this number readable.
    bool gt_scored = false;
    std::size_t gt_delta = 0;
    std::size_t gt_aln_len = 0;
    std::size_t ref_delta = 0;
    std::size_t ref_aln_len = 0;
    bool gt_called_carrier = false;   // the VCF genotypes this haplotype as carrying an alt here
    bool gt_true_carrier = false;     // its walk really does diverge, by at least min-sv-bp
};

struct EventRow {
    std::string sample;
    std::size_t bubble_id = 0;
    std::size_t ref_bp = 0, hap_bp = 0, size_bp = 0;
    const char* klass = "";
    std::string variant_id;
    std::string svtype;
};

struct HapResult {
    std::string sample;
    bool scored = false;
    bool joined = false;
    std::size_t sum_delta = 0;
    std::size_t sum_aln = 0;
    std::size_t sum_run_lt = 0;
    std::size_t sum_run_ge = 0;
    std::size_t called_sum_delta = 0;
    std::size_t called_sum_aln = 0;
    Ledger led;
    bool gt_scored = false;
    std::size_t gt_sum_delta = 0;
    std::size_t gt_sum_aln = 0;
    std::size_t ref_sum_delta = 0;
    std::size_t ref_sum_aln = 0;
    std::size_t coarse_blocks = 0;     // bubbles whose walks were too big to decompose
    std::size_t skipped_bubbles = 0;   // bubbles this haplotype does not traverse
    GtStats gt_stats;
    std::vector<BubbleObs> obs;
    std::vector<EventRow> events;
};

} // namespace

int run_benchmark_command(const std::vector<std::string>& args) {
    if (args.empty()) { print_benchmark_help(); return 0; }

    std::string gfa_path, bubble_prefix_in, bubbles_csv_in, variant_nodes_in, vcf_in, reference_path, out_prefix;
    std::size_t threads = 0;
    std::size_t min_sv_bp = 50;   // event-size threshold for the truth ledger (match the call run)
    bool quiet = false;
    // Every reference-traversed bubble is scored. Restricting to called bubbles cannot measure a miss:
    // the bubble a caller said nothing about is precisely the one a miss lives in.
    bool called_only = false;
    bool write_events = true;
    // How a copy-number record is laid down. `cnbp` -- the per-sample bp delta -- is the default
    // because it is measured: on acot it closes 95.7% of the gap against 51.5% for `cn`, since
    // CN x RU_LEN understates the true bp change by roughly half at a cyclic array.
    bool dup_cn_model = false;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) throw std::runtime_error("Missing value after " + flag);
            return args[++i];
        };
        if (arg == "-h" || arg == "--help") { print_benchmark_help(); return 0; }
        else if (arg == "-i" || arg == "--gfa") gfa_path = require_value(arg);
        else if (arg == "-b" || arg == "--bubble-prefix-in") bubble_prefix_in = require_value(arg);
        else if (arg == "-c" || arg == "--bubbles-csv-in") bubbles_csv_in = require_value(arg);
        else if (arg == "--variant-nodes") variant_nodes_in = require_value(arg);
        else if (arg == "--vcf") vcf_in = require_value(arg);
        else if (arg == "--dup-model") {
            const std::string v = require_value(arg);
            if (v == "cn") dup_cn_model = true;
            else if (v == "cnbp") dup_cn_model = false;
            else throw std::runtime_error("--dup-model must be cn or cnbp");
        }
        else if (arg == "-r" || arg == "--reference-path") reference_path = require_value(arg);
        else if (arg == "-o" || arg == "--out-prefix") out_prefix = require_value(arg);
        else if (arg == "--called-bubbles-only") called_only = true;
        // Retained: scoring every bubble is now the default, so this asks for what it already gets.
        else if (arg == "--all-bubbles") called_only = false;
        else if (arg == "--no-truth-events") write_events = false;
        else if (arg == "--min-sv-bp") min_sv_bp = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--threads") threads = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "-q" || arg == "--quiet") quiet = true;
        else throw std::runtime_error("Unknown option for benchmark: " + arg);
    }

    if (gfa_path.empty()) throw std::runtime_error("Missing required input: --gfa <path>");
    if (!bubble_prefix_in.empty()) {
        const std::string derived = bubble_prefix_in + ".bubbles.csv";
        if (bubbles_csv_in.empty()) bubbles_csv_in = derived;
        else if (bubbles_csv_in != derived)
            throw std::runtime_error("Conflicting bubble inputs: --bubble-prefix-in resolves to '" +
                                     derived + "' but --bubbles-csv-in is '" + bubbles_csv_in + "'");
    }
    if (bubbles_csv_in.empty())
        throw std::runtime_error("benchmark requires --bubble-prefix-in <prefix> or --bubbles-csv-in <path>");
    if (variant_nodes_in.empty()) throw std::runtime_error("benchmark requires --variant-nodes <path>");
    if (reference_path.empty()) throw std::runtime_error("--reference-path is required for module 'benchmark'");
    if (out_prefix.empty()) throw std::runtime_error("benchmark requires -o/--out-prefix <prefix>");
    if (min_sv_bp == 0) throw std::runtime_error("--min-sv-bp must be at least 1");

    // Every destination is resolved and checked against every input BEFORE anything is opened, then
    // written to a staged sibling and renamed in only once the run has succeeded.
    std::vector<std::string> finals = {out_prefix + ".qv.tsv",
                                       out_prefix + ".qv_by_haplotype.tsv",
                                       out_prefix + ".qv_summary.tsv"};
    if (write_events) finals.push_back(out_prefix + ".truth_events.tsv");
    {
        const std::vector<std::string> inputs = {gfa_path, bubbles_csv_in, variant_nodes_in, vcf_in};
        std::unordered_set<std::string> seen_out;
        for (const std::string& f : finals) {
            if (!seen_out.insert(f).second)
                throw std::runtime_error("benchmark: two outputs would be written to the same file: " + f);
            for (const std::string& in : inputs) {
                if (in.empty()) continue;
                std::error_code ec;
                if (f == in || std::filesystem::equivalent(f, in, ec))
                    throw std::runtime_error("benchmark: output '" + f + "' is also an input; refusing to "
                                             "overwrite what the run reads");
            }
        }
    }

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(gfa_path, parse_options);
    if (graph.paths.empty())
        throw std::runtime_error("Input GFA has no P/W paths; benchmarking requires paths");
    // benchmark spells every walk by concatenating whole segments and measures bp from them, so a step
    // naming a missing node, a step pair with no link, a duplicate path name and a non-zero overlap all
    // corrupt the numbers rather than failing.
    validate_graph_paths(graph, "benchmark", true, true);
    const std::string ref_name = resolve_reference_path_name(graph, reference_path, "benchmark");

    cli::RunLog log("benchmark", quiet);
    log.info("input " + gfa_path + " (" + std::to_string(graph.nodes.size()) + " nodes, " +
             std::to_string(graph.paths.size()) + " paths); reference " + ref_name);

    const std::vector<Bubble> all_bubbles_csv = read_bubbles_csv(bubbles_csv_in);
    std::unordered_map<std::size_t, std::unordered_set<std::string>> bubble_members;
    {
        for (const Bubble& b : all_bubbles_csv) {
            auto need = [&](const std::string& n) {
                if (!graph.nodes.count(n))
                    throw std::runtime_error("bubble " + std::to_string(b.id) + " names node '" + n +
                                             "', which is not in " + gfa_path +
                                             " -- the CSV and the graph do not describe the same data");
                return n;
            };
            std::unordered_set<std::string> members;
            members.insert(need(b.source));
            members.insert(need(b.sink));
            for (const std::string& n : b.inside) members.insert(need(n));
            if (!bubble_members.emplace(b.id, std::move(members)).second)
                throw std::runtime_error("bubbles CSV has a duplicate bubble id " + std::to_string(b.id) +
                                         ": which site a call or a score refers to would be undefined");
        }
    }
    const CalledVariants cv = load_variant_nodes(variant_nodes_in, bubble_members);
    const PathRecord* ref_path = nullptr;
    for (const PathRecord& p : graph.paths) if (p.name == ref_name) ref_path = &p;

    std::vector<Bubble> bubbles;
    for (const Bubble& b : all_bubbles_csv) if (!called_only || cv.bubble_ids.count(b.id)) bubbles.push_back(b);

    // Reference walk per scored bubble (the reconstruction backbone). A bubble the reference doesn't
    // traverse has no baseline, so it is dropped from scoring -- counted, never silent.
    const BubblePathIndex ref_index = build_bubble_path_index(*ref_path);
    std::unordered_map<std::size_t, std::vector<PathStep>> ref_steps;
    for (const Bubble& b : bubbles) {
        const auto s = bubble_steps(*ref_path, ref_index, b);
        if (s) ref_steps[b.id] = *s;
    }
    const std::size_t no_ref_bubbles = bubbles.size() - ref_steps.size();
    log.info("scoring " + std::to_string(ref_steps.size()) + " of " + std::to_string(bubbles.size()) +
             (called_only ? " called bubbles" : " bubbles") + " across " +
             std::to_string(graph.paths.size()) + " haplotypes (" + std::to_string(no_ref_bubbles) +
             " not traversed by the reference), event threshold " + std::to_string(min_sv_bp) + " bp");

    // Genotype-aware scoring, when a VCF is supplied. VCF POS lives in the reference path's own forward
    // coordinates, so each bubble's reference span is taken as a substring there and the reconstruction
    // is flipped back for bubbles the reference crosses sink->source.
    const bool do_gt = !vcf_in.empty();
    VcfData vcf;
    std::unordered_map<std::size_t, std::vector<const VcfRecord*>> recs_by_bubble;
    std::vector<int> hap_vcf(graph.paths.size(), -1);
    std::unordered_map<std::size_t, std::pair<std::size_t, std::size_t>> ref_span;  // bubble -> [bp start, bp len)
    std::unordered_map<std::size_t, bool> ref_flip;
    std::string ref_seq;
    std::size_t region_start = 1;
    std::size_t joined = 0, unjoined_samples = 0;
    if (do_gt) {
        vcf = load_vcf(vcf_in);
        for (const VcfRecord& r : vcf.records) recs_by_bubble[r.bubble_id].push_back(&r);
        std::unordered_map<std::string, std::size_t> sample_of;
        for (std::size_t i = 0; i < vcf.samples.size(); ++i) sample_of[vcf.samples[i]] = i;
        std::unordered_set<std::size_t> used;
        for (std::size_t i = 0; i < graph.paths.size(); ++i) {
            const auto it = sample_of.find(graph.paths[i].name);
            if (it != sample_of.end()) { hap_vcf[i] = static_cast<int>(it->second); used.insert(it->second); ++joined; }
        }
        unjoined_samples = vcf.samples.size() - used.size();
        if (joined == 0)
            throw std::runtime_error("No VCF sample column matches a graph path name, so no haplotype "
                                     "can be genotype-scored: " + vcf_in);
        region_start = parse_reference_path_label(ref_path->name).region_start_1based;

        std::vector<std::size_t> prefix(ref_path->steps.size() + 1, 0);
        for (std::size_t i = 0; i < ref_path->steps.size(); ++i) {
            const auto nit = graph.nodes.find(ref_path->steps[i].node_id);
            prefix[i + 1] = prefix[i] + (nit != graph.nodes.end() ? nit->second.sequence.size() : 0);
        }
        ref_seq = spell_path_steps_sequence(graph, ref_path->steps);
        for (const Bubble& b : bubbles) {
            if (!ref_steps.count(b.id)) continue;
            std::size_t lo = 0, hi = 0;
            bool flip = false, ok = false;
            const auto iv = find_best_bubble_path_interval(ref_index, b);
            if (iv.has_value()) {
                lo = std::min(iv->left, iv->right);
                hi = std::max(iv->left, iv->right);
                flip = !iv->source_to_sink;
                ok = true;
            } else {
                // The same fallback bubble_steps uses: an adjacent source/sink pair with no interior.
                const auto si = ref_index.positions.find(b.source);
                const auto ki = ref_index.positions.find(b.sink);
                if (si != ref_index.positions.end() && ki != ref_index.positions.end()) {
                    const std::unordered_set<std::size_t> sink_pos(ki->second.begin(), ki->second.end());
                    for (const std::size_t p : si->second)
                        if (sink_pos.count(p + 1)) { lo = p; hi = p + 1; ok = true; break; }
                    if (!ok) {
                        const std::unordered_set<std::size_t> src_pos(si->second.begin(), si->second.end());
                        for (const std::size_t p : ki->second)
                            if (src_pos.count(p + 1)) { lo = p; hi = p + 1; flip = true; ok = true; break; }
                    }
                }
            }
            if (!ok || hi + 1 >= prefix.size()) continue;
            ref_span[b.id] = { prefix[lo], prefix[hi + 1] - prefix[lo] };
            ref_flip[b.id] = flip;
        }
        log.info("genotype scoring from " + vcf_in + ": " + std::to_string(vcf.records.size()) +
                 " records, " + std::to_string(joined) + "/" + std::to_string(graph.paths.size()) +
                 " haplotypes joined to VCF samples, " + std::to_string(ref_span.size()) +
                 "/" + std::to_string(bubbles.size()) + " bubbles placed on the reference");
        // A partial join is a selected subset, and a QV computed over a subset is not the QV of the run.
        if (joined != graph.paths.size() || unjoined_samples)
            log.info("warning: the genotype score covers " + std::to_string(joined) + " of " +
                     std::to_string(graph.paths.size()) + " graph haplotypes; " +
                     std::to_string(graph.paths.size() - joined) + " path(s) have no VCF column and " +
                     std::to_string(unjoined_samples) + " VCF sample(s) have no path. It is a SUBSET score");
    }

    std::vector<HapResult> results(graph.paths.size());
    static const std::vector<CalledRecord> kNoCalled;
    run_parallel(graph.paths.size(), threads, [&](std::size_t pi) {
        const PathRecord& path = graph.paths[pi];
        HapResult& hr = results[pi];
        hr.sample = path.name;
        hr.joined = do_gt && hap_vcf[pi] >= 0;
        if (&path == ref_path) return;
        const BubblePathIndex idx = build_bubble_path_index(path);
        for (const Bubble& b : bubbles) {
            const auto rit = ref_steps.find(b.id);
            if (rit == ref_steps.end()) continue;                // reference doesn't traverse -> no baseline
            const auto steps = bubble_steps(path, idx, b);
            if (!steps) { ++hr.skipped_bubbles; continue; }      // haplotype doesn't traverse
            const std::string truth = spell_path_steps_sequence(graph, *steps);
            const auto crit = cv.records.find(b.id);
            const std::vector<CalledRecord>& recs = crit != cv.records.end() ? crit->second : kNoCalled;

            // ONE decomposition, three consumers: the ledger, the optimistic graph reconstruction and
            // the strict called-only one. They cannot disagree about what a block is.
            bool coarse = false;
            const auto anchors = walk_anchors(rit->second, *steps, coarse);
            if (coarse) ++hr.coarse_blocks;
            const std::vector<TruthBlock> blocks =
                divergent_blocks(graph, rit->second, *steps, anchors, recs);

            BubbleObs o;
            o.bubble_id = b.id;
            for (const TruthBlock& blk : blocks) {
                if (blk.size_bp == 0) continue;                  // identical content, nothing to find
                ++o.led.events;
                const char* klass;
                if (blk.size_bp < min_sv_bp) { ++o.led.below; o.led.below_bp += blk.size_bp; klass = "below_threshold"; }
                // `called` = a specific record shares at least one node with this block. That is
                // attribution, not coverage: it does not establish that the record spans the block or
                // represents it correctly.
                else if (blk.rec >= 0)       { ++o.led.called; o.led.called_bp += blk.size_bp; klass = "called"; }
                else                          { ++o.led.missed; o.led.missed_bp += blk.size_bp; klass = "missed"; }
                if (write_events) {
                    EventRow er;
                    er.sample = path.name;
                    er.bubble_id = b.id;
                    er.ref_bp = blk.ref_bp;
                    er.hap_bp = blk.hap_bp;
                    er.size_bp = blk.size_bp;
                    er.klass = klass;
                    er.variant_id = blk.rec >= 0 ? recs[static_cast<std::size_t>(blk.rec)].variant_id : ".";
                    er.svtype = blk.rec >= 0 ? recs[static_cast<std::size_t>(blk.rec)].svtype : ".";
                    hr.events.push_back(std::move(er));
                }
            }

            // Optimistic: substitute wherever the block shares a node with ANY call at this bubble.
            const auto nit = cv.nodes.find(b.id);
            static const std::unordered_set<std::string> kEmptyNodes;
            const std::unordered_set<std::string>& union_nodes =
                nit != cv.nodes.end() ? nit->second : kEmptyNodes;
            bool applied = false;
            const std::vector<PathStep> recon_steps =
                reconstruct(rit->second, *steps, anchors, blocks,
                            [&](const TruthBlock& blk) {
                                for (std::size_t k = blk.ri0; k < blk.ri1; ++k)
                                    if (union_nodes.count(rit->second[k].node_id)) return true;
                                for (std::size_t k = blk.hi0; k < blk.hi1; ++k)
                                    if (union_nodes.count((*steps)[k].node_id)) return true;
                                return false;
                            }, applied);
            const std::string recon = spell_path_steps_sequence(graph, recon_steps);
            const NwAlign nw = nw_edit_distance(recon, truth, min_sv_bp);

            // Strict: substitute only blocks that map to a specific emitted call AND reach the threshold.
            bool c_applied = false;
            const std::vector<PathStep> called_steps =
                reconstruct(rit->second, *steps, anchors, blocks,
                            [&](const TruthBlock& blk) { return blk.rec >= 0 && blk.size_bp >= min_sv_bp; },
                            c_applied);
            const NwAlign cnw = nw_edit_distance(spell_path_steps_sequence(graph, called_steps), truth);

            hr.scored = true;
            hr.sum_delta += nw.edits;
            hr.sum_aln += nw.aln_len;
            hr.sum_run_lt += nw.sub_threshold_bp;
            hr.sum_run_ge += nw.over_threshold_bp;
            hr.called_sum_delta += cnw.edits;
            hr.called_sum_aln += cnw.aln_len;
            hr.led.events += o.led.events;
            hr.led.called += o.led.called;   hr.led.called_bp += o.led.called_bp;
            hr.led.missed += o.led.missed;   hr.led.missed_bp += o.led.missed_bp;
            hr.led.below  += o.led.below;    hr.led.below_bp  += o.led.below_bp;

            o.carrier = applied;
            o.delta = nw.edits;
            o.aln_len = nw.aln_len;
            o.run_lt_bp = nw.sub_threshold_bp;
            o.run_ge_bp = nw.over_threshold_bp;
            o.called_delta = cnw.edits;
            o.called_aln_len = cnw.aln_len;
            // Truth for the carrier call comes from the WALKS, not from an alignment: this haplotype
            // carries something callable here iff one of its divergent blocks reaches the threshold.
            o.gt_true_carrier = o.led.called > 0 || o.led.missed > 0;

            const auto span = do_gt && hap_vcf[pi] >= 0 ? ref_span.find(b.id) : ref_span.end();
            if (span != ref_span.end()) {
                const std::size_t hv = static_cast<std::size_t>(hap_vcf[pi]);
                static const std::vector<const VcfRecord*> kNoRecs;
                const auto brit = recs_by_bubble.find(b.id);
                const std::vector<const VcfRecord*>& brecs =
                    brit != recs_by_bubble.end() ? brit->second : kNoRecs;
                const std::string ref_bubble = ref_seq.substr(span->second.first, span->second.second);
                std::string gt_recon = apply_genotype(ref_bubble, span->second.first, region_start,
                                                      brecs, hv, dup_cn_model, hr.gt_stats);
                std::string ref_recon = ref_bubble;
                if (ref_flip.at(b.id)) {
                    gt_recon = reverse_complement(gt_recon);
                    ref_recon = reverse_complement(ref_recon);
                }
                const NwAlign gnw = nw_edit_distance(gt_recon, truth);
                const NwAlign rnw = nw_edit_distance(ref_recon, truth);
                hr.gt_scored = true;
                hr.gt_sum_delta += gnw.edits;
                hr.gt_sum_aln += gnw.aln_len;
                hr.ref_sum_delta += rnw.edits;
                hr.ref_sum_aln += rnw.aln_len;
                o.gt_scored = true;
                o.gt_delta = gnw.edits;
                o.gt_aln_len = gnw.aln_len;
                o.ref_delta = rnw.edits;
                o.ref_aln_len = rnw.aln_len;
                for (const VcfRecord* r : brecs)
                    if (hv < r->gt.size() && r->gt[hv] >= 1) { o.gt_called_carrier = true; break; }
            }
            hr.obs.push_back(std::move(o));
        }
    });

    cli::StagedOutputs staged("benchmark");
    std::vector<std::string> staged_paths;
    for (const std::string& f : finals) staged_paths.push_back(staged.stage(f));
    auto staged_for = [&](const std::string& f) {
        for (std::size_t i = 0; i < finals.size(); ++i) if (finals[i] == f) return staged_paths[i];
        throw std::runtime_error("benchmark: unstaged output " + f);
    };

    std::map<std::string, std::size_t> band_count;
    std::map<std::string, std::size_t> quintile_count;
    std::map<std::string, std::size_t> gt_band_count;
    std::map<std::string, std::size_t> gt_quintile_count;
    std::map<std::string, std::size_t> called_quintile_count;
    std::size_t scored_haps = 0, gt_scored_haps = 0;
    double gt_closed_sum = 0.0;
    std::size_t gt_closed_n = 0, gt_closed_na = 0;
    std::size_t gt_worse_than_ref = 0;
    std::size_t gt_tot_ref_delta = 0, gt_tot_gt_delta = 0, gt_tot_graph_delta = 0;
    std::size_t tot_run_lt = 0, tot_run_ge = 0;
    std::size_t tot_called_delta = 0;
    std::size_t tot_coarse = 0, tot_skipped = 0;
    Ledger tot_led;
    GtStats tot_gt;
    {
        std::ofstream by_hap(staged_for(out_prefix + ".qv_by_haplotype.tsv"));
        if (!by_hap) throw std::runtime_error("Failed to write " + out_prefix + ".qv_by_haplotype.tsv");
        by_hap << "sample\tn_bubbles\tsum_delta\tsum_aln_len\tqv\tband\tqv_max\tqv_ratio\tquintile\tidentity"
                  "\tresid_run_lt_bp\tresid_run_ge_bp"
                  "\ttruth_events\ttruth_called\ttruth_missed\ttruth_below"
                  "\ttruth_called_bp\ttruth_missed_bp\ttruth_below_bp"
                  "\tcalled_sum_delta\tcalled_qv\tcalled_quintile"
                  "\tgt_sum_delta\tgt_sum_aln_len\tgt_qv\tgt_band\tgt_quintile\tgt_identity"
                  "\tref_sum_delta\tref_qv\tgap_closed\n";
        for (const HapResult& hr : results) {
            if (!hr.scored) continue;
            const double qv = qv_from(hr.sum_delta, hr.sum_aln);
            const double qmax = qv_max_from(hr.sum_aln);
            const double ratio = std::min(1.0, std::max(0.0, qmax > 0.0 ? qv / qmax : 1.0));
            // Linear reconstruction identity: 1 - error rate. Unlike QV it does not overweight a few
            // uncalled SNPs over a long region (a handful of edits in tens of kb stays ~1.0).
            const double identity = hr.sum_aln ? std::max(0.0, 1.0 - static_cast<double>(hr.sum_delta) / static_cast<double>(hr.sum_aln)) : 1.0;
            const std::string band = qv_band(qv);
            const std::string quintile = qv_ratio_quintile(ratio);
            const double cqv = qv_from(hr.called_sum_delta, hr.called_sum_aln);
            const double cmax = qv_max_from(hr.called_sum_aln);
            const std::string cquint = qv_ratio_quintile(std::min(1.0, std::max(0.0, cmax > 0.0 ? cqv / cmax : 1.0)));
            by_hap << hr.sample << '\t' << hr.obs.size() << '\t' << hr.sum_delta << '\t'
                   << hr.sum_aln << '\t' << qv << '\t' << band << '\t' << qmax << '\t'
                   << ratio << '\t' << quintile << '\t' << identity << '\t'
                   << hr.sum_run_lt << '\t' << hr.sum_run_ge << '\t'
                   << hr.led.events << '\t' << hr.led.called << '\t' << hr.led.missed << '\t'
                   << hr.led.below << '\t' << hr.led.called_bp << '\t' << hr.led.missed_bp << '\t'
                   << hr.led.below_bp << '\t'
                   << hr.called_sum_delta << '\t' << cqv << '\t' << cquint;
            if (hr.gt_scored) {
                const double gqv = qv_from(hr.gt_sum_delta, hr.gt_sum_aln);
                const double gmax = qv_max_from(hr.gt_sum_aln);
                const double gratio = std::min(1.0, std::max(0.0, gmax > 0.0 ? gqv / gmax : 1.0));
                const double gident = hr.gt_sum_aln ? std::max(0.0, 1.0 - static_cast<double>(hr.gt_sum_delta) / static_cast<double>(hr.gt_sum_aln)) : 1.0;
                const std::string gband = qv_band(gqv);
                const std::string gquint = qv_ratio_quintile(gratio);
                // How much of the distance between doing nothing and the graph bound the VCF closes.
                // 1 = everything the graph could offer, 0 = no better than plain reference, < 0 = the
                // genotypes put edits where they do not belong. When the baseline already equals the
                // graph bound there is no distance to close and the ratio is UNDEFINED -- reporting it
                // as 1.0 called a total miss a perfect score, since baseline == genotype == graph then
                // reads "closed everything".
                const double denom = static_cast<double>(hr.ref_sum_delta) - static_cast<double>(hr.sum_delta);
                const bool has_gap = std::fabs(denom) > 1e-9;
                const double closed = has_gap
                    ? (static_cast<double>(hr.ref_sum_delta) - static_cast<double>(hr.gt_sum_delta)) / denom
                    : 0.0;
                by_hap << '\t' << hr.gt_sum_delta << '\t' << hr.gt_sum_aln << '\t' << gqv << '\t'
                       << gband << '\t' << gquint << '\t' << gident << '\t'
                       << hr.ref_sum_delta << '\t' << qv_from(hr.ref_sum_delta, hr.ref_sum_aln) << '\t';
                if (has_gap) by_hap << closed; else by_hap << "NA";
                if (has_gap) { gt_closed_sum += closed; ++gt_closed_n; } else ++gt_closed_na;
                if (hr.gt_sum_delta > hr.ref_sum_delta) ++gt_worse_than_ref;
                gt_tot_ref_delta += hr.ref_sum_delta;
                gt_tot_gt_delta += hr.gt_sum_delta;
                gt_tot_graph_delta += hr.sum_delta;
                ++gt_band_count[gband];
                ++gt_quintile_count[gquint];
                ++gt_scored_haps;
                tot_gt.applied += hr.gt_stats.applied;
                tot_gt.unplaceable += hr.gt_stats.unplaceable;
                tot_gt.clamped += hr.gt_stats.clamped;
                tot_gt.unhandled += hr.gt_stats.unhandled;
                tot_gt.ref_mismatch += hr.gt_stats.ref_mismatch;
                tot_gt.heuristic += hr.gt_stats.heuristic;
            } else {
                by_hap << "\t.\t.\t.\t.\t.\t.\t.\t.\t.";
            }
            by_hap << '\n';
            ++band_count[band];
            ++quintile_count[quintile];
            ++called_quintile_count[cquint];
            ++scored_haps;
            tot_run_lt += hr.sum_run_lt;
            tot_run_ge += hr.sum_run_ge;
            tot_called_delta += hr.called_sum_delta;
            tot_coarse += hr.coarse_blocks;
            tot_skipped += hr.skipped_bubbles;
            tot_led.events += hr.led.events;
            tot_led.called += hr.led.called;   tot_led.called_bp += hr.led.called_bp;
            tot_led.missed += hr.led.missed;   tot_led.missed_bp += hr.led.missed_bp;
            tot_led.below  += hr.led.below;    tot_led.below_bp  += hr.led.below_bp;
        }
    }

    if (write_events) {
        std::ofstream ev(staged_for(out_prefix + ".truth_events.tsv"));
        if (!ev) throw std::runtime_error("Failed to write " + out_prefix + ".truth_events.tsv");
        ev << "sample\tbubble_id\tref_bp\thap_bp\tsize_bp\tclass\tvariant_id\tsvtype\n";
        for (const HapResult& hr : results)
            for (const EventRow& e : hr.events)
                ev << e.sample << '\t' << e.bubble_id << '\t' << e.ref_bp << '\t' << e.hap_bp << '\t'
                   << e.size_bp << '\t' << e.klass << '\t' << e.variant_id << '\t' << e.svtype << '\n';
    }

    std::map<std::string, std::map<std::string, std::size_t>> class_bands;
    // Carrier confusion at the (haplotype, bubble) level: the genotype call (does the VCF say this
    // haplotype carries an alt here) against the truth (does its walk diverge by at least the
    // threshold). Reported for ALL only -- a bubble judgement has no single svtype, and fanning it out
    // over every svtype the bubble contains reported one observation as several.
    std::array<std::size_t, 4> carrier{{0, 0, 0, 0}};   // {TP, FP, FN, TN}
    // The event ledger by svtype IS a partition: each event maps to at most one record.
    std::map<std::string, std::array<std::size_t, 3>> ledger_by_type;  // svtype -> {called, missed, below}
    {
        std::ofstream qv_out(staged_for(out_prefix + ".qv.tsv"));
        if (!qv_out) throw std::runtime_error("Failed to write " + out_prefix + ".qv.tsv");
        qv_out << "sample\tbubble_id\tsvtypes\tis_carrier\tdelta\taln_len\tqv\tresid_run_lt_bp\tresid_run_ge_bp"
                  "\ttruth_events\ttruth_called\ttruth_missed\ttruth_below\ttruth_missed_bp"
                  "\tcalled_delta\tcalled_qv"
                  "\tgt_delta\tgt_aln_len\tgt_qv\tref_delta\tgt_called_carrier\tgt_true_carrier\n";
        for (const HapResult& hr : results) {
            if (!hr.scored) continue;
            for (const BubbleObs& o : hr.obs) {
                std::string svt;
                const auto sit = cv.svtypes.find(o.bubble_id);
                if (sit != cv.svtypes.end())
                    for (const std::string& s : sit->second) { if (!svt.empty()) svt += ','; svt += s; }
                const double qv = qv_from(o.delta, o.aln_len);
                const std::string band = qv_band(qv);
                qv_out << hr.sample << '\t' << o.bubble_id << '\t' << (svt.empty() ? "." : svt) << '\t'
                       << (o.carrier ? 1 : 0) << '\t' << o.delta << '\t' << o.aln_len << '\t' << qv << '\t'
                       << o.run_lt_bp << '\t' << o.run_ge_bp << '\t'
                       << o.led.events << '\t' << o.led.called << '\t' << o.led.missed << '\t'
                       << o.led.below << '\t' << o.led.missed_bp << '\t'
                       << o.called_delta << '\t' << qv_from(o.called_delta, o.called_aln_len);
                if (o.gt_scored) {
                    qv_out << '\t' << o.gt_delta << '\t' << o.gt_aln_len << '\t'
                           << qv_from(o.gt_delta, o.gt_aln_len) << '\t' << o.ref_delta << '\t'
                           << (o.gt_called_carrier ? 1 : 0) << '\t' << (o.gt_true_carrier ? 1 : 0);
                    ++carrier[o.gt_called_carrier ? (o.gt_true_carrier ? 0 : 1)
                                                  : (o.gt_true_carrier ? 2 : 3)];
                } else {
                    qv_out << "\t.\t.\t.\t.\t.\t" << (o.gt_true_carrier ? 1 : 0);
                }
                qv_out << '\n';
                if (sit != cv.svtypes.end())
                    for (const std::string& s : sit->second) ++class_bands[s][band];
            }
            for (const EventRow& e : hr.events) {
                auto& c = ledger_by_type[e.svtype];
                if (std::string(e.klass) == "called") ++c[0];
                else if (std::string(e.klass) == "missed") ++c[1];
                else ++c[2];
            }
        }
    }

    static const char* kBands[] = {"<17", "17-23", "23-33", ">33"};
    static const char* kQuintiles[] = {"0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"};
    {
        std::ofstream sum(staged_for(out_prefix + ".qv_summary.tsv"));
        if (!sum) throw std::runtime_error("Failed to write " + out_prefix + ".qv_summary.tsv");
        sum << "scope\tkey\tband\tn\tpct\n";
        // Length-fair headline: % of haplotypes per qv/qv_max quintile (Q5 = 0.8-1.0 is best).
        for (const char* q : kQuintiles) {
            const std::size_t n = quintile_count.count(q) ? quintile_count[q] : 0;
            const double pct = scored_haps ? 100.0 * static_cast<double>(n) / static_cast<double>(scored_haps) : 0.0;
            sum << "quintile\tALL\t" << q << '\t' << n << '\t' << pct << '\n';
        }
        for (const char* band : kBands) {
            const std::size_t n = band_count.count(band) ? band_count[band] : 0;
            const double pct = scored_haps ? 100.0 * static_cast<double>(n) / static_cast<double>(scored_haps) : 0.0;
            sum << "haplotype\tALL\t" << band << '\t' << n << '\t' << pct << '\n';
        }
        for (const char* q : kQuintiles) {
            const std::size_t n = called_quintile_count.count(q) ? called_quintile_count[q] : 0;
            const double pct = scored_haps ? 100.0 * static_cast<double>(n) / static_cast<double>(scored_haps) : 0.0;
            sum << "called_quintile\tALL\t" << q << '\t' << n << '\t' << pct << '\n';
        }
        // The truth ledger. Every divergent block between the reference walk and a haplotype's walk is
        // one event, sized from the walks; `called` means it contains a node of a specific emitted
        // record, which is weaker than "that record reproduces it" and is labelled accordingly.
        {
            const std::size_t tot = tot_led.events ? tot_led.events : 1;
            sum << "truth_event\tALL\tevents\t" << tot_led.events << "\t100\n";
            sum << "truth_event\tALL\tcalled\t" << tot_led.called << '\t'
                << 100.0 * static_cast<double>(tot_led.called) / static_cast<double>(tot) << '\n';
            sum << "truth_event\tALL\tmissed\t" << tot_led.missed << '\t'
                << 100.0 * static_cast<double>(tot_led.missed) / static_cast<double>(tot) << '\n';
            sum << "truth_event\tALL\tbelow_threshold\t" << tot_led.below << '\t'
                << 100.0 * static_cast<double>(tot_led.below) / static_cast<double>(tot) << '\n';
            sum << "truth_bp\tALL\tcalled_bp\t" << tot_led.called_bp << "\t0\n";
            sum << "truth_bp\tALL\tmissed_bp\t" << tot_led.missed_bp << "\t0\n";
            sum << "truth_bp\tALL\tbelow_threshold_bp\t" << tot_led.below_bp << "\t0\n";
            for (const auto& [svt, c] : ledger_by_type) {
                sum << "truth_event\t" << svt << "\tcalled\t" << c[0] << "\t0\n";
                sum << "truth_event\t" << svt << "\tmissed\t" << c[1] << "\t0\n";
                sum << "truth_event\t" << svt << "\tbelow_threshold\t" << c[2] << "\t0\n";
            }
        }
        // Residual by contiguous ALIGNMENT RUN length. This is a property of the edit path edlib
        // returns, not of any event: a clean 60 bp deletion of non-repetitive sequence comes back as a
        // dozen short runs, because co-optimal paths distribute the gap across chance matches. Kept as
        // a description of the residual's shape; the called/missed question is the ledger's.
        {
            const std::size_t tot = tot_run_lt + tot_run_ge;
            const double lt_pct = tot ? 100.0 * static_cast<double>(tot_run_lt) / static_cast<double>(tot) : 0.0;
            const double ge_pct = tot ? 100.0 * static_cast<double>(tot_run_ge) / static_cast<double>(tot) : 0.0;
            sum << "residual_run\t<" << min_sv_bp << "bp\tshort_runs\t" << tot_run_lt << '\t' << lt_pct << '\n';
            sum << "residual_run\t>=" << min_sv_bp << "bp\tlong_runs\t" << tot_run_ge << '\t' << ge_pct << '\n';
        }
        // What was left out of the denominators, so no rate is read as if it covered everything.
        sum << "excluded\tALL\tbubbles_no_reference_walk\t" << no_ref_bubbles << "\t0\n";
        sum << "excluded\tALL\thaplotype_bubbles_not_traversed\t" << tot_skipped << "\t0\n";
        sum << "excluded\tALL\tbubbles_decomposed_coarsely\t" << tot_coarse << "\t0\n";
        sum << "excluded\tALL\thaplotypes_without_vcf_column\t"
            << (do_gt ? graph.paths.size() - joined : 0) << "\t0\n";
        sum << "excluded\tALL\tvcf_samples_without_path\t" << unjoined_samples << "\t0\n";
        sum << "called_recon\tALL\tdelta\t" << tot_called_delta << "\t0\n";
        // Genotype-aware reconstruction, reported on the same scale so the two are read side by side.
        if (gt_scored_haps) {
            for (const char* q : kQuintiles) {
                const std::size_t n = gt_quintile_count.count(q) ? gt_quintile_count[q] : 0;
                const double pct = 100.0 * static_cast<double>(n) / static_cast<double>(gt_scored_haps);
                sum << "gt_quintile\tALL\t" << q << '\t' << n << '\t' << pct << '\n';
            }
            for (const char* band : kBands) {
                const std::size_t n = gt_band_count.count(band) ? gt_band_count[band] : 0;
                const double pct = 100.0 * static_cast<double>(n) / static_cast<double>(gt_scored_haps);
                sum << "gt_haplotype\tALL\t" << band << '\t' << n << '\t' << pct << '\n';
            }
            {
                const std::size_t tp = carrier[0], fp = carrier[1], fn = carrier[2], tn = carrier[3];
                const double prec = (tp + fp) ? 100.0 * static_cast<double>(tp) / static_cast<double>(tp + fp) : 0.0;
                const double rec = (tp + fn) ? 100.0 * static_cast<double>(tp) / static_cast<double>(tp + fn) : 0.0;
                sum << "gt_carrier\tALL\tTP\t" << tp << '\t' << prec << '\n';
                sum << "gt_carrier\tALL\tFP\t" << fp << '\t' << prec << '\n';
                sum << "gt_carrier\tALL\tFN\t" << fn << '\t' << rec << '\n';
                sum << "gt_carrier\tALL\tTN\t" << tn << '\t' << rec << '\n';
            }
            {
                const double denom = static_cast<double>(gt_tot_ref_delta) - static_cast<double>(gt_tot_graph_delta);
                const bool has_gap = std::fabs(denom) > 1e-9;
                const double closed = has_gap
                    ? (static_cast<double>(gt_tot_ref_delta) - static_cast<double>(gt_tot_gt_delta)) / denom : 0.0;
                sum << "gt_gap\tALL\tbaseline_delta\t" << gt_tot_ref_delta << "\t0\n";
                sum << "gt_gap\tALL\tgenotype_delta\t" << gt_tot_gt_delta << "\t0\n";
                sum << "gt_gap\tALL\tgraph_delta\t" << gt_tot_graph_delta << "\t0\n";
                sum << "gt_gap\tALL\tgap_closed_pooled\t0\t";
                if (has_gap) sum << 100.0 * closed << '\n'; else sum << "NA\n";
                sum << "gt_gap\tALL\tgap_closed_mean\t" << gt_closed_n << '\t';
                if (gt_closed_n) sum << 100.0 * gt_closed_sum / static_cast<double>(gt_closed_n) << '\n';
                else sum << "NA\n";
                sum << "gt_gap\tALL\tgap_closed_undefined\t" << gt_closed_na << "\t0\n";
                // The same question against the other denominator: of the distance between the plain
                // reference and the truth, what fraction do the calls remove. gap_closed measures
                // against the ACHIEVABLE bound (the graph), this against ALL of it -- so a locus whose
                // graph cannot hold the haplotype scores low here and high there, and the pair says
                // which of the two is the limit. Undefined when there is no variation to recover.
                const double base = static_cast<double>(gt_tot_ref_delta);
                sum << "variation_recovered\tALL\tgenotype\t0\t";
                if (base > 0.0) sum << 100.0 * (base - static_cast<double>(gt_tot_gt_delta)) / base << '\n';
                else sum << "NA\n";
                sum << "variation_recovered\tALL\tcalled\t0\t";
                if (base > 0.0) sum << 100.0 * (base - static_cast<double>(tot_called_delta)) / base << '\n';
                else sum << "NA\n";
                sum << "variation_recovered\tALL\tgraph\t0\t";
                if (base > 0.0) sum << 100.0 * (base - static_cast<double>(gt_tot_graph_delta)) / base << '\n';
                else sum << "NA\n";
                sum << "gt_gap\tALL\tworse_than_baseline\t" << gt_worse_than_ref << "\t"
                    << (gt_scored_haps ? 100.0 * static_cast<double>(gt_worse_than_ref) / static_cast<double>(gt_scored_haps) : 0.0) << "\n";
            }
            // What the genotype QV was actually computed from, so it is never read without that context.
            sum << "gt_records\tALL\tapplied\t" << tot_gt.applied << "\t0\n";
            sum << "gt_records\tALL\tunplaceable\t" << tot_gt.unplaceable << "\t0\n";
            sum << "gt_records\tALL\tclamped\t" << tot_gt.clamped << "\t0\n";
            sum << "gt_records\tALL\tunhandled\t" << tot_gt.unhandled << "\t0\n";
            sum << "gt_records\tALL\tref_mismatch\t" << tot_gt.ref_mismatch << "\t0\n";
            sum << "gt_records\tALL\theuristic\t" << tot_gt.heuristic << "\t0\n";
        }
        for (const auto& [svclass, bands] : class_bands) {
            std::size_t total = 0;
            for (const auto& [b, n] : bands) total += n;
            for (const char* band : kBands) {
                const auto bit = bands.find(band);
                const std::size_t n = bit != bands.end() ? bit->second : 0;
                const double pct = total ? 100.0 * static_cast<double>(n) / static_cast<double>(total) : 0.0;
                sum << "svclass\t" << svclass << '\t' << band << '\t' << n << '\t' << pct << '\n';
            }
        }
    }
    staged.commit();

    auto quintile_recap = [&](const std::map<std::string, std::size_t>& counts, std::size_t denom) {
        std::string s;
        for (const char* q : kQuintiles) {
            const auto it = counts.find(q);
            const std::size_t n = it != counts.end() ? it->second : 0;
            const double pct = denom ? 100.0 * static_cast<double>(n) / static_cast<double>(denom) : 0.0;
            char buf[64];
            std::snprintf(buf, sizeof(buf), " %s=%.1f%%", q, pct);
            s += buf;
        }
        return s;
    };
    log.info("graph (optimistic upper bound) " + std::to_string(scored_haps) + " haplotypes, qv/qv_max quintiles:" +
             quintile_recap(quintile_count, scored_haps));
    log.info("called (retained calls only) " + std::to_string(scored_haps) + " haplotypes, qv/qv_max quintiles:" +
             quintile_recap(called_quintile_count, scored_haps));
    {
        char lb[256];
        std::snprintf(lb, sizeof(lb),
                      "truth events at >=%zu bp: %zu called, %zu MISSED (%zu bp); %zu below threshold",
                      min_sv_bp, tot_led.called, tot_led.missed, tot_led.missed_bp, tot_led.below);
        log.info(lb);
    }
    if (tot_coarse)
        log.info("note: " + std::to_string(tot_coarse) + " haplotype-bubble pair(s) were too large to "
                 "decompose and were treated as a single event");
    if (gt_scored_haps) {
        log.info("genotype (from VCF) " + std::to_string(gt_scored_haps) + " haplotypes, qv/qv_max quintiles:" +
                 quintile_recap(gt_quintile_count, gt_scored_haps));
        {
            const double denom = static_cast<double>(gt_tot_ref_delta) - static_cast<double>(gt_tot_graph_delta);
            const bool has_gap = std::fabs(denom) > 1e-9;
            const double closed = has_gap
                ? (static_cast<double>(gt_tot_ref_delta) - static_cast<double>(gt_tot_gt_delta)) / denom : 0.0;
            char gb[288];
            if (has_gap)
                std::snprintf(gb, sizeof(gb),
                              "residual bases: baseline(no edits)=%zu genotype=%zu graph=%zu -> gap closed %.1f%%"
                              " (%zu haplotypes worse than baseline)",
                              gt_tot_ref_delta, gt_tot_gt_delta, gt_tot_graph_delta, 100.0 * closed, gt_worse_than_ref);
            else
                std::snprintf(gb, sizeof(gb),
                              "residual bases: baseline(no edits)=%zu genotype=%zu graph=%zu -> gap closed UNDEFINED"
                              " (the graph bound equals the baseline, so there is no gap to close)",
                              gt_tot_ref_delta, gt_tot_gt_delta, gt_tot_graph_delta);
            log.info(gb);
        }
        const double prec = (carrier[0] + carrier[1]) ? 100.0 * static_cast<double>(carrier[0]) / static_cast<double>(carrier[0] + carrier[1]) : 0.0;
        const double rec = (carrier[0] + carrier[2]) ? 100.0 * static_cast<double>(carrier[0]) / static_cast<double>(carrier[0] + carrier[2]) : 0.0;
        char buf[192];
        std::snprintf(buf, sizeof(buf),
                      "carrier calls: TP=%zu FP=%zu FN=%zu TN=%zu (precision %.1f%%, recall %.1f%%)",
                      carrier[0], carrier[1], carrier[2], carrier[3], prec, rec);
        log.info(buf);
    } else if (do_gt) {
        log.info("genotype scoring produced no scored haplotype (no bubble placed on the reference)");
    }
    log.wrote(finals);
    log.done();
    return 0;
}

} // namespace panvar
