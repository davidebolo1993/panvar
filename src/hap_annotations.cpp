#include "panvar/hap_annotations.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/gz_reader.hpp"
#include "panvar/ref_path.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace panvar {
namespace {

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::size_t a = 0;
    while (true) {
        const std::size_t b = s.find(sep, a);
        out.push_back(s.substr(a, b == std::string::npos ? std::string::npos : b - a));
        if (b == std::string::npos) return out;
        a = b + 1;
    }
}

// A name without its PanSN prefix sample#hap#.
std::string strip_pansn(const std::string& s) {
    const std::size_t a = s.find('#');
    const std::size_t b = a == std::string::npos ? a : s.find('#', a + 1);
    return b == std::string::npos ? s : s.substr(b + 1);
}

std::string strip_version(const std::string& s) {
    const std::size_t d = s.rfind('.');
    return d != std::string::npos && d + 1 < s.size() &&
                   s.find_first_not_of("0123456789", d + 1) == std::string::npos ? s.substr(0, d) : s;
}

using Attrs = std::map<std::string, std::string>;

// key=value;... (GFF3, values percent-encoded) or key "value"; ... (GTF).
Attrs parse_attrs(const std::string& col, bool gtf) {
    Attrs a;
    for (std::string item : split(col, ';')) {
        const std::size_t b = item.find_first_not_of(' ');
        if (b == std::string::npos) continue;
        item = item.substr(b);
        if (gtf) {
            const std::size_t sp = item.find(' ');
            if (sp == std::string::npos) continue;
            std::string v = item.substr(sp + 1);
            if (v.size() >= 2 && v.front() == '"') v = v.substr(1, v.rfind('"') - 1);
            a.emplace(item.substr(0, sp), v);   // the first of a repeated key (e.g. tag) is kept
        } else {
            const std::size_t eq = item.find('=');
            if (eq != std::string::npos) a.emplace(item.substr(0, eq), item.substr(eq + 1));
        }
    }
    return a;
}

std::string first_of(const Attrs& a, std::initializer_list<const char*> keys) {
    for (const char* k : keys) {
        const auto it = a.find(k);
        if (it != a.end() && !it->second.empty()) return it->second;
    }
    return std::string();
}

struct PathSpan {
    std::size_t path = 0;
    std::string contig;
    std::size_t first = 0, len = 0;   // the path's first base on its contig (1-based), and its length
};

struct Feature {
    std::string seqid;
    char strand = '+';
    Attrs attrs;
    std::string type;
    std::vector<std::string> parents;
};

}  // namespace

std::vector<AnnotatedTranscript> load_hap_annotations(const Graph& graph, const std::string& table,
                                                      const std::string& alias, const std::string& module) {
    // Contig aliases: every name on a row names one contig.
    std::unordered_map<std::string, std::size_t> group;
    if (!alias.empty()) {
        GzLineReader in(alias);
        if (!in.ok()) throw std::runtime_error(module + ": cannot read --contig-alias " + alias);
        std::string line;
        std::size_t row = 0;
        while (in.getline(line)) {
            if (line.empty() || line[0] == '#') continue;
            ++row;
            for (const std::string& n : split(line, '\t'))
                if (!n.empty()) group.emplace(n, row);
        }
    }
    const auto same_contig = [&](const std::string& seqid, const std::string& contig) {
        if (seqid == contig || strip_pansn(seqid) == contig) return true;
        for (const std::string& s : {seqid, strip_pansn(seqid)}) {
            const auto a = group.find(s), b = group.find(contig);
            if (a != group.end() && b != group.end() && a->second == b->second) return true;
        }
        return false;
    };

    std::ifstream tf(table);
    if (!tf) throw std::runtime_error(module + ": cannot read --hap-annotations " + table);
    std::vector<AnnotatedTranscript> out;
    std::string row;
    while (std::getline(tf, row)) {
        if (!row.empty() && row.back() == '\r') row.pop_back();
        const std::vector<std::string> f = split(row, '\t');
        if (row.empty() || row[0] == '#' || f[0] == "haplotype") continue;
        if (f.size() < 2) throw std::runtime_error(module + ": --hap-annotations row without a path: " + row);
        const std::string& hap = f[0];
        // The haplotype's paths and the contig interval each covers.
        std::vector<PathSpan> spans;
        for (std::size_t p = 0; p < graph.paths.size(); ++p) {
            const PathRecord& pr = graph.paths[p];
            if (pr.name.compare(0, hap.size() + 1, hap + "#") != 0) continue;
            const std::vector<std::size_t> pref = path_prefix_bp(pr, graph.nodes);
            const ParsedReferencePath meta = parse_reference_path_label(pr.name, pref.back(), pr.type == 'W');
            if (!meta.has_interval) {
                std::cerr << cli::log_prefix(module) << "warning: path " << pr.name
                          << " carries no contig:start-end; its annotation is not placed\n";
                continue;
            }
            const std::string tail = pr.name.substr(hap.size() + 1);
            spans.push_back({p, tail.substr(0, tail.rfind(':')), meta.region_start_1based, pref.back()});
        }
        if (spans.empty()) {
            std::cerr << cli::log_prefix(module) << "warning: --hap-annotations: no path of haplotype " << hap << "\n";
            continue;
        }
        std::unordered_map<std::string, std::vector<std::size_t>> on_contig;   // seqid -> spans
        const auto spans_of = [&](const std::string& seqid) -> const std::vector<std::size_t>& {
            auto it = on_contig.find(seqid);
            if (it == on_contig.end()) {
                std::vector<std::size_t> v;
                for (std::size_t k = 0; k < spans.size(); ++k)
                    if (same_contig(seqid, spans[k].contig)) v.push_back(k);
                it = on_contig.emplace(seqid, std::move(v)).first;
            }
            return it->second;
        };

        // Exons by transcript, and the transcripts' and genes' own records (GFF3).
        GzLineReader in(f[1]);
        if (!in.ok()) throw std::runtime_error(module + ": cannot read annotation " + f[1]);
        std::unordered_map<std::string, Feature> features;
        std::map<std::string, std::vector<std::pair<std::size_t, std::size_t>>> exons;
        std::unordered_map<std::string, Feature> exon_tx;   // GTF: the transcript, from its exon rows
        std::string line;
        while (in.getline(line)) {
            if (line.empty() || line[0] == '#') continue;
            const std::vector<std::string> c = split(line, '\t');
            if (c.size() < 9 || spans_of(c[0]).empty()) continue;
            const bool gtf = c[8].find('"') != std::string::npos && c[8].find('=') == std::string::npos;
            std::size_t s = 0, e = 0;
            try {
                s = std::stoull(c[3]);
                e = std::stoull(c[4]);
            } catch (const std::exception&) {
                continue;
            }
            Attrs a = parse_attrs(c[8], gtf);
            if (gtf) {
                if (c[2] != "exon") continue;
                // A transcript is its gene record's: copies of one transcript id in two gene records are two.
                const std::string tid = first_of(a, {"gene_id"}) + '\t' + first_of(a, {"transcript_id"});
                if (tid.size() == 1) continue;
                exons[tid].emplace_back(s, e);
                Feature& t = exon_tx[tid];
                if (t.seqid.empty()) t = {c[0], c[6] == "-" ? '-' : '+', std::move(a), "transcript", {}};
                continue;
            }
            if (c[2] == "exon") {
                for (const std::string& p : split(first_of(a, {"Parent"}), ','))
                    if (!p.empty()) exons[p].emplace_back(s, e);
                continue;
            }
            const std::string id = first_of(a, {"ID"});
            if (id.empty()) continue;
            Feature ft{c[0], c[6] == "-" ? '-' : '+', {}, c[2], split(first_of(a, {"Parent"}), ',')};
            ft.attrs = std::move(a);
            features.emplace(id, std::move(ft));
        }

        std::vector<AnnotatedTranscript> mine;
        for (auto& [tid, ex] : exons) {
            const auto fit = features.find(tid);
            const Feature* tx = fit != features.end() ? &fit->second : nullptr;
            if (!tx) {
                const auto git = exon_tx.find(tid);
                if (git == exon_tx.end()) continue;   // exons whose transcript record is elsewhere
                tx = &git->second;
            }
            const bool gtf = fit == features.end();
            const Feature* gene = nullptr;
            if (!gtf && !tx->parents.empty() && !tx->parents[0].empty()) {
                const auto g = features.find(tx->parents[0]);
                if (g != features.end()) gene = &g->second;
            }
            std::sort(ex.begin(), ex.end());
            ex.erase(std::unique(ex.begin(), ex.end()), ex.end());
            // The path holding most of the transcript's exonic bases.
            std::size_t best = SIZE_MAX, best_bp = 0;
            for (const std::size_t k : spans_of(tx->seqid)) {
                const std::size_t lo = spans[k].first, hi = spans[k].first + spans[k].len - 1;
                std::size_t bp = 0;
                for (const auto& [s, e] : ex)
                    if (e >= lo && s <= hi) bp += std::min(e, hi) - std::max(s, lo) + 1;
                if (bp > best_bp) { best = k; best_bp = bp; }
            }
            if (best == SIZE_MAX) continue;
            const PathSpan& ps = spans[best];
            AnnotatedTranscript at;
            at.path = ps.path;
            std::size_t left = 0, right = 0;
            for (const auto& [s, e] : ex) {
                const std::size_t lo = ps.first, hi = ps.first + ps.len - 1;
                if (s < lo) left += std::min(e + 1, lo) - s;
                if (e > hi) right += e - std::max(s, hi + 1) + 1;
                if (e < lo || s > hi) continue;
                at.t.exons.emplace_back(std::max(s, lo) - ps.first, std::min(e, hi) - ps.first + 1);
            }
            at.t.strand = tx->strand;
            at.bp_5p = tx->strand == '-' ? right : left;
            at.bp_3p = tx->strand == '-' ? left : right;
            const Attrs& ta = tx->attrs;
            const Attrs none;
            const Attrs& ga = gene ? gene->attrs : none;
            const std::string src = first_of(ta, {"source_transcript"});
            at.t.transcript_id = src.rfind("ENS", 0) == 0 ? strip_version(src) : first_of(ta, {"transcript_id", "ID"});
            if (at.t.transcript_id.empty()) at.t.transcript_id = tid;
            at.t.transcript_name = first_of(ta, {"source_transcript_name", "transcript_name"});
            if (at.t.transcript_name.empty()) at.t.transcript_name = at.t.transcript_id;
            at.t.gene_name = first_of(ga, {"gene_name", "source_gene_common_name", "gene", "Name"});
            if (at.t.gene_name.empty()) at.t.gene_name = first_of(ta, {"gene_name", "source_gene_common_name", "gene"});
            const std::string sg = first_of(ga, {"source_gene"}).empty() ? first_of(ta, {"source_gene"}) : first_of(ga, {"source_gene"});
            at.t.gene_id = sg.rfind("ENS", 0) == 0 ? strip_version(sg) : first_of(ga, {"gene_id", "ID"});
            if (at.t.gene_id.empty()) at.t.gene_id = first_of(ta, {"gene_id"});
            if (at.t.gene_name.empty()) at.t.gene_name = at.t.gene_id.empty() ? at.t.transcript_id : at.t.gene_id;
            if (at.t.gene_id.empty()) at.t.gene_id = at.t.gene_name;
            at.t.biotype = first_of(ta, {"transcript_biotype", "transcript_type", "biotype"});
            if (at.t.biotype.empty()) at.t.biotype = tx->type;
            at.gene_record = gtf ? first_of(ta, {"gene_id"}) : gene ? tx->parents[0] : tid;
            for (const char* k : {"source_transcript", "source_gene", "transcript_class", "alignment_mode"}) {
                const std::string v = first_of(ta, {k});
                if (!v.empty()) at.detail += (at.detail.empty() ? "" : ";") + std::string(k) + ":" + v;
            }
            if (!at.t.exons.empty()) mine.push_back(std::move(at));
        }
        // Copies: each gene's records on a path, ranked by their first exon along the path.
        std::map<std::pair<std::size_t, std::string>, std::map<std::string, std::size_t>> records;
        for (const AnnotatedTranscript& at : mine) {
            std::size_t& first = records[{at.path, at.t.gene_name}].emplace(at.gene_record, SIZE_MAX).first->second;
            first = std::min(first, at.t.exons.front().first);
        }
        for (AnnotatedTranscript& at : mine) {
            const auto& recs = records[{at.path, at.t.gene_name}];
            std::vector<std::pair<std::size_t, std::string>> order;
            for (const auto& [rec, first] : recs) order.emplace_back(first, rec);
            std::sort(order.begin(), order.end());
            at.copies = static_cast<int>(order.size());
            for (std::size_t k = 0; k < order.size(); ++k)
                if (order[k].second == at.gene_record) at.copy = static_cast<int>(k + 1);
        }
        out.insert(out.end(), std::make_move_iterator(mine.begin()), std::make_move_iterator(mine.end()));
    }
    return out;
}

}  // namespace panvar
