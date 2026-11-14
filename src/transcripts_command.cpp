#include "panvar/transcripts_command.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/transcripts.hpp"

#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {
namespace {

void print_transcripts_help() {
    std::cout
        << "Usage:\n"
        << "  panvar transcripts -i <graph.gfa> -r <name> --gtf <genes.gtf> -o <prefix> [options]\n"
        << "  panvar transcripts -i <graph.gfa> --hap-annotations <table.tsv> -o <prefix> [options]\n\n"
        << "Project a reference annotation onto every haplotype walk: one haplotype-specific transcript\n"
        << "per transcript, haplotype and copy of its gene.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                  Locus GFA with haplotype paths (P or W lines) (required)\n"
        << "  -r, --reference-path <name>       Reference path, exact name or unique substring; its name\n"
        << "                                    carries contig:start-end (required with --gtf)\n"
        << "      --gtf <path>                  GTF in reference coordinates, plain or gzipped, to project\n"
        << "      --hap-annotations <tsv>       Haplotypes' own annotations: haplotype (sample#hap) and the\n"
        << "                                    path of its GFF3 or GTF, in the assembly's coordinates\n"
        << "      --contig-alias <tsv>          Contig names, one contig per row, matching the annotations'\n"
        << "                                    contigs to the paths' (default: exact, or without sample#hap#)\n"
        << "  -o, --out-prefix <prefix>         Output prefix (required)\n"
        << "      --min-exon-share <F>          Share of the gene's marker exons a pass must walk to be a\n"
        << "                                    copy (default 0.8)\n"
        << "      --min-identity <F>            Identity of the marker exons for a copy placed by sequence,\n"
        << "                                    where its nodes fall short of the share (default 0.93)\n"
        << "  -t, --threads <N>                 Worker threads (default 0 = all cores)\n"
        << "  -q, --quiet                       No progress logs\n"
        << "  -h, --help                        Show this help\n";
}

std::string join(const std::vector<std::size_t>& v) {
    std::string s;
    for (const std::size_t x : v) s += (s.empty() ? "" : ",") + std::to_string(x);
    return s.empty() ? "." : s;
}

std::string orientation(const HstRecord& r) { return r.status == HstStatus::NoCopy ? "." : r.reverse ? "-" : "+"; }

// Writes the FASTA and the table as haplotypes arrive, and tallies the per-gene summary.
class FileSink : public HstSink {
public:
    FileSink(const std::string& fa, const std::string& tsv) : fa_(fa), tsv_(tsv), fa_path_(fa), tsv_path_(tsv) {
        if (!fa_) throw std::runtime_error("transcripts: cannot write " + fa);
        if (!tsv_) throw std::runtime_error("transcripts: cannot write " + tsv);
        tsv_ << "transcript_id\ttranscript_name\tgene_id\tgene_name\tbiotype\tstrand\thaplotype\tcopy\tcopies\tsource\trecord"
                "\tstatus\torientation\tpath_start\tpath_end\tlength\texons\tref_length\tref_exons\tsame_as_ref"
                "\tmoved_boundaries\tchanged_splice_sites\texon_starts\texon_lengths\tdetail\n";
    }

    void add(const TranscriptsResult& res, const std::vector<HstRecord>& recs) override {
        if (!indexed_) {   // projected transcripts by id, and gene ids by name, for the annotated ones
            for (std::size_t k = 0; k < res.transcripts.size(); ++k) {
                ref_of_.emplace(res.transcripts[k].transcript_id, k);
                gene_of_.emplace(res.transcripts[k].gene_name, res.transcripts[k].gene_id);
                gene_ids_.insert(res.transcripts[k].gene_id);
            }
            for (const AnnotatedTranscript& at : res.annotated) annotated_genes_.insert(gene_id(at.t));
            indexed_ = true;
        }
        std::map<std::string, std::pair<int, int>> copies;   // per gene on this haplotype: projected, annotated
        for (const HstRecord& r : recs) {
            const bool annotated_only = r.source == kAnnotated;
            const TranscriptFeature& t = annotated_only ? res.annotated[r.transcript].t : res.transcripts[r.transcript];
            const auto ri = annotated_only ? ref_of_.find(t.transcript_id) : ref_of_.end();
            const long ref = !annotated_only ? static_cast<long>(r.transcript) : ri == ref_of_.end() ? -1 : static_cast<long>(ri->second);
            const std::string gid = annotated_only ? gene_id(t) : t.gene_id;
            const std::string& hap = res.haplotypes[r.haplotype];
            const bool has = !r.sequence.empty();
            const bool same = has && ref >= 0 && r.sequence == res.reference_sequences[static_cast<std::size_t>(ref)];
            const char* source = r.source == kProjected ? "projected" : annotated_only ? "annotated" : "projected+annotated";
            tsv_ << t.transcript_id << '\t' << t.transcript_name << '\t' << gid << '\t' << t.gene_name << '\t'
                 << (t.biotype.empty() ? "." : t.biotype) << '\t' << t.strand << '\t' << hap << '\t' << r.copy << '\t'
                 << r.copies << '\t' << source << '\t' << (r.record.empty() ? "." : r.record) << '\t'
                 << hst_status_name(r.status) << '\t'
                 << (annotated_only ? "." : orientation(r)) << '\t';
            if (has) tsv_ << r.path_start << '\t' << r.path_end << '\t' << r.sequence.size() << '\t' << r.exons;
            else tsv_ << ".\t.\t.\t.";
            if (ref >= 0)
                tsv_ << '\t' << res.reference_sequences[static_cast<std::size_t>(ref)].size() << '\t'
                     << res.transcripts[static_cast<std::size_t>(ref)].exons.size();
            else
                tsv_ << "\t.\t.";
            tsv_ << '\t' << (has && ref >= 0 ? (same ? "yes" : "no") : ".") << '\t' << r.moved << '\t' << r.splice_changed
                 << '\t' << join(r.exon_starts) << '\t' << join(r.exon_lengths) << '\t' << (r.detail.empty() ? "." : r.detail)
                 << '\n';
            if (has)
                fa_ << '>' << t.transcript_id << '|' << t.gene_name << '|' << hap << '|' << r.copy
                    << (annotated_only ? "|annotated:" + r.record : "") << " source=" << source << " status=" << hst_status_name(r.status)
                    << '\n' << r.sequence << '\n';
            Tally& g = genes_[gid];
            g.name = t.gene_name;
            ++g.status[static_cast<int>(r.status)];
            g.hst += has;
            g.same += same;
            g.annotated += has && (r.source & kAnnotated);
            auto& c = copies.emplace(gid, std::make_pair(0, 0)).first->second;
            if (r.source & kProjected) c.first = r.copies;
            if (r.annotated_copies >= 0) c.second = r.annotated_copies;
        }
        for (const auto& [gid, c] : copies) {
            ++genes_[gid].by_copies[c.first];
            // Only genes an annotation has somewhere: one it never annotates says nothing about copies.
            if (res.has_annotation[recs.front().haplotype] && c.first != c.second && annotated_genes_.count(gid)) {
                const std::string& hap = res.haplotypes[recs.front().haplotype];
                genes_[gid].differ.push_back(hap.substr(0, hap.find('#', hap.find('#') + 1)) + ":" +
                                             std::to_string(c.first) + "/" + std::to_string(c.second));
            }
        }
    }

    void finish(const TranscriptsResult& res, const std::string& summary_path) {
        cli::close_output(fa_, fa_path_);
        cli::close_output(tsv_, tsv_path_);
        std::map<std::string, std::pair<std::string, std::size_t>> tx, out;
        std::map<std::string, std::size_t> clipped;
        for (std::size_t k = 0; k < res.transcripts.size(); ++k) {
            const TranscriptFeature& t = res.transcripts[k];
            tx[t.gene_id].first = t.gene_name;
            ++tx[t.gene_id].second;
            clipped[t.gene_id] += res.clips[k].bp_5p + res.clips[k].bp_3p > 0;
        }
        for (const TranscriptFeature& t : res.outside) { out[t.gene_id].first = t.gene_name; ++out[t.gene_id].second; }
        std::vector<std::string> order;   // projected genes in reference order, then genes only partly inside
        for (const TranscriptFeature& t : res.transcripts)
            if (order.empty() || order.back() != t.gene_id) order.push_back(t.gene_id);
        for (const auto& [gid, v] : out) if (!tx.count(gid)) order.push_back(gid);
        for (const auto& [gid, g] : genes_)   // genes only the haplotypes' annotations have
            if (!tx.count(gid) && !out.count(gid)) { order.push_back(gid); tx[gid].first = g.name; }
        std::ofstream s(summary_path);
        if (!s) throw std::runtime_error("transcripts: cannot write " + summary_path);
        s << "gene_id\tgene_name\ttranscripts\tclipped\toutside\thaplotypes_by_copies\thsts\tsame_as_ref\tannotated";
        for (int k = 0; k < kHstStatuses; ++k) s << '\t' << hst_status_name(static_cast<HstStatus>(k));
        s << "\tcopies_differ\n";
        for (const std::string& gid : order) {
            const Tally& g = genes_[gid];
            const std::string name = tx.count(gid) ? tx[gid].first : out[gid].first;
            s << gid << '\t' << name << '\t' << (tx.count(gid) ? tx[gid].second : 0) << '\t' << clipped[gid] << '\t'
              << (out.count(gid) ? out[gid].second : 0) << '\t';
            std::string bc;
            for (const auto& [c, n] : g.by_copies) bc += (bc.empty() ? "" : ",") + std::to_string(c) + ":" + std::to_string(n);
            s << (bc.empty() ? "." : bc) << '\t' << g.hst << '\t' << g.same << '\t' << g.annotated;
            for (int k = 0; k < kHstStatuses; ++k) s << '\t' << g.status[k];
            std::string d;
            for (const std::string& x : g.differ) d += (d.empty() ? "" : ",") + x;
            s << '\t' << (d.empty() ? "." : d) << '\n';
            hsts_ += g.hst;
        }
        cli::close_output(s, summary_path);
    }

    std::size_t hsts() const { return hsts_; }

private:
    struct Tally {
        std::string name;
        std::size_t status[kHstStatuses] = {};
        std::size_t hst = 0, same = 0, annotated = 0;
        std::map<int, std::size_t> by_copies;   // projected copies -> haplotypes
        std::vector<std::string> differ;        // sample#hap:projected/annotated copies, where they differ
    };
    std::ofstream fa_, tsv_;
    std::string fa_path_, tsv_path_;
    std::map<std::string, Tally> genes_;
    bool indexed_ = false;
    std::unordered_map<std::string, std::size_t> ref_of_;
    std::unordered_map<std::string, std::string> gene_of_;
    std::set<std::string> annotated_genes_;   // gene ids any haplotype's annotation has
    std::set<std::string> gene_ids_;          // projected genes
    // An annotated transcript's gene, as a projected gene where its id or else its name is one.
    std::string gene_id(const TranscriptFeature& t) const {
        if (gene_ids_.count(t.gene_id)) return t.gene_id;
        const auto gi = gene_of_.find(t.gene_name);
        return gi != gene_of_.end() ? gi->second : t.gene_id;
    }
    std::size_t hsts_ = 0;
};

}  // namespace

int run_transcripts_command(const std::vector<std::string>& args) {
    if (args.empty()) {
        print_transcripts_help();
        return 0;
    }
    TranscriptsOptions opt;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) throw std::runtime_error("Missing value after " + flag);
            return args[++i];
        };
        if (arg == "-h" || arg == "--help") { print_transcripts_help(); return 0; }
        else if (arg == "-i" || arg == "--gfa") opt.gfa_path = require_value(arg);
        else if (arg == "-r" || arg == "--reference-path") opt.reference_path = require_value(arg);
        else if (arg == "--gtf") opt.gtf_path = require_value(arg);
        else if (arg == "--hap-annotations") opt.hap_annotations = require_value(arg);
        else if (arg == "--contig-alias") opt.contig_alias = require_value(arg);
        else if (arg == "-o" || arg == "--out-prefix") opt.out_prefix = require_value(arg);
        else if (arg == "--min-exon-share") opt.min_exon_share = cli::parse_unit_fraction_arg(arg, require_value(arg));
        else if (arg == "--min-identity") opt.min_identity = cli::parse_unit_fraction_arg(arg, require_value(arg));
        else if (arg == "-t" || arg == "--threads") opt.threads = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "-q" || arg == "--quiet") opt.quiet = true;
        else throw std::runtime_error("Unknown option: " + arg);
    }
    if (opt.gfa_path.empty()) throw std::runtime_error("Missing required input: --gfa <path>");
    if (opt.gtf_path.empty() && opt.hap_annotations.empty())
        throw std::runtime_error("Missing required input: --gtf <path> or --hap-annotations <tsv>");
    if (!opt.gtf_path.empty() && opt.reference_path.empty())
        throw std::runtime_error("Missing required input: --reference-path <name> (with --gtf)");
    if (opt.out_prefix.empty()) throw std::runtime_error("Missing required input: --out-prefix <prefix>");
    const std::string fa = opt.out_prefix + ".hst.fa", tsv = opt.out_prefix + ".hst.tsv",
                      summary = opt.out_prefix + ".hst.summary.tsv";
    cli::reject_output_collisions("transcripts", {fa, tsv, summary},
                                  {opt.gfa_path, opt.gtf_path, opt.hap_annotations, opt.contig_alias});

    cli::RunLog log("transcripts", opt.quiet);
    log.info("reading: " + opt.gfa_path);
    Graph graph = parse_gfa(opt.gfa_path);
    validate_graph_paths(graph, "transcripts", true, true);
    for (auto& [id, node] : graph.nodes) {   // soft-masked graphs spell the same transcripts
        (void)id;
        for (char& c : node.sequence) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    log.info("input: haplotypes: " + std::to_string(graph.paths.size()) + "; nodes: " +
             std::to_string(graph.nodes.size()));

    cli::StagedOutputs staged("transcripts");
    FileSink sink(staged.stage(fa), staged.stage(tsv));
    const TranscriptsResult res = project_transcripts(graph, opt, sink);
    sink.finish(res, staged.stage(summary));
    if (!res.reference.empty()) log.info("reference: " + res.reference);
    if (!res.not_given_back.empty()) {
        std::string ids;
        for (std::size_t k = 0; k < res.not_given_back.size() && k < 10; ++k) ids += (k ? "," : "") + res.not_given_back[k];
        std::cerr << cli::log_prefix("transcripts") << "warning: the reference walk does not give back "
                  << res.not_given_back.size() << " transcript(s) as annotated: " << ids
                  << (res.not_given_back.size() > 10 ? ",..." : "") << '\n';
    }
    if (!opt.hap_annotations.empty())
        log.info("haplotypes' own transcripts on their paths: " + std::to_string(res.annotated.size()));
    std::size_t genes = 0;
    for (std::size_t k = 0; k < res.transcripts.size(); ++k)
        genes += k == 0 || res.transcripts[k].gene_id != res.transcripts[k - 1].gene_id;
    std::size_t clipped = 0;
    for (const TranscriptsResult::EdgeClip& c : res.clips) clipped += c.bp_5p + c.bp_3p > 0;
    log.info("transcripts: " + std::to_string(res.transcripts.size()) + " in " + std::to_string(genes) +
             " genes, " + std::to_string(clipped) + " cut at the reference path's ends; with no exon on the path: " +
             std::to_string(res.outside.size()));
    log.info("haplotype-specific transcripts: " + std::to_string(sink.hsts()));
    staged.commit();
    log.wrote({fa, tsv, summary});
    log.done();
    return 0;
}

}  // namespace panvar
