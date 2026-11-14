#include "panvar/transcripts_command.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/transcripts.hpp"

#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {
namespace {

void print_transcripts_help() {
    std::cout
        << "Usage:\n"
        << "  panvar transcripts -i <graph.gfa> -r <name> --gtf <genes.gtf> -o <prefix> [options]\n\n"
        << "Project a reference annotation onto every haplotype walk: one haplotype-specific transcript\n"
        << "per transcript, haplotype and copy of its gene.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                  Locus GFA with haplotype paths (P or W lines) (required)\n"
        << "  -r, --reference-path <name>       Reference path, exact name or unique substring; its name\n"
        << "                                    carries contig:start-end (required)\n"
        << "      --gtf <path>                  GTF in reference coordinates, plain or gzipped (required)\n"
        << "  -o, --out-prefix <prefix>         Output prefix (required)\n"
        << "      --min-exon-share <F>          Share of the gene's marker exons a pass must walk to be a\n"
        << "                                    copy (default 0.8)\n"
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
        tsv_ << "transcript_id\ttranscript_name\tgene_id\tgene_name\tbiotype\tstrand\thaplotype\tcopy\tcopies\tsource"
                "\tstatus\torientation\tpath_start\tpath_end\tlength\texons\tref_length\tref_exons\tsame_as_ref"
                "\tmoved_boundaries\tchanged_splice_sites\texon_starts\texon_lengths\tdetail\n";
    }

    void add(const TranscriptsResult& res, const std::vector<HstRecord>& recs) override {
        std::map<std::string, int> copies;   // per gene on this haplotype
        for (const HstRecord& r : recs) {
            const TranscriptFeature& t = res.transcripts[r.transcript];
            const std::string& hap = res.haplotypes[r.haplotype];
            const bool has = !r.sequence.empty();
            const bool same = has && r.sequence == res.reference_sequences[r.transcript];
            tsv_ << t.transcript_id << '\t' << t.transcript_name << '\t' << t.gene_id << '\t' << t.gene_name << '\t'
                 << (t.biotype.empty() ? "." : t.biotype) << '\t' << t.strand << '\t' << hap << '\t' << r.copy << '\t'
                 << r.copies << "\tprojected\t" << hst_status_name(r.status) << '\t' << orientation(r) << '\t';
            if (has) tsv_ << r.path_start << '\t' << r.path_end << '\t' << r.sequence.size() << '\t' << r.exons;
            else tsv_ << ".\t.\t.\t.";
            tsv_ << '\t' << res.reference_sequences[r.transcript].size() << '\t' << t.exons.size() << '\t'
                 << (has ? (same ? "yes" : "no") : ".") << '\t' << r.moved << '\t' << r.splice_changed << '\t'
                 << join(r.exon_starts) << '\t' << join(r.exon_lengths) << '\t' << (r.detail.empty() ? "." : r.detail)
                 << '\n';
            if (has)
                fa_ << '>' << t.transcript_id << '|' << t.gene_name << '|' << hap << '|' << r.copy
                    << " source=projected status=" << hst_status_name(r.status) << '\n' << r.sequence << '\n';
            Tally& g = genes_[t.gene_id];
            g.name = t.gene_name;
            ++g.status[static_cast<int>(r.status)];
            g.hst += has;
            g.same += same;
            copies[t.gene_id] = r.copies;
        }
        for (const auto& [gid, c] : copies) ++genes_[gid].by_copies[c];
    }

    void finish(const TranscriptsResult& res, const std::string& summary_path) {
        cli::close_output(fa_, fa_path_);
        cli::close_output(tsv_, tsv_path_);
        std::map<std::string, std::pair<std::string, std::size_t>> tx, out;
        for (const TranscriptFeature& t : res.transcripts) { tx[t.gene_id].first = t.gene_name; ++tx[t.gene_id].second; }
        for (const TranscriptFeature& t : res.partial) { out[t.gene_id].first = t.gene_name; ++out[t.gene_id].second; }
        std::vector<std::string> order;   // projected genes in reference order, then genes only partly inside
        for (const TranscriptFeature& t : res.transcripts)
            if (order.empty() || order.back() != t.gene_id) order.push_back(t.gene_id);
        for (const auto& [gid, v] : out) if (!tx.count(gid)) order.push_back(gid);
        std::ofstream s(summary_path);
        if (!s) throw std::runtime_error("transcripts: cannot write " + summary_path);
        s << "gene_id\tgene_name\ttranscripts\toutside\thaplotypes_by_copies\thsts\tsame_as_ref\tok"
             "\tsplice_site_changed\tboundary_moved\tmissing_exon\tno_copy\n";
        for (const std::string& gid : order) {
            const Tally& g = genes_[gid];
            const std::string name = tx.count(gid) ? tx[gid].first : out[gid].first;
            s << gid << '\t' << name << '\t' << (tx.count(gid) ? tx[gid].second : 0) << '\t'
              << (out.count(gid) ? out[gid].second : 0) << '\t';
            std::string bc;
            for (const auto& [c, n] : g.by_copies) bc += (bc.empty() ? "" : ",") + std::to_string(c) + ":" + std::to_string(n);
            s << (bc.empty() ? "." : bc) << '\t' << g.hst << '\t' << g.same;
            for (int k = 0; k < 5; ++k) s << '\t' << g.status[k];
            s << '\n';
            hsts_ += g.hst;
        }
        cli::close_output(s, summary_path);
    }

    std::size_t hsts() const { return hsts_; }

private:
    struct Tally {
        std::string name;
        std::size_t status[5] = {0, 0, 0, 0, 0};
        std::size_t hst = 0, same = 0;
        std::map<int, std::size_t> by_copies;   // copies -> haplotypes
    };
    std::ofstream fa_, tsv_;
    std::string fa_path_, tsv_path_;
    std::map<std::string, Tally> genes_;
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
        else if (arg == "-o" || arg == "--out-prefix") opt.out_prefix = require_value(arg);
        else if (arg == "--min-exon-share") opt.min_exon_share = cli::parse_unit_fraction_arg(arg, require_value(arg));
        else if (arg == "-t" || arg == "--threads") opt.threads = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "-q" || arg == "--quiet") opt.quiet = true;
        else throw std::runtime_error("Unknown option: " + arg);
    }
    if (opt.gfa_path.empty()) throw std::runtime_error("Missing required input: --gfa <path>");
    if (opt.reference_path.empty()) throw std::runtime_error("Missing required input: --reference-path <name>");
    if (opt.gtf_path.empty()) throw std::runtime_error("Missing required input: --gtf <path>");
    if (opt.out_prefix.empty()) throw std::runtime_error("Missing required input: --out-prefix <prefix>");
    const std::string fa = opt.out_prefix + ".hst.fa", tsv = opt.out_prefix + ".hst.tsv",
                      summary = opt.out_prefix + ".hst.summary.tsv";
    cli::reject_output_collisions("transcripts", {fa, tsv, summary}, {opt.gfa_path, opt.gtf_path});

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
    log.info("reference: " + res.reference);
    std::size_t genes = 0;
    for (std::size_t k = 0; k < res.transcripts.size(); ++k)
        genes += k == 0 || res.transcripts[k].gene_id != res.transcripts[k - 1].gene_id;
    log.info("transcripts: " + std::to_string(res.transcripts.size()) + " in " + std::to_string(genes) +
             " genes; reaching past the reference path, not projected: " + std::to_string(res.partial.size()));
    log.info("haplotype-specific transcripts: " + std::to_string(sink.hsts()));
    staged.commit();
    log.wrote({fa, tsv, summary});
    log.done();
    return 0;
}

}  // namespace panvar
