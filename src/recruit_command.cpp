#include "panvar/recruit_command.hpp"

#include "panvar/background.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/kmer.hpp"
#include "panvar/locus_index.hpp"
#include "panvar/read_input.hpp"
#include "panvar/recruit_reads.hpp"

#include <algorithm>
#include <chrono>
#include <map>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace panvar {

namespace {

namespace fs = std::filesystem;

void print_recruit_help() {
    std::cout
        << "panvar recruit - recruit a sample's reads to locus graphs without aligning\n\n"
        << "Usage:\n"
        << "  panvar recruit --locus <graph.gfa> [--locus ...] --background <store.bgk> -r <reference path> --index-only\n"
        << "  panvar recruit --locus <graph.gfa | index.plx> [...] [--background ... -r ...] <reads> -o <out.tsv>\n"
        << "  panvar recruit --locus ... --samples <sheet.tsv> -o <out dir>\n"
        << "  panvar recruit <reads> --dry-run\n\n"
        << "Reads: -1/-2 FASTQ pairs, -i any one file (FASTQ, SAM, BAM, CRAM), each repeatable; one sample's files\n"
        << "are read in turn. Secondary, supplementary and QC-failed records are skipped, unmapped and duplicate\n"
        << "ones kept; mates are paired in any sort order (coordinate-sorted: mates far apart are fetched at the\n"
        << "end, through the index when there is one).\n\n"
        << "Options:\n"
        << "      --locus <gfa>        A locus graph; repeatable\n"
        << "      --background <dir>   The background store (its .sshash beside it is used when current)\n"
        << "  -r, --reference <name>   The graph path carrying the coordinates (e.g. grch38)\n"
        << "      --coordinate-reference NAME  The background's coordinate reference (default: its only one)\n"
        << "      --index-only         Build or check the locus indexes, read nothing\n"
        << "      --save-index <dir>   Where indexes are kept and reused (default: beside each graph)\n"
        << "      --informative-max <N>  A k-mer seen up to N times per background haplotype outside the locus is\n"
        << "                           shared, more is repeat (default 8)\n"
        << "      --pad <bp>           Decoy windows around homologous hits (default 1500)\n"
        << "      --dust <score>       Private k-mers this low-complexity become shared (default 2.5)\n"
        << "      --decoy-genomes A,B  Genomes whose homologous windows are decoys (default: the coordinate reference)\n"
        << "      --decoy-min-hits <N> Shared k-mer hits, each within a pad of the next, that make a homologous window (default 20)\n"
        << "      --no-calibration     Build indexes without calibration (reads simulated from held-out samples,\n"
        << "                           decoy windows and the reference's 5 Mb flanks, scored at each site cut)\n"
        << "  -1, --reads <fq>         Reads (FASTQ, plain or gzip); repeatable\n"
        << "  -2, --mates <fq>         Their mates, in the same order (one per -1)\n"
        << "  -i, --input <file>       One file of reads: FASTQ (single-end, or --interleaved), SAM, BAM or CRAM;\n"
        << "                           repeatable\n"
        << "      --interleaved        Single FASTQ files hold mates as consecutive records\n"
        << "      --cram-reference <fa>  The FASTA (with .fai) CRAM files were compressed against\n"
        << "      --samples <tsv>      Several samples: lines of sample, file[, mate file] (a sample may take\n"
        << "                           several lines); -o is then a directory with <sample>.tsv per sample\n"
        << "      --dry-run            Count each input's records by kind and the time to read them; nothing else\n"
        << "  -o, --out <tsv>          One row per screened fragment and locus: fate (certain, ambiguous,\n"
        << "                           unrecruited), reason, each mate's locus segment (read bases) and the\n"
        << "                           evidence. A certain fragment is partial when a read runs out of the\n"
        << "                           locus, rescued when one mate alone decides the pair\n"
        << "      --method <m>         runs (default: >= --min-sites independent private sites, decoy sites at\n"
        << "                           most --decoy-share of them, at least --min-density of them),\n"
        << "                           kmers (>= 3 private k-mers, no decoy k-mer), fraction (>= --min-fraction\n"
        << "                           of k-mers in the locus, one private)\n"
        << "      --false-certain-rate <r>  Per locus, the lowest private-site cut (at least 2) whose simulated\n"
        << "                           off-target reads of the run's read class are certain at most at this\n"
        << "                           rate (default 1e-6; an index without calibration uses 2)\n"
        << "      --min-sites <N>      One private-site cut for every locus instead (--method runs)\n"
        << "      --min-fraction <f>   For --method fraction (default 0.5)\n"
        << "      --decoy-share <f>    Decoy sites tolerated as a share of the private sites (default 0.05)\n"
        << "      --min-density <f>    Private sites needed per 31 callable bases (default 0.05): a long read\n"
        << "                           elsewhere can carry a few private k-mers by sequencing error\n"
        << "      --error-rate <e>     The reads' error rate per base: more mismatch runs than errors and new\n"
        << "                           alleles explain make a read divergent (default: estimated from the base\n"
        << "                           qualities of the first reads; logged and in the summary)\n"
        << "      --preset <p>         The error rate of a technology: illumina 0.002, hifi 0.001, ont 0.03,\n"
        << "                           clr 0.1 (--error-rate wins)\n"
        << "      --reads-out <dir>    Write the recruited reads as gzip FASTQ, per locus and fate\n"
        << "                           (<locus>.certain_1/_2.fq.gz for pairs, <locus>.certain.fq.gz for single\n"
        << "                           reads; likewise ambiguous), as sequenced, with their qualities\n"
        << "      --certain-only       With --reads-out: certain reads only\n"
        << "      --summary <tsv>      The run's summary (default: beside -o, .summary.tsv): input records by\n"
        << "                           kind, fragments by locus, fate and reason, error rate, backbone depth\n"
        << "  -t, --threads <N>        Scoring threads; a reader and N/4 decompression threads run beside them\n"
        << "  -q, --quiet              No progress logs\n"
        << "  -h, --help               Show this help\n";
}

double parse_number(const std::string& name, const std::string& v) {
    std::size_t used = 0;
    double x = 0;
    try { x = std::stod(v, &used); } catch (const std::exception&) { used = 0; }
    if (used == 0 || used != v.size()) throw std::runtime_error(name + " takes a number, got '" + v + "'");
    return x;
}

}  // namespace

int run_recruit_command(const std::vector<std::string>& args) {
    if (args.empty()) { print_recruit_help(); return 0; }
    recruit::BuildOptions o;
    recruit::ReadOptions ro;
    std::string save_dir, sheet, summary;
    double error_rate = -1, preset_rate = -1;
    std::vector<std::string> reads1, reads2, singles;
    bool index_only = false, quiet = false, dry_run = false, rate_given = false;
    std::size_t threads = 0;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& a = args[i];
        auto value = [&]() -> const std::string& {
            if (i + 1 >= args.size()) throw std::runtime_error("Missing value after " + a);
            return args[++i];
        };
        if (a == "-h" || a == "--help") { print_recruit_help(); return 0; }
        else if (a == "--locus") o.loci.push_back(value());
        else if (a == "--background") o.background = value();
        else if (a == "-r" || a == "--reference") o.reference_query = value();
        else if (a == "--coordinate-reference") o.coordinate_reference = value();
        else if (a == "--index-only") index_only = true;
        else if (a == "--save-index") save_dir = value();
        else if (a == "--informative-max") o.informative_max = static_cast<std::uint32_t>(cli::parse_size_arg(a, value()));
        else if (a == "--pad") o.pad = cli::parse_size_arg(a, value());
        else if (a == "--dust") o.dust = parse_number(a, value());
        else if (a == "--no-calibration") o.calibrate = false;
        else if (a == "--false-certain-rate") { ro.false_certain_rate = parse_number(a, value()); rate_given = true; }
        else if (a == "--decoy-min-hits") o.decoy_min_hits = static_cast<std::uint32_t>(cli::parse_size_arg(a, value()));
        else if (a == "--decoy-genomes") {
            std::stringstream ss(value());
            for (std::string d; std::getline(ss, d, ',');) if (!d.empty()) o.decoy_genomes.push_back(d);
        } else if (a == "--dump-index") {   // PLX=TSV: every k-mer of a locus index, for checks
            const std::string v = value();
            const auto eq = v.find('=');
            if (eq == std::string::npos) throw std::runtime_error("--dump-index takes PLX=TSV");
            const recruit::LocusIndex x = recruit::load(v.substr(0, eq));
            std::ofstream out(v.substr(eq + 1));
            out << "kmer\tclass\toutside\tmultiplicity\tbackbone\tscreen\tnodes\n";
            for (std::size_t k = 0; k < x.codes.size(); ++k) {
                const bool screened = std::binary_search(x.screen.begin(), x.screen.end(), x.codes[k]);
                out << kmer::decode(x.codes[k], 31) << '\t' << "PSR"[x.cls[k]] << '\t' << x.outside[k] << '\t' << x.mult[k] << '\t'
                    << int(x.backbone[k]) << '\t' << int(screened) << '\t';
                for (std::uint32_t i = x.nodeset_first[x.nodeset[k]]; i < x.nodeset_first[x.nodeset[k] + 1]; ++i)
                    out << (i > x.nodeset_first[x.nodeset[k]] ? "," : "") << x.nodes[x.nodeset_nodes[i]];
                out << '\n';
            }
            for (const std::uint64_t d : x.decoy) out << kmer::decode(d, 31) << "\tD\t.\t.\t.\t0\t.\n";
            if (!out) throw std::runtime_error("recruit: cannot write " + v.substr(eq + 1));
            return 0;
        } else if (a == "-1" || a == "--reads") reads1.push_back(value());
        else if (a == "-2" || a == "--mates") reads2.push_back(value());
        else if (a == "-i" || a == "--input") singles.push_back(value());
        else if (a == "--interleaved") ro.interleaved = true;
        else if (a == "--cram-reference") ro.cram_reference = value();
        else if (a == "--samples") sheet = value();
        else if (a == "--dry-run") dry_run = true;
        else if (a == "-o" || a == "--out") ro.out = value();
        else if (a == "--method") {
            const std::string m = value();
            if (m == "runs") ro.method = recruit::Method::Runs;
            else if (m == "kmers") ro.method = recruit::Method::Kmers;
            else if (m == "fraction") ro.method = recruit::Method::Fraction;
            else throw std::runtime_error("--method is runs, kmers or fraction");
        } else if (a == "--min-sites") { ro.min_sites = static_cast<std::uint32_t>(cli::parse_size_arg(a, value())); ro.calibrated_cut = false; }
        else if (a == "--min-fraction") ro.min_fraction = parse_number(a, value());
        else if (a == "--decoy-share") ro.decoy_share = parse_number(a, value());
        else if (a == "--min-density") ro.min_density = parse_number(a, value());
        else if (a == "--error-rate") error_rate = parse_number(a, value());
        else if (a == "--preset") {
            const std::string p = value();
            if (p == "illumina") preset_rate = 0.002;
            else if (p == "hifi") preset_rate = 0.001;
            else if (p == "ont") preset_rate = 0.03;
            else if (p == "clr") preset_rate = 0.1;
            else throw std::runtime_error("--preset is illumina, hifi, ont or clr");
        } else if (a == "--reads-out") ro.reads_out = value();
        else if (a == "--certain-only") ro.reads_ambiguous = false;
        else if (a == "--summary") summary = value();
        else if (a == "-t" || a == "--threads") threads = cli::parse_size_arg(a, value());
        else if (a == "-q" || a == "--quiet") quiet = true;
        else throw std::runtime_error("Unknown option: " + a);
    }
    // one sample's inputs, or the sample sheet's samples in order
    if (!reads2.empty() && reads2.size() != reads1.size()) throw std::runtime_error("recruit: give one -2 per -1");
    std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string>>>> samples;
    if (!sheet.empty()) {
        if (!reads1.empty() || !singles.empty()) throw std::runtime_error("recruit: give reads by --samples or by -1/-2/-i, not both");
        std::ifstream in(sheet);
        if (!in) throw std::runtime_error("recruit: cannot read " + sheet);
        std::map<std::string, std::size_t> at;
        const fs::path base = fs::path(sheet).parent_path();
        auto resolve = [&](const std::string& p) { return p.empty() || fs::path(p).is_absolute() ? p : (base / p).string(); };
        for (std::string line; std::getline(in, line);) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty() || line[0] == '#') continue;
            std::vector<std::string> f;
            std::stringstream ss(line);
            for (std::string x; std::getline(ss, x, '\t');) f.push_back(x);
            if (f.size() < 2 || f.size() > 3 || f[0].empty() || f[1].empty())
                throw std::runtime_error("recruit: " + sheet + ": a line is not sample, file[, mate file]: " + line);
            if (!at.count(f[0])) { at[f[0]] = samples.size(); samples.push_back({f[0], {}}); }
            samples[at[f[0]]].second.emplace_back(resolve(f[1]), f.size() == 3 ? resolve(f[2]) : std::string());
        }
        if (samples.empty()) throw std::runtime_error("recruit: " + sheet + " lists no sample");
    } else if (!reads1.empty() || !singles.empty()) {
        std::vector<std::pair<std::string, std::string>> files;
        for (std::size_t i = 0; i < reads1.size(); ++i) files.emplace_back(reads1[i], reads2.empty() ? std::string() : reads2[i]);
        for (const auto& f : singles) files.emplace_back(f, std::string());
        samples.push_back({"", files});
    }
    if (dry_run) {
        if (samples.empty()) throw std::runtime_error("recruit --dry-run requires reads (-1/-2, -i or --samples)");
        recruit::InputOptions io;
        io.cram_reference = ro.cram_reference;
        io.threads = threads == 0 ? bg::default_threads() : threads;
        io.interleaved = ro.interleaved;
        std::cout << "sample\tfile\trecords\tsecondary\tsupplementary\tqc_fail\tduplicates\tunmapped\tsingle_reads\tpaired_reads\tseconds\n";
        for (const auto& [name, files] : samples)
            for (const auto& [a, b] : files)
                for (const std::string* f : {&a, &b}) {
                    if (f->empty()) continue;
                    const auto t0 = std::chrono::steady_clock::now();
                    const recruit::InputStats s = recruit::count_records(*f, io);
                    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                    std::cout << (name.empty() ? "." : name) << '\t' << *f << '\t' << s.records << '\t' << s.secondary << '\t' << s.supplementary
                              << '\t' << s.qc_fail << '\t' << s.duplicates << '\t' << s.unmapped << '\t' << s.fragments << '\t' << s.pairs
                              << '\t' << sec << '\n';
                }
        return 0;
    }
    if (o.loci.empty()) throw std::runtime_error("recruit requires --locus <graph.gfa | index.plx>");
    auto is_plx = [](const std::string& p) { return p.size() > 4 && p.compare(p.size() - 4, 4, ".plx") == 0; };
    std::vector<std::string> given_indexes, graphs;
    for (const auto& l : o.loci) {
        if (!fs::is_regular_file(l)) throw std::runtime_error("recruit: no file at " + l);
        (is_plx(l) ? given_indexes : graphs).push_back(l);
    }
    o.loci = graphs;
    for (const auto& plx : given_indexes) {   // a given index is checked now, not first when reads arrive
        try {
            (void)recruit::load(plx);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string(e.what()) + " (rebuild it from its graph)");
        }
    }
    if (!graphs.empty() && o.background.empty()) throw std::runtime_error("recruit requires --background <store.bgk> to index a graph");
    if (!graphs.empty() && o.reference_query.empty()) throw std::runtime_error("recruit requires -r <reference path> to index a graph");
    if (!index_only && (samples.empty() || ro.out.empty()))
        throw std::runtime_error("recruit requires reads (-1/-2, -i or --samples) and -o (or --index-only)");
    if (ro.min_sites == 0) throw std::runtime_error("--min-sites must be at least 1");
    if (rate_given && !ro.calibrated_cut) throw std::runtime_error("give --false-certain-rate or --min-sites, not both");
    if (!(ro.false_certain_rate >= 0 && ro.false_certain_rate < 1)) throw std::runtime_error("--false-certain-rate must be in [0, 1)");
    if (ro.decoy_share < 0 || ro.decoy_share > 1) throw std::runtime_error("--decoy-share must be between 0 and 1");
    if (ro.min_density < 0 || ro.min_density > 1) throw std::runtime_error("--min-density must be between 0 and 1");
    if (error_rate >= 0.5) throw std::runtime_error("--error-rate must be between 0 and 0.5");
    ro.error_rate = error_rate >= 0 ? error_rate : preset_rate;   // below 0: estimated from base qualities
    if (!(ro.min_fraction > 0 && ro.min_fraction <= 1)) throw std::runtime_error("--min-fraction must be in (0, 1]");
    cli::RunLog log("recruit", quiet);
    o.threads = threads == 0 ? bg::default_threads() : threads;
    o.log = [&](const std::string& s) { log.info(s); };

    auto index_path = [&](const std::string& gfa, const std::string& ext) {
        std::string n = fs::path(gfa).filename().string();
        for (const char* e : {".gz", ".gfa"}) {
            const std::string x(e);
            if (n.size() > x.size() && n.compare(n.size() - x.size(), x.size(), x) == 0) n.resize(n.size() - x.size());
        }
        const fs::path dir = save_dir.empty() ? fs::path(gfa).parent_path() : fs::path(save_dir);
        fs::create_directories(dir.empty() ? fs::path(".") : dir);
        return (dir / (n + ext)).string();
    };
    {   // each locus writes <name>.plx: two graphs with one name would overwrite each other
        std::map<std::string, std::string> taken;
        for (const auto& l : o.loci) {
            const std::string plx = fs::absolute(index_path(l, ".plx")).lexically_normal().string();
            const auto [it, fresh] = taken.emplace(plx, l);
            if (!fresh) throw std::runtime_error("recruit: " + it->second + " and " + l + " would both be indexed as " + plx + "; rename one graph");
        }
    }
    recruit::BuildOptions todo = o;
    todo.loci.clear();
    std::vector<std::string> wrote;
    for (const auto& l : o.loci) {
        const std::string plx = index_path(l, ".plx");
        if (fs::exists(plx)) {
            try {
                const recruit::LocusIndex old = recruit::load(plx);
                if (old.key() == recruit::expected_key(o, l)) { log.info(l + ": index " + plx + " is current, reused"); continue; }
                log.info(l + ": index " + plx + " was built from another graph, background or settings; rebuilding");
            } catch (const std::exception& e) {
                log.info(l + ": " + e.what() + "; rebuilding");
            }
        }
        todo.loci.push_back(l);
    }
    if (!todo.loci.empty()) {
        const auto built = recruit::build(todo);
        for (std::size_t i = 0; i < built.size(); ++i) {
            const std::string plx = index_path(todo.loci[i], ".plx"), tsv = index_path(todo.loci[i], ".index.tsv");
            recruit::save(built[i], plx);
            recruit::write_report(built[i], tsv);
            wrote.push_back(plx);
            wrote.push_back(tsv);
        }
    }
    if (!index_only) {
        ro.indexes = given_indexes;
        for (const auto& l : o.loci) ro.indexes.push_back(index_path(l, ".plx"));
        ro.threads = o.threads;
        ro.log = o.log;
        const std::string out = ro.out;
        if (!sheet.empty()) fs::create_directories(out);
        const std::string reads_out = ro.reads_out;
        auto beside = [](const std::string& table) {   // x.tsv -> x.summary.tsv
            const std::string stem = table.size() > 4 && table.compare(table.size() - 4, 4, ".tsv") == 0 ? table.substr(0, table.size() - 4) : table;
            return stem + ".summary.tsv";
        };
        for (const auto& [name, files] : samples) {
            ro.inputs = files;
            ro.out = sheet.empty() ? out : (fs::path(out) / (name + ".tsv")).string();
            ro.summary = !summary.empty() && sheet.empty() ? summary : beside(ro.out);
            if (!reads_out.empty()) ro.reads_out = sheet.empty() ? reads_out : (fs::path(reads_out) / name).string();
            if (!name.empty()) log.info(name + ": " + std::to_string(files.size()) + " input(s)");
            const recruit::ReadStats st = recruit::recruit_reads(ro);
            log.info((name.empty() ? std::string() : name + ": ") + std::to_string(st.fragments) + " fragments, " + std::to_string(st.screened) +
                     " screened; per locus: " + std::to_string(st.certain) + " certain, " + std::to_string(st.ambiguous) + " ambiguous, " +
                     std::to_string(st.unrecruited) + " unrecruited");
            wrote.push_back(ro.out);
            wrote.push_back(ro.summary);
        }
    }
    log.wrote(wrote);
    log.done();
    return 0;
}

} // namespace panvar
