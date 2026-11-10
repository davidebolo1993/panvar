#include "panvar/recruit_command.hpp"

#include "panvar/background.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/kmer.hpp"
#include "panvar/locus_index.hpp"

#include <algorithm>
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
        << "  panvar recruit --locus <graph.gfa> [--locus ...] --background <store.bgk> -r <reference path> --index-only\n\n"
        << "Options:\n"
        << "      --locus <gfa>        A locus graph; repeatable\n"
        << "      --background <dir>   The background store (its .sshash beside it is used when current)\n"
        << "  -r, --reference <name>   The graph path carrying the coordinates (e.g. grch38)\n"
        << "      --coordinate-reference NAME  The background's coordinate reference (default: its only one)\n"
        << "      --index-only         Build or check the locus indexes, read nothing\n"
        << "      --save-index <dir>   Where indexes are kept and reused (default: beside each graph)\n"
        << "      --informative-max <N>  A k-mer seen up to N times outside the locus is shared, more is repeat (default 8)\n"
        << "      --pad <bp>           Decoy windows around homologous hits (default 1500)\n"
        << "      --dust <score>       Private k-mers this low-complexity become shared (default 2.5)\n"
        << "      --decoy-genomes A,B  Genomes whose homologous windows are decoys (default: the coordinate reference)\n"
        << "      --decoy-min-hits <N> Shared k-mer hits, each within a pad of the next, that make a homologous window (default 20)\n"
        << "  -t, --threads <N>        Worker threads\n"
        << "  -q, --quiet              No progress logs\n"
        << "  -h, --help               Show this help\n";
}

}  // namespace

int run_recruit_command(const std::vector<std::string>& args) {
    if (args.empty()) { print_recruit_help(); return 0; }
    recruit::BuildOptions o;
    std::string save_dir;
    bool index_only = false, quiet = false;
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
        else if (a == "--dust") o.dust = std::stod(value());
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
        } else if (a == "-t" || a == "--threads") threads = cli::parse_size_arg(a, value());
        else if (a == "-q" || a == "--quiet") quiet = true;
        else throw std::runtime_error("Unknown option: " + a);
    }
    if (o.loci.empty()) throw std::runtime_error("recruit requires --locus <graph.gfa>");
    if (o.background.empty()) throw std::runtime_error("recruit requires --background <store.bgk>");
    if (o.reference_query.empty()) throw std::runtime_error("recruit requires -r <reference path>");
    if (!index_only) throw std::runtime_error("recruit: reading reads is not available yet; run with --index-only");
    for (const auto& l : o.loci) if (!fs::is_regular_file(l)) throw std::runtime_error("recruit: no graph at " + l);
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
    log.wrote(wrote);
    log.done();
    return 0;
}

} // namespace panvar
