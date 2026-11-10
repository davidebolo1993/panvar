#include "panvar/background_command.hpp"

#include "panvar/background.hpp"
#include "panvar/cli_utils.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace panvar {

namespace {

namespace fs = std::filesystem;

void print_background_help() {
    std::cout
        << "panvar background - count every 31-mer of a set of genomes, for read recruitment\n\n"
        << "Usage:\n"
        << "  panvar background -o <store.bgk> --genome NAME=FASTA [--genome ...] [options]\n"
        << "  panvar background -o <store.bgk> --dump <out.tsv>\n\n"
        << "Options:\n"
        << "  -o, --out <dir>         The store (a folder)\n"
        << "      --genome NAME=PATH  A genome, FASTA plain or gzip/bgzip; repeatable\n"
        << "      --genomes <tsv>     Genomes as name and path per line\n"
        << "      --partitions <N>    Partitions, a power of two (default 4096)\n"
        << "      --memory <size>     Memory for counting, across threads, e.g. 2G, 512M (default 2G)\n"
        << "      --tmp-dir <dir>     Where the temporary spill files go (default: beside the store)\n"
        << "      --dump <path>       Write every k-mer and its count (all segments summed) as TSV\n"
        << "  -t, --threads <N>       Worker threads (default 0 = all cores)\n"
        << "  -q, --quiet             No progress logs\n"
        << "  -h, --help              Show this help\n";
}

std::size_t parse_memory(const std::string& v) {
    if (v.empty()) throw std::runtime_error("--memory is empty");
    std::size_t mult = 1;
    std::string digits = v;
    const char last = static_cast<char>(std::toupper(static_cast<unsigned char>(v.back())));
    if (last == 'K' || last == 'M' || last == 'G' || last == 'T') {
        digits.pop_back();
        mult = last == 'K' ? 1ULL << 10 : last == 'M' ? 1ULL << 20 : last == 'G' ? 1ULL << 30 : 1ULL << 40;
    }
    return cli::parse_size_arg("--memory", digits) * mult;
}

std::vector<std::string> segment_dirs(const fs::path& store) {
    std::vector<std::string> out;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(store, ec))
        if (e.is_directory() && e.path().filename().string().rfind("seg-", 0) == 0) out.push_back(e.path().string());
    std::sort(out.begin(), out.end());
    return out;
}

std::string json_escape(const std::string& s) {
    std::string o;
    for (char c : s) {
        if (c == '"' || c == '\\') o += '\\';
        o += c;
    }
    return o;
}

} // namespace

int run_background_command(const std::vector<std::string>& args) {
    if (args.empty()) { print_background_help(); return 0; }
    std::string out, dump, tmp_dir;
    std::vector<bg::Genome> genomes;
    std::size_t threads = 0, memory = std::size_t{2} << 30, partitions = 4096;
    bool quiet = false;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& a = args[i];
        auto value = [&]() -> const std::string& {
            if (i + 1 >= args.size()) throw std::runtime_error("Missing value after " + a);
            return args[++i];
        };
        if (a == "-h" || a == "--help") { print_background_help(); return 0; }
        else if (a == "-o" || a == "--out") out = value();
        else if (a == "--genome") {
            const std::string v = value();
            const auto eq = v.find('=');
            if (eq == std::string::npos || eq == 0 || eq + 1 == v.size())
                throw std::runtime_error("--genome takes NAME=PATH, got " + v);
            genomes.push_back({v.substr(0, eq), v.substr(eq + 1)});
        } else if (a == "--genomes") {
            const std::string path = value();
            std::ifstream in(path);
            if (!in) throw std::runtime_error("cannot read " + path);
            std::string line;
            while (std::getline(in, line)) {
                line = cli::trim_ascii(line);
                if (line.empty() || line[0] == '#') continue;
                const auto tab = line.find_first_of("\t ");
                if (tab == std::string::npos) throw std::runtime_error(path + ": expected name and path, got " + line);
                genomes.push_back({line.substr(0, tab), cli::trim_ascii(line.substr(tab + 1))});
            }
        } else if (a == "--partitions") partitions = cli::parse_size_arg(a, value());
        else if (a == "--memory") memory = parse_memory(value());
        else if (a == "--tmp-dir") tmp_dir = value();
        else if (a == "--dump") dump = value();
        else if (a == "-t" || a == "--threads") threads = cli::parse_size_arg(a, value());
        else if (a == "-q" || a == "--quiet") quiet = true;
        else throw std::runtime_error("Unknown option: " + a);
    }
    if (out.empty()) throw std::runtime_error("background requires -o/--out <store.bgk>");
    if (threads == 0) threads = std::max(1u, std::thread::hardware_concurrency());

    cli::RunLog log("background", quiet);
    if (!dump.empty()) {
        if (!genomes.empty()) throw std::runtime_error("--dump reads an existing store: give no --genome");
        const auto segs = segment_dirs(out);
        if (segs.empty()) throw std::runtime_error(out + " holds no segment");
        std::vector<std::unique_ptr<bg::Segment>> open;
        for (const auto& s : segs) open.push_back(std::make_unique<bg::Segment>(s));
        if (open.size() > 1) throw std::runtime_error("--dump of several segments comes with merge");
        std::FILE* f = std::fopen(dump.c_str(), "w");
        if (f == nullptr) throw std::runtime_error("cannot write " + dump);
        std::uint64_t n = 0;
        for (std::uint32_t p = 0; p < open[0]->partitions(); ++p)
            open[0]->for_each(p, [&](std::uint64_t code, std::uint64_t count) {
                std::fprintf(f, "%s\t%llu\n", kmer::decode(code, bg::kK).c_str(), static_cast<unsigned long long>(count));
                ++n;
            });
        if (std::fclose(f) != 0) throw std::runtime_error("cannot write " + dump);
        log.info("dumped " + std::to_string(n) + " k-mers");
        log.wrote({dump});
        log.done();
        return 0;
    }

    if (genomes.empty()) throw std::runtime_error("background: give --genome NAME=PATH (or --genomes)");
    if (fs::exists(fs::path(out) / "manifest.json"))
        throw std::runtime_error(out + " already holds a store");
    bg::BuildOptions o;
    o.genomes = genomes;
    o.segment_dir = (fs::path(out) / "seg-0").string();
    o.tmp_dir = tmp_dir.empty() ? out : tmp_dir;
    o.partitions = static_cast<std::uint32_t>(partitions);
    o.threads = threads;
    o.memory_bytes = memory;
    o.log = [&](const std::string& m) { log.info(m); };
    fs::create_directories(out);
    const bg::SegmentStats st = bg::build_segment(o);
    log.info("segment seg-0: " + std::to_string(st.distinct) + " distinct 31-mers of " + std::to_string(st.kmers) +
             " occurrences; " + std::to_string(st.bytes) + " bytes; spill " + std::to_string(st.spill_bytes) +
             " bytes in " + std::to_string(st.spill_files) + " files");

    const std::string manifest = (fs::path(out) / "manifest.json").string();
    std::ofstream m(manifest);
    m << "{\n  \"format\": \"panvar-background\",\n  \"version\": 1,\n  \"k\": " << bg::kK << ",\n  \"minimizer\": "
      << bg::kM << ",\n  \"partitions\": " << partitions << ",\n  \"segments\": [\n    {\"dir\": \"seg-0\", \"distinct\": "
      << st.distinct << ", \"occurrences\": " << st.kmers << ", \"genomes\": [";
    for (std::size_t g = 0; g < genomes.size(); ++g)
        m << (g ? ", " : "") << "{\"name\": \"" << json_escape(genomes[g].name) << "\", \"path\": \""
          << json_escape(fs::absolute(genomes[g].path).string()) << "\"}";
    m << "]}\n  ]\n}\n";
    if (!m) throw std::runtime_error("cannot write " + manifest);
    log.wrote({manifest, o.segment_dir});
    log.done();
    return 0;
}

} // namespace panvar
