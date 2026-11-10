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
#include <set>
#include <sstream>
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
        << "  panvar background -o <store.bgk> --genome NAME=FASTA [--genome ...] [options]   build\n"
        << "  panvar background -o <store.bgk> --add --genome NAME=FASTA [...]                add a segment\n"
        << "  panvar background -o <store.bgk> --merge                                        merge segments\n"
        << "  panvar background -o <store.bgk> --dump <out.tsv>                               every k-mer\n\n"
        << "Options:\n"
        << "  -o, --out <dir>         The store (a folder of segments)\n"
        << "      --genome NAME=PATH  A genome, FASTA plain or gzip/bgzip; repeatable\n"
        << "      --genomes <tsv>     Genomes as name and path per line\n"
        << "      --add               Add the genomes to an existing store as a new segment\n"
        << "      --merge             Merge the store's segments into one\n"
        << "      --dump <path>       Write every k-mer and its count (segments summed) as TSV\n"
        << "      --partitions <N>    Partitions of a new store, a power of two (default 4096)\n"
        << "      --memory <size>     Memory for counting, across threads, e.g. 2G, 512M (default 2G)\n"
        << "      --tmp-dir <dir>     Where the temporary spill files go (default: inside the store)\n"
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

std::string json_escape(const std::string& s) {
    std::string o;
    for (char c : s) {
        if (c == '"' || c == '\\') o += '\\';
        o += c;
    }
    return o;
}

// A segment's own record: its genomes and counts, one key per line.
struct SegmentInfo {
    std::vector<bg::Genome> genomes;
    std::uint64_t distinct = 0, occurrences = 0;
};

void write_segment_info(const std::string& dir, const SegmentInfo& info) {
    const std::string path = (fs::path(dir) / "segment.tsv").string();
    std::ofstream f(path);
    for (const auto& g : info.genomes) f << "genome\t" << g.name << '\t' << g.path << '\n';
    f << "distinct\t" << info.distinct << "\noccurrences\t" << info.occurrences << '\n';
    if (!f) throw std::runtime_error("cannot write " + path);
}

SegmentInfo read_segment_info(const std::string& dir) {
    SegmentInfo info;
    std::ifstream f(fs::path(dir) / "segment.tsv");
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        std::string key;
        std::getline(ss, key, '\t');
        if (key == "genome") {
            bg::Genome g;
            std::getline(ss, g.name, '\t');
            std::getline(ss, g.path);
            info.genomes.push_back(g);
        } else if (key == "distinct") ss >> info.distinct;
        else if (key == "occurrences") ss >> info.occurrences;
    }
    return info;
}

// manifest.json, regenerated from the segment folders.
void write_manifest(const std::string& store, std::uint32_t partitions) {
    const std::string path = (fs::path(store) / "manifest.json").string();
    std::ofstream m(path);
    m << "{\n  \"format\": \"panvar-background\",\n  \"version\": 1,\n  \"k\": " << bg::kK << ",\n  \"minimizer\": "
      << bg::kM << ",\n  \"partitions\": " << partitions << ",\n  \"segments\": [";
    const auto dirs = bg::segment_dirs_of(store);
    for (std::size_t s = 0; s < dirs.size(); ++s) {
        const SegmentInfo info = read_segment_info(dirs[s]);
        m << (s ? "," : "") << "\n    {\"dir\": \"" << fs::path(dirs[s]).filename().string() << "\", \"distinct\": "
          << info.distinct << ", \"occurrences\": " << info.occurrences << ", \"genomes\": [";
        for (std::size_t g = 0; g < info.genomes.size(); ++g)
            m << (g ? ", " : "") << "{\"name\": \"" << json_escape(info.genomes[g].name) << "\", \"path\": \""
              << json_escape(info.genomes[g].path) << "\"}";
        m << "]}";
    }
    m << "\n  ]\n}\n";
    if (!m) throw std::runtime_error("cannot write " + path);
}

std::string next_segment_dir(const std::string& store) {
    long next = 0;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(store, ec)) {
        const std::string name = e.path().filename().string();
        if (name.rfind("seg-", 0) == 0 && name.size() > 4 && name.find_first_not_of("0123456789", 4) == std::string::npos)
            next = std::max(next, std::stol(name.substr(4)) + 1);
    }
    return (fs::path(store) / ("seg-" + std::to_string(next))).string();
}

std::string human_bytes(std::uint64_t b) {
    char buf[32];
    if (b >= (1ULL << 30)) std::snprintf(buf, sizeof buf, "%.2f GB", static_cast<double>(b) / (1ULL << 30));
    else std::snprintf(buf, sizeof buf, "%.1f MB", static_cast<double>(b) / (1ULL << 20));
    return buf;
}

} // namespace

int run_background_command(const std::vector<std::string>& args) {
    if (args.empty()) { print_background_help(); return 0; }
    std::string out, dump, tmp_dir;
    std::vector<bg::Genome> genomes;
    std::size_t threads = 0, memory = std::size_t{2} << 30, partitions = 4096;
    bool quiet = false, add = false, merge = false, partitions_given = false;
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
        } else if (a == "--add") add = true;
        else if (a == "--merge") merge = true;
        else if (a == "--dump") dump = value();
        else if (a == "--partitions") { partitions = cli::parse_size_arg(a, value()); partitions_given = true; }
        else if (a == "--memory") memory = parse_memory(value());
        else if (a == "--tmp-dir") tmp_dir = value();
        else if (a == "-t" || a == "--threads") threads = cli::parse_size_arg(a, value());
        else if (a == "-q" || a == "--quiet") quiet = true;
        else throw std::runtime_error("Unknown option: " + a);
    }
    if (out.empty()) throw std::runtime_error("background requires -o/--out <store.bgk>");
    if (static_cast<int>(add) + static_cast<int>(merge) + static_cast<int>(!dump.empty()) > 1)
        throw std::runtime_error("background: --add, --merge and --dump are separate runs");
    if (threads == 0) threads = std::max(1u, std::thread::hardware_concurrency());
    const bool exists = fs::exists(fs::path(out) / "manifest.json");
    cli::RunLog log("background", quiet);

    if (!dump.empty() || merge) {
        if (!genomes.empty()) throw std::runtime_error("background: --dump and --merge take no --genome");
        if (!exists) throw std::runtime_error(out + " is not a background store (no manifest.json)");
    }
    if (!dump.empty()) {
        const bg::Store store(out);
        std::FILE* f = std::fopen(dump.c_str(), "w");
        if (f == nullptr) throw std::runtime_error("cannot write " + dump);
        std::uint64_t n = 0;
        for (std::uint32_t p = 0; p < store.partitions(); ++p)
            store.for_each(p, [&](std::uint64_t code, std::uint64_t count) {
                std::fprintf(f, "%s\t%llu\n", kmer::decode(code, bg::kK).c_str(), static_cast<unsigned long long>(count));
                ++n;
            });
        if (std::fclose(f) != 0) throw std::runtime_error("cannot write " + dump);
        log.info("dumped " + std::to_string(n) + " k-mers from " + std::to_string(store.segments()) + " segment(s)");
        log.wrote({dump});
        log.done();
        return 0;
    }

    if (merge) {
        const auto dirs = bg::segment_dirs_of(out);
        if (dirs.empty()) throw std::runtime_error(out + " holds no segment");
        if (dirs.size() == 1) { log.info("one segment: nothing to merge"); log.done(); return 0; }
        SegmentInfo merged;
        for (const auto& d : dirs) {
            const SegmentInfo info = read_segment_info(d);
            merged.genomes.insert(merged.genomes.end(), info.genomes.begin(), info.genomes.end());
            merged.occurrences += info.occurrences;
        }
        const std::string target = next_segment_dir(out);
        bg::MergeOptions mo;
        mo.segment_dirs = dirs;
        mo.out_segment_dir = target;
        mo.threads = threads;
        log.info("merging " + std::to_string(dirs.size()) + " segments into " + fs::path(target).filename().string());
        bg::SegmentStats st;
        try {
            st = bg::merge_segments(mo);
        } catch (...) {
            std::error_code ec;
            fs::remove_all(target, ec);
            throw;
        }
        merged.distinct = st.distinct;
        write_segment_info(target, merged);
        // The merged segment is complete: only now do the old ones go.
        for (const auto& d : dirs) fs::remove_all(d);
        const bg::Store store(out);
        write_manifest(out, store.partitions());
        log.info("merged: " + std::to_string(st.distinct) + " distinct 31-mers, " + human_bytes(st.bytes));
        log.wrote({(fs::path(out) / "manifest.json").string(), target});
        log.done();
        return 0;
    }

    if (genomes.empty()) throw std::runtime_error("background: give --genome NAME=PATH (or --genomes)");
    if (exists && !add) throw std::runtime_error(out + " already holds a store: --add adds a segment to it");
    if (!exists && add) throw std::runtime_error(out + " is not a background store (no manifest.json) to add to");
    std::set<std::string> names;
    if (exists)
        for (const auto& d : bg::segment_dirs_of(out))
            for (const auto& g : read_segment_info(d).genomes) names.insert(g.name);
    for (auto& g : genomes) {
        if (!names.insert(g.name).second) throw std::runtime_error("background: genome name " + g.name + " is already used");
        g.path = fs::absolute(g.path).string();
    }
    if (exists) {
        const bg::Store store(out);
        if (partitions_given && partitions != store.partitions())
            throw std::runtime_error("background: the store has " + std::to_string(store.partitions()) + " partitions");
        partitions = store.partitions();
    }
    fs::create_directories(out);
    bg::BuildOptions o;
    o.genomes = genomes;
    o.segment_dir = next_segment_dir(out);
    o.tmp_dir = tmp_dir.empty() ? out : tmp_dir;
    o.partitions = static_cast<std::uint32_t>(partitions);
    o.threads = threads;
    o.memory_bytes = memory;
    o.log = [&](const std::string& m) { log.info(m); };
    bg::SegmentStats st;
    try {
        st = bg::build_segment(o);
    } catch (...) {
        std::error_code ec;
        fs::remove_all(o.segment_dir, ec);
        throw;
    }
    write_segment_info(o.segment_dir, SegmentInfo{genomes, st.distinct, st.kmers});
    write_manifest(out, o.partitions);
    log.info(fs::path(o.segment_dir).filename().string() + ": " + std::to_string(st.distinct) + " distinct 31-mers of " +
             std::to_string(st.kmers) + " occurrences, " + human_bytes(st.bytes) + "; spill " + human_bytes(st.spill_bytes) +
             " in " + std::to_string(st.spill_files) + " files");
    log.wrote({(fs::path(out) / "manifest.json").string(), o.segment_dir});
    log.done();
    return 0;
}

} // namespace panvar
