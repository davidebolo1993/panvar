#include "panvar/background_command.hpp"

#include "panvar/anchors.hpp"
#include "panvar/background.hpp"
#include "panvar/background_graph.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gz_reader.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <unistd.h>

namespace panvar {

namespace {

namespace fs = std::filesystem;

void print_background_help() {
    std::cout
        << "Usage:\n"
        << "  panvar background -o <store.bgk> [--add] --genome NAME=PATH [...] [--coordinate-reference NAME]\n"
        << "  panvar background -o <store.bgk> --segment-only <dir> --genome NAME=PATH [...]   (then --attach <dir>)\n"
        << "  panvar background -o <store.bgk> --merge [--range I/N] | --compact | --clean\n\n"
        << "Count every 31-mer of a set of genomes once, with anchors that find any interval of the coordinate\n"
        << "reference in every genome, for recruit.\n\n"
        << "Building:\n"
        << "  -o, --out <dir>                 The store (a folder; manifest.json lists its segments)\n"
        << "      --genome NAME=PATH          A genome: FASTA (plain, gzip or bgzip) or a GFA of haplotypes; repeatable\n"
        << "      --genomes <tsv>             Genomes as lines of name and path\n"
        << "      --coordinate-reference NAME The genome whose coordinates the locus graphs carry: its anchors are\n"
        << "                                  placed in every genome (on an existing store: adds them)\n"
        << "      --add                       Count the genomes into a new segment of an existing store\n"
        << "      --segment-only <dir>        Count the genomes into <dir> without committing (a parallel job)\n"
        << "      --attach <dir>...           Commit segments built with --segment-only\n"
        << "      --partitions <N>            Partitions of a new store, a power of two (default 4096)\n"
        << "      --memory <size>             Memory bound for counting, e.g. 2G, 512M (default 2G)\n"
        << "      --tmp-dir <dir>             Temporary files (default: $TMPDIR, else inside the store)\n\n"
        << "Finishing:\n"
        << "      --merge                     Merge the store's segments into one\n"
        << "      --range I/N                 With --merge: only the I-th of N partition ranges (0-based), one job each\n"
        << "      --compact                   Write the merged store as an SSHash dictionary, <store>.sshash\n"
        << "      --clean                     Remove what no committed segment or running job uses\n\n"
        << "Inspecting:\n"
        << "      --locate <region>           Where each genome holds a region of the coordinate reference\n"
        << "      --dump <path>               Every k-mer and its count (segments summed), as TSV\n"
        << "      --verify <path>             Check that another background (store or .sshash) holds the same counts\n\n"
        << "  -t, --threads <N>               Worker threads (default: Slurm's allocation or the CPUs available)\n"
        << "  -q, --quiet                     No progress logs\n"
        << "  -h, --help                      Show this help\n";
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

// Genome and reference names become file names: letters, digits, '.', '_', '-', not starting with '.'.
void check_name(const std::string& what, const std::string& n) {
    if (n.empty() || n.size() > 100 || n[0] == '.' ||
        n.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-") != std::string::npos)
        throw std::runtime_error(what + " '" + n + "': use letters, digits, '.', '_' and '-' (not starting with '.'), at most 100");
}

// Test hooks: PANVAR_TEST_FAIL_AT=<point> ends the process abruptly there, as a kill would;
// PANVAR_TEST_WAIT_AT=<point>:<file> waits there (up to two minutes) until <file> exists, so a test can
// slot another job in at that moment.
void fail_point(const char* point) {
    const char* at = std::getenv("PANVAR_TEST_FAIL_AT");
    if (at != nullptr && std::string(at) == point) {
        std::fflush(nullptr);
        std::_Exit(137);
    }
    const char* wait = std::getenv("PANVAR_TEST_WAIT_AT");
    if (wait != nullptr) {
        const std::string w(wait);
        const auto colon = w.find(':');
        if (colon != std::string::npos && w.substr(0, colon) == point)
            for (int i = 0; i < 1200 && !fs::exists(w.substr(colon + 1)); ++i) std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}

std::string human_bytes(std::uint64_t b) {
    char buf[32];
    if (b >= (1ULL << 30)) std::snprintf(buf, sizeof buf, "%.2f GiB", static_cast<double>(b) / (1ULL << 30));
    else std::snprintf(buf, sizeof buf, "%.1f MiB", static_cast<double>(b) / (1ULL << 20));
    return buf;
}

std::string random_id() {
    std::random_device rd;
    const std::uint64_t x = kmer::mix64((static_cast<std::uint64_t>(rd()) << 32) ^ rd() ^
                                        static_cast<std::uint64_t>(std::chrono::steady_clock::now().time_since_epoch().count()) ^
                                        static_cast<std::uint64_t>(::getpid()));
    char buf[17];
    std::snprintf(buf, sizeof buf, "%016llx", static_cast<unsigned long long>(x));
    return buf;
}

std::string default_tmp(const std::string& given, const std::string& fallback) {
    if (!given.empty()) return given;
    if (const char* t = std::getenv("TMPDIR")) if (*t != '\0' && fs::is_directory(t)) return t;
    return fallback;
}

bool is_gfa(const std::string& path) {
    GzLineReader in(path);
    std::string line;
    while (in.ok() && in.getline(line)) {
        if (line.empty()) continue;
        return line.size() > 1 && line[1] == '\t' && std::string("HSLPWCJ#").find(line[0]) != std::string::npos;
    }
    return false;
}

fs::path ref_path(const std::string& store, const std::string& ref) { return fs::path(store) / "anchors" / ref / "reference.anc"; }

// The store's coordinate references, loaded.
std::vector<std::unique_ptr<anchors::ReferenceAnchors>> load_references(const std::string& store, const bg::Manifest& m) {
    std::vector<std::unique_ptr<anchors::ReferenceAnchors>> out;
    for (const auto& [name, fp] : m.references) {
        out.push_back(std::make_unique<anchors::ReferenceAnchors>(anchors::load_reference(ref_path(store, name).string())));
        if (out.back()->fingerprint() != fp) throw std::runtime_error("background: the anchor set of " + name + " does not match the manifest");
    }
    return out;
}

// Maps the genomes of a segment (FASTA or graph) on a reference: the segment's anchor file for it.
void map_segment(const bg::SegmentInfo& info, const anchors::ReferenceAnchors& ref, const std::string& out_path, std::size_t threads,
                 cli::RunLog& log) {
    const anchors::AnchorLookup lookup(ref);
    std::vector<anchors::Hit> hits;
    std::vector<std::string> names;
    std::vector<std::vector<std::string>> contigs(info.genomes.size());
    std::vector<std::vector<std::uint64_t>> lengths(info.genomes.size());
    for (std::size_t g = 0; g < info.genomes.size(); ++g) {
        const auto& ge = info.genomes[g];
        names.push_back(ge.name);
        std::error_code ec;
        if (fs::file_size(ge.path, ec) != ge.bytes || ec)
            throw std::runtime_error("background: " + ge.path + " (" + ge.name + ") is not the file that was counted (size differs)");
        if (is_gfa(ge.path)) {
            auto h = bg::graph::graph_anchors(ge.path, lookup, static_cast<std::uint32_t>(g), threads);
            hits.insert(hits.end(), h.begin(), h.end());
            contigs[g] = bg::graph::graph_contigs(ge.path, &lengths[g]);   // as counting saw them
        } else {
            auto h = anchors::map_fasta(lookup, ge.path, static_cast<std::uint32_t>(g), threads, contigs[g], lengths[g]);
            hits.insert(hits.end(), h.begin(), h.end());
        }
        log.info("anchors: " + ref.genome + "; mapped in: " + ge.name);
    }
    anchors::write_segment_anchors(out_path, ref, names, contigs, lengths, hits);
}

// The next seg-<n> name that the manifest has not used and no leftover folder holds.
std::string take_segment_name(const std::string& store, bg::Manifest& m) {
    while (fs::exists(fs::path(store) / ("seg-" + std::to_string(m.next_segment)))) ++m.next_segment;
    return "seg-" + std::to_string(m.next_segment++);
}

// Records a coordinate reference's anchor file in a committed segment (replacing an earlier record of it).
void add_anchor_record(const fs::path& seg_dir, const std::string& ref, std::uint64_t fp) {
    bg::SegmentInfo info = bg::read_segment_info(seg_dir.string());
    info.anchor_sets.erase(std::remove_if(info.anchor_sets.begin(), info.anchor_sets.end(), [&](const auto& x) { return x.first == ref; }),
                           info.anchor_sets.end());
    info.anchor_sets.emplace_back(ref, fp);
    bg::write_segment_info(seg_dir.string(), info);
}

// A copy of a segment folder at dst: files hard-linked where the filesystem allows, else copied.
void link_or_copy_dir(const fs::path& src, const fs::path& dst) {
    fs::create_directories(dst);
    for (const auto& e : fs::directory_iterator(src)) {
        if (!e.is_regular_file() || e.path().filename() == "job.lock") continue;
        std::error_code ec;
        fs::create_hard_link(e.path(), dst / e.path().filename(), ec);
        if (ec) fs::copy_file(e.path(), dst / e.path().filename(), fs::copy_options::overwrite_existing);
    }
}

// Maps genomes of committed segments on a reference into anchor files under `dir` (one per segment, named
// after it), skipping segments already mapped there.
void map_segments(const std::string& out, const std::vector<bg::ManifestSegment>& segs, const anchors::ReferenceAnchors& ref,
                  const fs::path& dir, const std::string& suffix, std::size_t threads, cli::RunLog& log) {
    for (const auto& s : segs) {
        const fs::path f = dir / (s.dir + suffix);
        if (fs::exists(f)) continue;
        map_segment(bg::read_segment_info((fs::path(out) / s.dir).string()), ref, f.string(), threads, log);
    }
}

// Genome identity: a name, or a sequence, already in the store.
void check_new_genomes(const bg::Manifest& m, const std::vector<bg::Genome>& genomes) {
    std::map<std::string, std::string> by_checksum;
    std::set<std::string> names;
    for (const auto& s : m.segments)
        for (const auto& g : s.genomes) { names.insert(g.name); if (!g.checksum.empty()) by_checksum[g.checksum] = g.name; }
    for (const auto& g : genomes) {
        if (names.count(g.name)) throw std::runtime_error("background: genome " + g.name + " is already in the store");
        if (!g.checksum.empty()) {
            const auto it = by_checksum.find(g.checksum);
            if (it != by_checksum.end())
                throw std::runtime_error("background: genome " + g.name + " has the same sequence as " + it->second + ", already counted");
            by_checksum[g.checksum] = g.name;
        }
        names.insert(g.name);
    }
}

struct Staging {
    fs::path dir;
    std::unique_ptr<bg::JobLock> lock;
    ~Staging() { lock.reset(); std::error_code ec; if (!dir.empty()) fs::remove_all(dir, ec); }
};

}  // namespace

int run_background_command(const std::vector<std::string>& args) {
    if (args.empty()) { print_background_help(); return 0; }
    std::string out, dump, tmp_dir, segment_only, verify, coordinate_reference, locate_region, range, spell;
    std::vector<std::string> attach;
    std::vector<bg::Genome> genomes;
    std::size_t threads = 0, memory = std::size_t{2} << 30, partitions = 4096, minimizer = 0;
    bool quiet = false, add = false, merge = false, compact = false, partitions_given = false, clean = false, local = false, dedup = true;
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
            if (eq == std::string::npos || eq == 0 || eq + 1 == v.size()) throw std::runtime_error("--genome takes NAME=PATH, got " + v);
            genomes.push_back({v.substr(0, eq), v.substr(eq + 1), 0, ""});
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
                genomes.push_back({line.substr(0, tab), cli::trim_ascii(line.substr(tab + 1)), 0, ""});
            }
        } else if (a == "--add") add = true;
        else if (a == "--segment-only") segment_only = value();
        else if (a == "--attach") {
            attach.push_back(value());
            while (i + 1 < args.size() && args[i + 1].rfind("-", 0) != 0) attach.push_back(args[++i]);
        } else if (a == "--merge") merge = true;
        else if (a == "--range") range = value();
        else if (a == "--compact") compact = true;
        else if (a == "--local-strings") local = true;
        else if (a == "--minimizer") minimizer = cli::parse_size_arg(a, value());
        else if (a == "--verify") verify = value();
        else if (a == "--coordinate-reference") coordinate_reference = value();
        else if (a == "--locate") locate_region = value();
        else if (a == "--dump") dump = value();
        else if (a == "--clean") clean = true;
        else if (a == "--spell-graph") spell = value();
        else if (a == "--partitions") { partitions = cli::parse_size_arg(a, value()); partitions_given = true; }
        else if (a == "--memory") memory = parse_memory(value());
        else if (a == "--no-dedup") dedup = false;
        else if (a == "--tmp-dir") tmp_dir = value();
        else if (a == "-t" || a == "--threads") threads = cli::parse_size_arg(a, value());
        else if (a == "-q" || a == "--quiet") quiet = true;
        else throw std::runtime_error("Unknown option: " + a);
    }
    cli::RunLog log("background", quiet);
    if (threads == 0) threads = bg::default_threads();

    if (!spell.empty()) {
        const auto eq = spell.find('=');
        if (eq == std::string::npos) throw std::runtime_error("--spell-graph takes GFA=FASTA");
        bg::graph::spell_graph(spell.substr(0, eq), spell.substr(eq + 1));
        log.wrote({spell.substr(eq + 1)});
        log.done();
        return 0;
    }
    if (out.empty()) throw std::runtime_error("background requires -o/--out <store.bgk>");
    const int modes = static_cast<int>(add) + static_cast<int>(merge) + static_cast<int>(compact) + static_cast<int>(!dump.empty()) +
        static_cast<int>(!verify.empty()) + static_cast<int>(!locate_region.empty()) + static_cast<int>(clean) +
        static_cast<int>(!segment_only.empty()) + static_cast<int>(!attach.empty());
    if (modes > 1) throw std::runtime_error("background: --add, --segment-only, --attach, --merge, --compact, --dump, --verify, --locate and --clean are separate runs");
    if (!range.empty() && !merge) throw std::runtime_error("--range goes with --merge");
    if (local && !compact) throw std::runtime_error("--local-strings goes with --compact");
    if (partitions_given && (partitions == 0 || partitions > bg::kMaxPartitions || (partitions & (partitions - 1)) != 0))
        throw std::runtime_error("--partitions must be a power of two, at most 65536");
    for (auto& g : genomes) {
        check_name("genome name", g.name);
        if (!fs::is_regular_file(g.path)) throw std::runtime_error("background: genome " + g.name + ": no file at " + g.path);
        g.path = fs::absolute(g.path).string();
    }
    if (!coordinate_reference.empty()) check_name("coordinate reference", coordinate_reference);
    const bool exists = bg::store_exists(out);
    const bool reads_store = !dump.empty() || merge || compact || !verify.empty() || !locate_region.empty() || clean || !attach.empty();
    if (reads_store && !genomes.empty()) throw std::runtime_error("background: this run takes no --genome");
    if ((reads_store || add || !segment_only.empty()) && !exists) throw std::runtime_error(out + " is not a background store (no manifest.json)");

    // ---- read-only runs
    if (!dump.empty()) {
        const bg::Store store(out);
        std::FILE* f = std::fopen(dump.c_str(), "w");
        if (f == nullptr) throw std::runtime_error("cannot write " + dump);
        std::uint64_t n = 0;
        bool ok = true;
        for (std::uint32_t p = 0; p < store.partitions(); ++p)
            store.for_each(p, [&](std::uint64_t code, std::uint64_t count) {
                if (std::fprintf(f, "%s\t%llu\n", kmer::decode(code, bg::kK).c_str(), static_cast<unsigned long long>(count)) < 0) ok = false;
                ++n;
            });
        if (std::fclose(f) != 0 || !ok) throw std::runtime_error("cannot write " + dump);
        log.info("dumped: k-mers: " + std::to_string(n) + "; segments: " + std::to_string(store.segments()));
        log.wrote({dump});
        log.done();
        return 0;
    }

    if (!locate_region.empty()) {
        const auto colon = locate_region.rfind(':'), dash = locate_region.rfind('-');
        if (colon == std::string::npos || dash == std::string::npos || dash < colon)
            throw std::runtime_error("--locate takes CONTIG:START-END, got " + locate_region);
        const std::string contig = locate_region.substr(0, colon);
        const std::uint64_t start = cli::parse_size_arg("--locate", locate_region.substr(colon + 1, dash - colon - 1));
        const std::uint64_t end = cli::parse_size_arg("--locate", locate_region.substr(dash + 1));
        if (start >= end) throw std::runtime_error("--locate: the start must be below the end");
        const bg::Manifest m = bg::read_manifest(out);
        std::string ref_name = coordinate_reference;
        if (ref_name.empty()) {
            if (m.references.size() != 1)
                throw std::runtime_error("--locate: give --coordinate-reference (the store has " + std::to_string(m.references.size()) + " coordinate references)");
            ref_name = m.references[0].first;
        }
        const auto rit = std::find_if(m.references.begin(), m.references.end(), [&](const auto& r) { return r.first == ref_name; });
        if (rit == m.references.end()) throw std::runtime_error("--locate: " + ref_name + " is not a coordinate reference of the store");
        const anchors::ReferenceAnchors ref = anchors::load_reference(ref_path(out, ref_name).string());
        std::vector<std::unique_ptr<anchors::SegmentAnchors>> segs;
        std::vector<const anchors::SegmentAnchors*> raw;
        for (const auto& s : m.segments) {
            segs.push_back(std::make_unique<anchors::SegmentAnchors>((fs::path(out) / s.dir / ("anchors-" + ref_name + ".bin")).string()));
            raw.push_back(segs.back().get());
        }
        std::uint64_t bytes = 0;
        const auto spans = anchors::locate(ref, raw, contig, start, end, {}, &bytes);
        std::cout << "genome\tcopy\tcontig\tstart\tend\tstrand\tanchors_inside\tanchors_flank\tchains\n";
        for (const auto& s : m.segments)
            for (const auto& g : s.genomes) {
                std::size_t copy = 0;
                for (const auto& sp : spans) {
                    if (sp.genome != g.name) continue;
                    std::cout << g.name << '\t' << ++copy << '\t' << sp.contig << '\t' << sp.start << '\t' << sp.end << '\t'
                              << (sp.forward ? '+' : '-') << '\t' << sp.anchors_inside << '\t' << sp.anchors_flank << '\t' << sp.chains << '\n';
                }
                if (copy == 0) std::cout << g.name << "\t0\t.\t.\t.\t.\t0\t0\t0\n";
            }
        log.info("anchor places: bytes read: " + std::to_string(bytes));
        log.done();
        return 0;
    }

    if (!verify.empty()) {
        const bg::Store store(out);
        const auto other = bg::Background::open(verify);
        log.info("verify: " + verify + " (" + other->kind() + "); against: " + out);
        std::uint64_t checked = 0, wrong = 0;
        std::vector<std::uint64_t> codes, counts;
        for (std::uint32_t p = 0; p < store.partitions(); ++p) {
            codes.clear(); counts.clear();
            store.for_each(p, [&](std::uint64_t code, std::uint64_t count) { codes.push_back(code); counts.push_back(count); });
            const auto got = other->count_batch(codes);
            for (std::size_t i = 0; i < codes.size(); ++i) {
                ++checked;
                if (got[i] != counts[i] && wrong++ < 5)
                    log.info("differs: " + kmer::decode(codes[i], bg::kK) + "; here: " + std::to_string(counts[i]) + "; there: " + std::to_string(got[i]));
            }
        }
        // random k-mers absent from this store must be absent there too
        std::uint64_t x = 0x243f6a8885a308d3ULL, absent = 0, absent_wrong = 0;
        std::vector<std::uint64_t> probe;
        for (int t = 0; t < 1000000; ++t) {
            x = kmer::mix64(x + static_cast<std::uint64_t>(t));
            probe.push_back(kmer::canonical(x & kmer::mask_of(bg::kK), bg::kK));
        }
        const auto here = store.count_batch(probe), there = other->count_batch(probe);
        for (std::size_t i = 0; i < probe.size(); ++i) if (here[i] == 0) { ++absent; absent_wrong += there[i] != 0; }
        log.info("random absent k-mers: " + std::to_string(absent) + "; found there: " + std::to_string(absent_wrong));
        wrong += absent_wrong;
        const std::uint64_t other_distinct = other->distinct();
        log.info("checked: k-mers: " + std::to_string(checked) + "; differ: " + std::to_string(wrong) + "; distinct there: " +
                 std::to_string(other_distinct) + "; distinct here: " + std::to_string(checked));
        log.done();
        return wrong == 0 && other_distinct == checked ? 0 : 1;
    }

    if (compact) {
        bg::CompactOptions co;
        co.store_dir = out;
        std::string base = out;
        while (base.size() > 1 && base.back() == '/') base.pop_back();
        co.out_path = base + ".sshash";
        const std::string tmp_root = default_tmp(tmp_dir, out);
        co.tmp_dir = (fs::path(tmp_root) / ("panvar-compact-" + random_id())).string();
        co.threads = threads;
        co.minimizer = minimizer;
        co.ram_gb = std::max<std::size_t>(1, memory >> 30);
        co.local_strings = local;
        co.quiet = quiet;
        co.log = [&](const std::string& m) { log.info(m); };
        bg::CompactStats st;
        {
            const bg::JobLock one_at_a_time((fs::path(out) / ".compact").string());   // throws when another --compact runs
            fail_point("compact-locked");
            Staging tmp{co.tmp_dir, std::make_unique<bg::JobLock>(co.tmp_dir)};
            st = bg::compact_store(co);
        }
        {
            char b[256];
            std::snprintf(b, sizeof b, "compacted: k-mers: %llu; strings: %llu; size: %s; bits per k-mer: %.2f; seconds: strings %.1f, SSHash %.1f",
                          static_cast<unsigned long long>(st.kmers), static_cast<unsigned long long>(st.strings), human_bytes(st.bytes).c_str(),
                          static_cast<double>(st.bytes) * 8 / static_cast<double>(std::max<std::uint64_t>(1, st.kmers)), st.strings_seconds, st.sshash_seconds);
            log.info(b);
        }
        log.wrote({co.out_path});
        log.done();
        return 0;
    }

    if (clean) {
        const bg::StoreLock lock(out);
        const bg::Manifest m = bg::read_manifest(out);
        std::set<std::string> listed;
        for (const auto& s : m.segments) listed.insert(s.dir);
        std::uint64_t removed = 0;
        auto remove = [&](const fs::path& p) { std::error_code ec; removed += fs::remove_all(p, ec) > 0; log.info("removed: " + p.filename().string()); };
        std::error_code ec;
        for (const auto& e : fs::directory_iterator(out, ec)) {
            const std::string name = e.path().filename().string();
            if (e.is_directory() && name.rfind("seg-", 0) == 0 && !listed.count(name)) remove(e.path());
            else if (e.is_directory() && (name.rfind("panvar-background-", 0) == 0 || name.rfind("panvar-compact-", 0) == 0) && !bg::job_alive(e.path().string())) remove(e.path());
            else if (e.is_directory() && name.rfind(".merge-", 0) == 0) {
                // a merge stays while every segment it plans is listed, or while one of its range jobs runs
                bool live = false, plan_ok = true;
                std::ifstream plan(e.path() / "plan.txt");
                std::string d;
                std::size_t planned = 0;
                while (std::getline(plan, d)) if (!d.empty()) { ++planned; plan_ok = plan_ok && listed.count(d); }
                std::error_code ec2;
                for (const auto& r : fs::directory_iterator(e.path(), ec2)) live = live || (r.is_directory() && bg::job_alive(r.path().string()));
                if (!live && (planned == 0 || !plan_ok)) remove(e.path());
            }
            else if (!e.is_directory() && name.size() > 4 && name.compare(name.size() - 4, 4, ".tmp") == 0) remove(e.path());
        }
        for (const auto& e : fs::directory_iterator(fs::path(out) / ".staging", ec))
            if (!bg::job_alive(e.path().string())) remove(e.path());
        log.info("clean: items removed: " + std::to_string(removed));
        log.done();
        return 0;
    }

    if (merge) {
        bg::Manifest m = bg::read_manifest(out);
        if (m.segments.size() <= 1) { log.info("merge: one segment or none, nothing to merge"); log.done(); return 0; }
        std::uint32_t ranges = 1, mine = 0;
        bool all_ranges = true;
        if (!range.empty()) {
            const auto slash = range.find('/');
            if (slash == std::string::npos) throw std::runtime_error("--range takes I/N");
            mine = static_cast<std::uint32_t>(cli::parse_size_arg("--range", range.substr(0, slash)));
            const std::size_t n = cli::parse_size_arg("--range", range.substr(slash + 1));
            if (n == 0 || n > m.partitions || mine >= n) throw std::runtime_error("--range: I must be below N, N at most the partition count");
            ranges = static_cast<std::uint32_t>(n);
            all_ranges = false;
        }
        // A merge of N ranges already under way, whose planned segments are all still listed, is joined: its
        // plan stands even if segments were added since. Otherwise a new merge plans the segments listed now.
        fs::path mdir;
        {
            std::map<std::string, bg::ManifestSegment> listed;
            for (const auto& s : m.segments) listed[s.dir] = s;
            std::error_code ec;
            for (const auto& e : fs::directory_iterator(out, ec)) {
                const std::string name = e.path().filename().string();
                const std::string suffix = "-" + std::to_string(ranges);
                if (!e.is_directory() || name.rfind(".merge-", 0) != 0 || name.size() <= suffix.size() ||
                    name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
                    continue;
                std::ifstream in(e.path() / "plan.txt");
                std::vector<bg::ManifestSegment> planned;
                bool ok = static_cast<bool>(in);
                for (std::string d; ok && std::getline(in, d);) {
                    if (d.empty()) continue;
                    const auto it = listed.find(d);
                    if (it == listed.end()) ok = false;
                    else planned.push_back(it->second);
                }
                if (ok && planned.size() > 1) { mdir = e.path(); m.segments = planned; break; }
            }
        }
        if (mdir.empty()) mdir = fs::path(out) / (".merge-" + std::to_string(m.generation) + "-" + std::to_string(ranges));
        fs::create_directories(mdir);
        std::vector<std::string> seg_dirs;
        std::string plan;
        for (const auto& s : m.segments) { seg_dirs.push_back((fs::path(out) / s.dir).string()); plan += s.dir + "\n"; }
        {
            const fs::path pp = mdir / "plan.txt";
            if (fs::exists(pp)) {
                std::ifstream in(pp);
                std::stringstream ss;
                ss << in.rdbuf();
                if (ss.str() != plan) throw std::runtime_error("background: " + mdir.string() + " plans another merge");
            } else {
                const std::string tmp = (mdir / ("plan.txt." + random_id())).string();
                { std::ofstream o(tmp); o << plan; }
                fs::rename(tmp, pp);
            }
        }
        auto part_name = [&](std::uint32_t i) { return "parts-" + std::to_string(i) + ".bin"; };
        for (std::uint32_t i = 0; i < ranges; ++i) {
            if (!all_ranges && i != mine) continue;
            const fs::path pf = mdir / part_name(i);
            if (fs::exists(pf)) { log.info("merge: range: " + std::to_string(i) + "; already merged"); continue; }
            const bg::JobLock rl((mdir / ("range-" + std::to_string(i))).string());
            bg::MergeRangeOptions mo;
            mo.segment_dirs = seg_dirs;
            mo.out_path = pf.string();
            mo.first = static_cast<std::uint32_t>(static_cast<std::uint64_t>(m.partitions) * i / ranges);
            mo.end = static_cast<std::uint32_t>(static_cast<std::uint64_t>(m.partitions) * (i + 1) / ranges);
            mo.threads = threads;
            log.info("merge: partitions: " + std::to_string(mo.first) + "-" + std::to_string(mo.end) + "; segments: " + std::to_string(m.segments.size()));
            mo.progress = !quiet;
            bg::merge_range(mo);
            fail_point("merge-range");
        }
        for (std::uint32_t i = 0; i < ranges; ++i)
            if (!fs::exists(mdir / part_name(i))) {
                log.info("merge: range: " + std::to_string(mine) + " of " + std::to_string(ranges) + "; done (the last range to finish commits)");
                log.done();
                return 0;
            }
        // every range is done: anchors, the segment record, and the commit (one job does it)
        std::unique_ptr<bg::JobLock> committing;
        try {
            committing = std::make_unique<bg::JobLock>((mdir / "commit").string());
        } catch (const std::exception&) {
            log.info("merge: another job is committing");
            log.done();
            return 0;
        }
        bg::SegmentInfo info;
        info.partitions = m.partitions;
        for (std::uint32_t i = 0; i < ranges; ++i) {
            info.parts.push_back({part_name(i), static_cast<std::uint32_t>(static_cast<std::uint64_t>(m.partitions) * i / ranges),
                                  static_cast<std::uint32_t>(static_cast<std::uint64_t>(m.partitions) * (i + 1) / ranges)});
            info.distinct += bg::parts_distinct((mdir / part_name(i)).string());
        }
        for (const auto& s : m.segments) {
            info.genomes.insert(info.genomes.end(), s.genomes.begin(), s.genomes.end());
            info.occurrences += s.occurrences;
        }
        const auto refs = load_references(out, m);
        for (std::size_t r = 0; r < refs.size(); ++r) {
            const std::string name = m.references[r].first;
            std::vector<std::string> inputs;
            for (const auto& s : m.segments) inputs.push_back((fs::path(out) / s.dir / ("anchors-" + name + ".bin")).string());
            anchors::merge_segment_anchors(inputs, *refs[r], (mdir / ("anchors-" + name + ".bin")).string());
            info.anchor_sets.emplace_back(name, m.references[r].second);
        }
        bg::write_segment_info(mdir.string(), info);
        std::vector<std::string> old;
        {
            const bg::StoreLock lock(out);
            bg::Manifest now = bg::read_manifest(out);
            if (!fs::exists(mdir)) { log.info("merge: another job committed it"); log.done(); return 0; }
            // the merge replaces the segments it planned; segments committed since stay as they are
            std::set<std::string> planned, listed;
            for (const auto& s : m.segments) planned.insert(s.dir);
            for (const auto& s : now.segments) listed.insert(s.dir);
            for (const auto& d : planned)
                if (!listed.count(d)) throw std::runtime_error("background: " + d + " left the store during the merge; run --merge again (and --clean)");
            if (now.references != m.references) throw std::runtime_error("background: the coordinate references changed during the merge; run --merge again");
            std::error_code ec;
            for (std::uint32_t i = 0; i < ranges; ++i) fs::remove_all(mdir / ("range-" + std::to_string(i)), ec);
            fs::remove(mdir / "plan.txt", ec);
            const std::string dir = take_segment_name(out, now);
            fs::rename(mdir, fs::path(out) / dir);
            fs::remove_all(fs::path(out) / dir / "commit", ec);
            bg::ManifestSegment ms;
            ms.dir = dir;
            ms.distinct = info.distinct;
            ms.occurrences = info.occurrences;
            ms.genomes = info.genomes;
            std::vector<bg::ManifestSegment> kept = {ms};
            for (const auto& s : now.segments) {
                if (planned.count(s.dir)) old.push_back(s.dir);
                else kept.push_back(s);
            }
            now.segments = kept;
            ++now.generation;
            fail_point("merge-before-manifest");
            bg::write_manifest(out, now);
            log.info("merged: segment: " + dir + "; distinct 31-mers: " + std::to_string(info.distinct));
        }
        fail_point("merge-after-commit");
        // the merged-away segments stay on disk: jobs still reading them would fail; --clean removes them
        log.info("merge: merged-away segments: " + std::to_string(old.size()) + "; kept until --clean (run it when no job reads the store)");
        log.wrote({(fs::path(out) / "manifest.json").string()});
        log.done();
        return 0;
    }

    if (!attach.empty()) {
        // everything is checked first, then each segment is linked (or copied) into the store, the manifest
        // committed, and only then are the sources removed: a failure or a kill leaves them where they were
        std::vector<bg::SegmentInfo> infos;
        for (const auto& a : attach) {
            infos.push_back(bg::read_segment_info(a));
            const bg::Segment check(a);   // complete and readable
        }
        const bg::StoreLock lock(out);
        bg::Manifest m = bg::read_manifest(out);
        std::vector<bg::Genome> incoming;
        for (std::size_t i = 0; i < attach.size(); ++i) {
            const auto& info = infos[i];
            if (bg::job_alive(attach[i])) throw std::runtime_error("background: " + attach[i] + " is still being built");
            if (info.partitions != m.partitions) throw std::runtime_error("background: " + attach[i] + " has " + std::to_string(info.partitions) + " partitions, the store " + std::to_string(m.partitions));
            for (const auto& ref_entry : m.references) {
                const std::string& name = ref_entry.first;
                const std::uint64_t fp = ref_entry.second;
                const auto it = std::find_if(info.anchor_sets.begin(), info.anchor_sets.end(), [&](const auto& x) { return x.first == name; });
                if (it == info.anchor_sets.end() || it->second != fp)
                    throw std::runtime_error("background: " + attach[i] + " lacks the anchors of " + name + " (build it with --segment-only against this store)");
                anchors::SegmentAnchors((fs::path(attach[i]) / ("anchors-" + name + ".bin")).string());
            }
            incoming.insert(incoming.end(), info.genomes.begin(), info.genomes.end());
        }
        check_new_genomes(m, incoming);
        std::vector<std::string> made;
        try {
            for (std::size_t i = 0; i < attach.size(); ++i) {
                const std::string dir = take_segment_name(out, m);
                link_or_copy_dir(attach[i], fs::path(out) / dir);
                made.push_back(dir);
                m.segments.push_back({dir, infos[i].distinct, infos[i].occurrences, infos[i].genomes});
                log.info("attach: " + attach[i] + "; as: " + dir);
            }
            ++m.generation;
            fail_point("attach-before-manifest");
            bg::write_manifest(out, m);
        } catch (...) {
            for (const auto& d : made) { std::error_code ec; fs::remove_all(fs::path(out) / d, ec); }
            throw;
        }
        for (const auto& a : attach) { std::error_code ec; fs::remove_all(a, ec); }
        log.wrote({(fs::path(out) / "manifest.json").string()});
        log.done();
        return 0;
    }

    // ---- adding a coordinate reference to an existing store: every committed genome is mapped
    if (exists && genomes.empty() && !coordinate_reference.empty()) {
        const bg::Manifest m = bg::read_manifest(out);
        for (const auto& r : m.references) if (r.first == coordinate_reference) throw std::runtime_error("background: " + coordinate_reference + " is already a coordinate reference");
        const bg::Genome* rg = nullptr;
        for (const auto& s : m.segments) for (const auto& g : s.genomes) if (g.name == coordinate_reference) rg = &g;
        if (rg == nullptr) throw std::runtime_error("--coordinate-reference: " + coordinate_reference + " is not a genome of the store");
        if (is_gfa(rg->path)) throw std::runtime_error("--coordinate-reference: " + coordinate_reference + " is a graph; a coordinate reference must be a FASTA");
        Staging st{fs::path(out) / ".staging" / random_id(), nullptr};
        st.lock = std::make_unique<bg::JobLock>(st.dir.string());
        const anchors::ReferenceAnchors ref = anchors::build_reference(coordinate_reference, rg->path, threads, [&](const std::string& s) { log.info(s); });
        anchors::save(ref, (st.dir / "reference.anc").string());
        // segments committed while mapping are mapped in turn before the commit
        std::unique_ptr<bg::StoreLock> lock;
        bg::Manifest now = m;
        for (int round = 0;; ++round) {
            map_segments(out, now.segments, ref, st.dir, ".bin", threads, log);
            fail_point("reference-before-commit");
            lock = std::make_unique<bg::StoreLock>(out);
            now = bg::read_manifest(out);
            for (const auto& r : now.references) if (r.first == coordinate_reference) throw std::runtime_error("background: another job added " + coordinate_reference);
            bool all = true;
            for (const auto& s : now.segments) all = all && fs::exists(st.dir / (s.dir + ".bin"));
            if (all) break;
            if (round == 20) throw std::runtime_error("background: the store keeps changing; run again later");
            lock.reset();
        }
        for (const auto& s : now.segments) {
            const fs::path sd = fs::path(out) / s.dir;
            fs::rename(st.dir / (s.dir + ".bin"), sd / ("anchors-" + coordinate_reference + ".bin"));
            add_anchor_record(sd, coordinate_reference, ref.fingerprint());
        }
        fs::create_directories(ref_path(out, coordinate_reference).parent_path());
        fs::rename(st.dir / "reference.anc", ref_path(out, coordinate_reference));
        now.references.emplace_back(coordinate_reference, ref.fingerprint());
        ++now.generation;
        bg::write_manifest(out, now);
        log.info("coordinate reference: " + coordinate_reference + "; anchors: " + std::to_string(ref.code.size()) + "; mapped in every genome");
        log.done();
        return 0;
    }

    // ---- counting genomes: a new store, --add, or --segment-only
    if (genomes.empty()) throw std::runtime_error("background: give --genome NAME=PATH (or --genomes)");
    std::set<std::string> given;
    for (const auto& g : genomes) if (!given.insert(g.name).second) throw std::runtime_error("background: genome name " + g.name + " is given twice");
    if (!exists) {
        if (add) throw std::runtime_error(out + " is not a background store to add to");
        fs::create_directories(out);
        const bg::StoreLock lock(out);
        if (!bg::store_exists(out)) {
            bg::Manifest m;
            m.store_id = random_id();
            m.partitions = static_cast<std::uint32_t>(partitions);
            bg::write_manifest(out, m);
        }
    }
    bg::Manifest m = bg::read_manifest(out);
    if (exists && !add && segment_only.empty() && !m.segments.empty()) throw std::runtime_error(out + " already holds a store: --add adds a segment to it");
    if (partitions_given && partitions != m.partitions) throw std::runtime_error("background: the store has " + std::to_string(m.partitions) + " partitions");
    check_new_genomes(m, genomes);   // names now; sequences once counted
    for (const auto& s : m.segments)
        for (const auto& g : s.genomes)
            for (const auto& n : genomes)
                if (fs::path(g.path) == fs::path(n.path)) throw std::runtime_error("background: " + n.path + " is already counted as " + g.name);
    // --segment-only claims its folder: it must not exist; a job lock inside keeps a second job out
    std::unique_ptr<bg::JobLock> claim;
    bool claimed = false;
    if (!segment_only.empty()) {
        std::error_code ec;
        if (!fs::create_directory(fs::absolute(segment_only), ec) || ec)
            throw std::runtime_error("background: --segment-only " + segment_only + " exists already (or cannot be made); give a new folder");
        claimed = true;
        try { claim = std::make_unique<bg::JobLock>(fs::absolute(segment_only).string()); }
        catch (...) { fs::remove(fs::absolute(segment_only), ec); throw; }
    }

    // anchors: the store's references, plus a new one among these genomes
    auto refs = load_references(out, m);
    std::vector<bg::AnchorSet> sets;
    for (std::size_t r = 0; r < refs.size(); ++r) sets.push_back({m.references[r].first, refs[r].get()});
    std::unique_ptr<anchors::ReferenceAnchors> new_ref;
    Staging st{segment_only.empty() ? fs::path(out) / ".staging" / random_id() : fs::path(), nullptr};
    if (!st.dir.empty()) st.lock = std::make_unique<bg::JobLock>(st.dir.string());
    if (!coordinate_reference.empty()) {
        if (!segment_only.empty()) throw std::runtime_error("--segment-only builds with the store's coordinate references; add a new one to the store itself");
        for (const auto& r : m.references) if (r.first == coordinate_reference) throw std::runtime_error("background: " + coordinate_reference + " is already a coordinate reference");
        const auto it = std::find_if(genomes.begin(), genomes.end(), [&](const bg::Genome& g) { return g.name == coordinate_reference; });
        if (it == genomes.end()) throw std::runtime_error("--coordinate-reference " + coordinate_reference + " is not among the genomes given");
        if (is_gfa(it->path)) throw std::runtime_error("--coordinate-reference: a coordinate reference must be a FASTA");
        new_ref = std::make_unique<anchors::ReferenceAnchors>(anchors::build_reference(coordinate_reference, it->path, threads, [&](const std::string& s) { log.info(s); }));
        anchors::save(*new_ref, (st.dir / "reference.anc").string());
        sets.push_back({coordinate_reference, new_ref.get()});
        // committed genomes must be anchored on it too
        for (const auto& s : m.segments)
            map_segment(bg::read_segment_info((fs::path(out) / s.dir).string()), *new_ref, (st.dir / (s.dir + ".anchors")).string(), threads, log);
    }

    bg::BuildOptions o;
    o.genomes = genomes;
    o.segment_dir = segment_only.empty() ? (st.dir / "seg").string() : fs::absolute(segment_only).string();
    const std::string tmp_root = default_tmp(tmp_dir, segment_only.empty() ? st.dir.string() : fs::path(o.segment_dir).parent_path().string());
    o.tmp_dir = tmp_root;
    o.partitions = m.partitions;
    o.threads = threads;
    o.memory_bytes = memory;
    o.dedup = dedup;
    o.anchor_sets = sets;
    o.log = [&](const std::string& s) { log.info(s); };
    o.progress = !quiet;
    bg::SegmentStats stats;
    try {
        stats = bg::build_segment(o);
    } catch (...) {
        std::error_code ec;
        if (claimed) { claim.reset(); fs::remove_all(o.segment_dir, ec); }   // only the folder this job made
        throw;
    }
    const std::string summary = "distinct 31-mers: " + std::to_string(stats.distinct) + "; occurrences: " + std::to_string(stats.kmers) + "; size: " +
        human_bytes(stats.bytes) + "; spill: " + human_bytes(stats.spill_bytes) + " in " + std::to_string(stats.spill_files) + " files" +
        (stats.run_bytes > 0 ? "; sorted runs: " + human_bytes(stats.run_bytes) : "");
    if (!segment_only.empty()) {
        fail_point("segment-built");
        { std::error_code ec; claim.reset(); fs::remove(fs::path(o.segment_dir) / "job.lock", ec); }
        log.info("segment: " + o.segment_dir + "; " + summary + " (commit it with --attach)");
        log.wrote({o.segment_dir});
        log.done();
        return 0;
    }
    fail_point("count-before-commit");
    {
        // References added to the store meanwhile are mapped in this segment; with a new reference, segments
        // committed meanwhile are mapped on it. Then the commit is tried again.
        std::unique_ptr<bg::StoreLock> lock;
        bg::Manifest now;
        bg::SegmentInfo mine = bg::read_segment_info(o.segment_dir);
        for (int round = 0;; ++round) {
            lock = std::make_unique<bg::StoreLock>(out);
            now = bg::read_manifest(out);
            check_new_genomes(now, mine.genomes);
            std::vector<std::pair<std::string, std::uint64_t>> missing;
            for (const auto& r : now.references) {
                const bool have = std::any_of(mine.anchor_sets.begin(), mine.anchor_sets.end(), [&](const auto& x) { return x == r; });
                if (!have) missing.push_back(r);
            }
            if (new_ref)
                for (const auto& r : now.references) if (r.first == coordinate_reference) throw std::runtime_error("background: another job added " + coordinate_reference);
            bool unmapped = false;
            if (new_ref) for (const auto& s : now.segments) unmapped = unmapped || !fs::exists(st.dir / (s.dir + ".anchors"));
            if (missing.empty() && !unmapped) break;
            if (round == 20) throw std::runtime_error("background: the store keeps changing; run again later");
            lock.reset();
            for (const auto& [name, fp] : missing) {
                const anchors::ReferenceAnchors r = anchors::load_reference(ref_path(out, name).string());
                if (r.fingerprint() != fp) throw std::runtime_error("background: the anchor set of " + name + " does not match the manifest");
                map_segment(mine, r, (fs::path(o.segment_dir) / ("anchors-" + name + ".bin")).string(), threads, log);
                mine.anchor_sets.emplace_back(name, fp);
                bg::write_segment_info(o.segment_dir, mine);
                log.info("anchors: " + name + " (added to the store meanwhile); mapped in this segment");
            }
            if (new_ref) map_segments(out, now.segments, *new_ref, st.dir, ".anchors", threads, log);
        }
        if (new_ref) {
            for (const auto& s : now.segments) {
                const fs::path sd = fs::path(out) / s.dir;
                fs::rename(st.dir / (s.dir + ".anchors"), sd / ("anchors-" + coordinate_reference + ".bin"));
                add_anchor_record(sd, coordinate_reference, new_ref->fingerprint());
            }
            fs::create_directories(ref_path(out, coordinate_reference).parent_path());
            fs::rename(st.dir / "reference.anc", ref_path(out, coordinate_reference));
            now.references.emplace_back(coordinate_reference, new_ref->fingerprint());
        }
        const std::string dir = take_segment_name(out, now);
        fs::rename(o.segment_dir, fs::path(out) / dir);
        now.segments.push_back({dir, mine.distinct, mine.occurrences, mine.genomes});
        ++now.generation;
        fail_point("commit-before-manifest");
        bg::write_manifest(out, now);
        log.info("segment: " + dir + "; " + summary);
    }
    log.wrote({(fs::path(out) / "manifest.json").string()});
    log.done();
    return 0;
}

} // namespace panvar
