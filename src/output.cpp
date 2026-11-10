#include "panvar/output.hpp"

#include "panvar/cli_utils.hpp"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace panvar {

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::string current;
    current.reserve(line.size());

    bool in_quotes = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
                current.push_back('"');
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
            continue;
        }
        if (c == ',' && !in_quotes) {
            fields.push_back(current);
            current.clear();
            continue;
        }
        current.push_back(c);
    }
    // An unterminated quote is a truncated or hand-damaged line: read on, it would give a shorter list
    // of interior nodes.
    if (in_quotes) throw std::runtime_error("unterminated quoted field");
    fields.push_back(current);
    return fields;
}

namespace {

// close_output names the destination, not a staging file; its hint is for a failed write, so a
// failed open gives the system's reason instead (permission denied, no such directory).
[[noreturn]] void fail_open(std::ofstream& out, const std::string& path, int err) {
    try {
        cli::close_output(out, path);
    } catch (const std::runtime_error& e) {
        std::string msg = e.what();
        const std::size_t hint = msg.find(" (disk full or I/O error?)");
        if (err != 0 && hint != std::string::npos)
            msg.replace(hint, std::string::npos, std::string(": ") + std::strerror(err));
        throw std::runtime_error(msg);
    }
    throw std::runtime_error("Failed writing " + path);
}

std::string join_nodes(const std::vector<std::string>& nodes) {
    std::ostringstream oss;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (i > 0) {
            oss << ';';
        }
        oss << nodes[i];
    }
    return oss.str();
}

std::vector<std::string> split_semicolon(const std::string& input) {
    std::vector<std::string> parts;
    std::string item;
    std::stringstream ss(input);
    while (std::getline(ss, item, ';')) {
        if (!item.empty()) {
            parts.push_back(item);
        }
    }
    return parts;
}

// Digits only, as stoull alone would read "1abc" as 1 and "-1" as 2^64 - 1.
std::size_t parse_size_field(const std::string& value, const std::string& field_name) {
    const bool digits = !value.empty() && value.size() <= 19 &&
                        std::all_of(value.begin(), value.end(),
                                    [](char c) { return std::isdigit(static_cast<unsigned char>(c)); });
    if (!digits) throw std::runtime_error("Invalid numeric field '" + field_name + "': " + value);
    return static_cast<std::size_t>(std::stoull(value));
}

// The writer's field and list separators, and the quote, cannot appear in a node id: the CSV has no
// escaping for them inside inside_nodes, so "n;2" was read back as two nodes, "n" and "2".
void check_csv_node_id(const std::string& id) {
    if (id.find_first_of(",;\"") != std::string::npos)
        throw std::runtime_error("node id '" + id + "' contains ',', ';' or '\"', which the bubbles CSV "
                                 "cannot represent");
}

void ensure_parent_dir(const std::string& output_path) {
    const std::filesystem::path p(output_path);
    const auto parent = p.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
}

std::size_t required_column(
    const std::unordered_map<std::string, std::size_t>& index_by_name,
    const std::string& name) {
    const auto it = index_by_name.find(name);
    if (it == index_by_name.end()) {
        throw std::runtime_error("Missing required bubbles CSV column: " + name);
    }
    return it->second;
}

std::optional<std::size_t> optional_column(
    const std::unordered_map<std::string, std::size_t>& index_by_name,
    const std::string& name) {
    const auto it = index_by_name.find(name);
    if (it == index_by_name.end()) {
        return std::nullopt;
    }
    return it->second;
}

} // namespace

void write_bubbles_csv(
    const std::string& output_path,
    const std::vector<Bubble>& bubbles) {

    ensure_parent_dir(output_path);
    errno = 0;
    std::ofstream out(output_path);
    if (!out) fail_open(out, output_path, errno);

    // path_support is TRAVERSAL support -- how many paths cross the site at all, which on a fully-typed
    // panel is close to the panel size for most bubbles. The allele columns say what those traversals
    // contain, which is what "is this variant supported" actually means. min/max_inside_bp are interior
    // SPAN, not the size of the difference between alleles.
    for (const Bubble& b : bubbles) {
        check_csv_node_id(b.source);
        check_csv_node_id(b.sink);
        for (const std::string& n : b.inside) check_csv_node_id(n);
    }
    out << "bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,"
           "distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,"
           "min_inside_bp,max_inside_bp,inside_nodes,kind,parent_id\n";

    for (const auto& bubble : bubbles) {
        const std::size_t total_nodes = bubble.inside.size() + 2;
        std::vector<std::string> inside = bubble.inside;
        std::sort(inside.begin(), inside.end());

        out << bubble.id << ','
            << bubble.source << ','
            << (bubble.source_reverse ? '-' : '+') << ','
            << bubble.sink << ','
            << (bubble.sink_reverse ? '-' : '+') << ','
            << bubble.inside.size() << ','
            << total_nodes << ','
            << bubble.path_support << ','
            << bubble.distinct_alleles << ','
            << bubble.ref_allele_support << ','
            << bubble.alt_allele_support_max << ','
            << bubble.alt_allele_support_min << ','
            << bubble.min_inside_bp << ','
            << bubble.max_inside_bp << ','
            << '"' << join_nodes(inside) << '"' << ','
            << (bubble.inversion_parent ? "inversion" : "site") << ','
            << (bubble.parent_id == 0 ? std::string() : std::to_string(bubble.parent_id))
            << '\n';
    }
    cli::close_output(out, output_path);
}

std::vector<Bubble> read_bubbles_csv(const std::string& input_path, bool with_inversion_sites) {
    std::error_code dir_ec;
    if (std::filesystem::is_directory(input_path, dir_ec))
        throw std::runtime_error("Bubbles CSV " + input_path + " is a directory, not a file");
    std::ifstream in(input_path);
    if (!in) {
        throw std::runtime_error("Failed to read bubbles CSV: " + input_path);
    }
    // Every module writes the table plain; a gzipped one would otherwise fail as a damaged quoted field.
    if (in.peek() == 0x1f) {
        char magic[2] = {0, 0};
        in.read(magic, 2);
        if (static_cast<unsigned char>(magic[1]) == 0x8b)
            throw std::runtime_error("Bubbles CSV " + input_path + " is gzip-compressed; decompress it first");
        in.clear();
        in.seekg(0);
    }

    std::string header_line;
    if (!std::getline(in, header_line)) {
        throw std::runtime_error("Bubbles CSV is empty: " + input_path);
    }
    if (!header_line.empty() && header_line.back() == '\r') {
        header_line.pop_back();
    }

    const std::vector<std::string> header_fields = split_csv_line(header_line);
    std::unordered_map<std::string, std::size_t> index_by_name;
    index_by_name.reserve(header_fields.size() * 2);
    for (std::size_t i = 0; i < header_fields.size(); ++i) {
        index_by_name[header_fields[i]] = i;
    }

    const std::size_t idx_id = required_column(index_by_name, "bubble_id");
    const std::size_t idx_source = required_column(index_by_name, "source");
    const std::size_t idx_sink = required_column(index_by_name, "sink");
    const std::size_t idx_path_support = required_column(index_by_name, "path_support");
    const std::size_t idx_min_bp = required_column(index_by_name, "min_inside_bp");
    const std::size_t idx_max_bp = required_column(index_by_name, "max_inside_bp");
    const auto idx_src_orient = optional_column(index_by_name, "source_orient");
    const auto idx_snk_orient = optional_column(index_by_name, "sink_orient");
    const auto idx_distinct = optional_column(index_by_name, "distinct_alleles");
    const auto idx_ref_sup = optional_column(index_by_name, "ref_allele_support");
    const auto idx_alt_max = optional_column(index_by_name, "alt_allele_support_max");
    const auto idx_alt_min = optional_column(index_by_name, "alt_allele_support_min");
    const auto idx_long_path = optional_column(index_by_name, "long_path_support");
    const auto idx_inversion = optional_column(index_by_name, "inversion_signal");
    const std::size_t idx_inside_nodes = required_column(index_by_name, "inside_nodes");
    const auto idx_kind = optional_column(index_by_name, "kind");
    const auto idx_parent = optional_column(index_by_name, "parent_id");

    std::vector<Bubble> bubbles;
    std::string line;
    std::size_t line_no = 1;
    while (std::getline(in, line)) {
        ++line_no;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }

        std::vector<std::string> fields;
        try {
            fields = split_csv_line(line);
        } catch (const std::exception& e) {
            throw std::runtime_error("Malformed bubbles CSV line " + std::to_string(line_no) + ": " + e.what());
        }
        const auto require_field = [&](std::size_t idx, const std::string& name) -> const std::string& {
            if (idx >= fields.size()) {
                throw std::runtime_error(
                    "Malformed bubbles CSV line " + std::to_string(line_no) + ": missing field " + name);
            }
            return fields[idx];
        };

        Bubble bubble;
        bubble.id = parse_size_field(require_field(idx_id, "bubble_id"), "bubble_id");
        // Site ids start at 1: 0 stands for "no site" (an empty parent_id, a record outside any site).
        if (bubble.id == 0)
            throw std::runtime_error("Malformed bubbles CSV line " + std::to_string(line_no) +
                                     ": bubble_id 0; site ids start at 1");
        bubble.source = require_field(idx_source, "source");
        bubble.sink = require_field(idx_sink, "sink");
        // Absent in a CSV written before boundaries carried orientation; forward is the safe default
        // because that is what an unoriented pair was always implicitly read as.
        const auto read_orient = [&](const std::optional<std::size_t>& idx, const char* name) {
            if (!idx.has_value() || *idx >= fields.size() || fields[*idx].empty() || fields[*idx] == "+")
                return false;
            if (fields[*idx] == "-") return true;
            throw std::runtime_error("Malformed bubbles CSV line " + std::to_string(line_no) + ": " + name +
                                     " is '" + fields[*idx] + "', expected '+' or '-'");
        };
        bubble.source_reverse = read_orient(idx_src_orient, "source_orient");
        bubble.sink_reverse = read_orient(idx_snk_orient, "sink_orient");
        bubble.orientation_known = idx_src_orient.has_value() && idx_snk_orient.has_value();
        bubble.path_support = parse_size_field(require_field(idx_path_support, "path_support"), "path_support");
        bubble.min_inside_bp = parse_size_field(require_field(idx_min_bp, "min_inside_bp"), "min_inside_bp");
        bubble.max_inside_bp = parse_size_field(require_field(idx_max_bp, "max_inside_bp"), "max_inside_bp");
        // Allele support counts, when the table has them, so a consumer can re-apply --min-alt-support.
        const auto read_optional = [&](const std::optional<std::size_t>& idx, const char* name,
                                       std::size_t& target) {
            if (idx.has_value() && *idx < fields.size() && !fields[*idx].empty()) {
                target = parse_size_field(fields[*idx], name);
            }
        };
        read_optional(idx_distinct, "distinct_alleles", bubble.distinct_alleles);
        read_optional(idx_ref_sup, "ref_allele_support", bubble.ref_allele_support);
        read_optional(idx_alt_max, "alt_allele_support_max", bubble.alt_allele_support_max);
        read_optional(idx_alt_min, "alt_allele_support_min", bubble.alt_allele_support_min);
        if (idx_long_path.has_value()) {
            bubble.long_path_support =
                parse_size_field(require_field(*idx_long_path, "long_path_support"), "long_path_support");
        } else {
            bubble.long_path_support = 0;
        }
        if (idx_inversion.has_value()) {
            const std::string inversion = require_field(*idx_inversion, "inversion_signal");
            bubble.inversion_signal = (inversion == "1" || inversion == "true" || inversion == "TRUE");
        } else {
            bubble.inversion_signal = false;
        }
        bubble.inside = split_semicolon(require_field(idx_inside_nodes, "inside_nodes"));
        std::sort(bubble.inside.begin(), bubble.inside.end());
        if (idx_kind.has_value() && *idx_kind < fields.size() && !fields[*idx_kind].empty()) {
            const std::string& kind = fields[*idx_kind];
            if (kind == "inversion") bubble.inversion_parent = true;
            else if (kind != "site")
                throw std::runtime_error("Malformed bubbles CSV line " + std::to_string(line_no) + ": kind is '" +
                                         kind + "', expected 'site' or 'inversion'");
        }
        read_optional(idx_parent, "parent_id", bubble.parent_id);
        if (bubble.inversion_parent && !with_inversion_sites) continue;

        bubbles.push_back(std::move(bubble));
    }

    std::sort(bubbles.begin(), bubbles.end(), [](const Bubble& lhs, const Bubble& rhs) {
        return lhs.id < rhs.id;
    });

    return bubbles;
}

void write_bandage_node_colors_csv(
    const std::string& output_path,
    const std::vector<Bubble>& retained_bubbles,
    const std::vector<Bubble>& non_snp_bubbles) {

    ensure_parent_dir(output_path);
    errno = 0;
    std::ofstream out(output_path);
    if (!out) fail_open(out, output_path, errno);

    // Only interiors are coloured: red for kept sites, blue for the other candidates, red winning on a node
    // in both. Boundaries stay uncoloured, so long backbone nodes shared by neighbouring sites do not read as
    // variation.
    constexpr const char* kRed = "#E41A1C";
    constexpr const char* kBlue = "#377EB8";
    std::unordered_map<std::string, int> priority;   // 2 red, 1 blue
    const auto paint = [&](const std::string& node, int p) {
        int& have = priority[node];
        have = std::max(have, p);
    };
    for (const auto& bubble : non_snp_bubbles)
        for (const auto& node : bubble.inside) paint(node, 1);
    for (const auto& bubble : retained_bubbles) {
        if (bubble.inversion_parent) continue;   // it would paint over the sites inside it
        for (const auto& node : bubble.inside) paint(node, 2);
    }

    std::vector<std::string> node_ids;
    node_ids.reserve(priority.size());
    for (const auto& kv : priority) node_ids.push_back(kv.first);
    std::sort(node_ids.begin(), node_ids.end());

    out << "Name,Colour\n";
    for (const auto& node_id : node_ids) {
        const int p = priority[node_id];
        out << node_id << ',' << (p == 2 ? kRed : kBlue) << '\n';
    }
    cli::close_output(out, output_path);
}

void write_snarl_debug_tsv(
    const std::string& output_path,
    const std::vector<SnarlDebugEntry>& entries) {

    ensure_parent_dir(output_path);
    errno = 0;
    std::ofstream out(output_path);
    if (!out) fail_open(out, output_path, errno);

    out << "candidate_id\tsource\tsink\tinside_node_count\tn_paths\tmin_inside_bp\t"
        << "long_path_support\tinversion_signal\taccepted\n";

    for (const auto& e : entries) {
        out << e.candidate_id << '\t'
            << e.source << '\t'
            << e.sink << '\t'
            << e.inside_node_count << '\t'
            << e.n_paths << '\t'
            << e.min_inside_bp << '\t'
            << e.long_path_support << '\t'
            << (e.inversion_signal ? 1 : 0) << '\t'
            << (e.accepted ? 1 : 0) << '\n';
    }
    cli::close_output(out, output_path);
}

} // namespace panvar
