#include "panvar/output.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace panvar {
namespace {

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

std::string color_for_level(const std::unordered_map<std::size_t, std::string>& palette, std::size_t level) {
    if (palette.empty()) {
        return "#808080";
    }
    auto it = palette.find(level);
    if (it != palette.end()) {
        return it->second;
    }

    std::size_t max_level = 0;
    for (const auto& [k, _v] : palette) {
        max_level = std::max(max_level, k);
    }
    if (max_level == 0) {
        return "#808080";
    }

    const std::size_t folded = ((level - 1) % max_level) + 1;
    auto fallback = palette.find(folded);
    return fallback == palette.end() ? "#808080" : fallback->second;
}

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
    fields.push_back(current);
    return fields;
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

std::size_t parse_size_field(const std::string& value, const std::string& field_name) {
    try {
        return static_cast<std::size_t>(std::stoull(value));
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid numeric field '" + field_name + "': " + value);
    }
}

BubbleType parse_bubble_type_field(const std::string& value) {
    if (value == "simple") {
        return BubbleType::Simple;
    }
    if (value == "insertion") {
        return BubbleType::Insertion;
    }
    if (value == "super") {
        return BubbleType::Super;
    }
    throw std::runtime_error("Invalid bubble type in CSV: " + value);
}

SiteMode parse_site_mode_field(const std::string& value) {
    if (value == "snarl") {
        return SiteMode::Snarl;
    }
    if (value == "superbubble") {
        return SiteMode::Superbubble;
    }
    throw std::runtime_error("Invalid site_mode in CSV: " + value);
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

} // namespace

std::string site_mode_to_string(SiteMode mode) {
    switch (mode) {
        case SiteMode::Snarl:
            return "snarl";
        case SiteMode::Superbubble:
        default:
            return "superbubble";
    }
}

std::unordered_map<std::size_t, std::string> level_palette() {
    return {
        {1, "#E41A1C"},
        {2, "#377EB8"},
        {3, "#4DAF4A"},
        {4, "#FF7F00"},
        {5, "#A65628"},
        {6, "#F781BF"},
        {7, "#999999"},
        {8, "#66C2A5"},
        {9, "#FC8D62"},
        {10, "#8DA0CB"},
        {11, "#E78AC3"},
        {12, "#A6D854"},
    };
}

void write_bubbles_csv(
    const std::string& output_path,
    const std::vector<Bubble>& bubbles) {

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write bubbles CSV: " + output_path);
    }

    out << "bubble_id,source,sink,site_mode,type,inside_node_count,total_node_count,nesting_level,parent_id,path_support,min_inside_bp,max_inside_bp,long_path_support,inversion_signal,inside_nodes\n";

    for (const auto& bubble : bubbles) {
        const std::size_t total_nodes = bubble.inside.size() + 2;
        std::vector<std::string> inside = bubble.inside;
        std::sort(inside.begin(), inside.end());

        out << bubble.id << ','
            << bubble.source << ','
            << bubble.sink << ','
            << site_mode_to_string(bubble.site_mode) << ','
            << bubble_type_to_string(bubble.type) << ','
            << bubble.inside.size() << ','
            << total_nodes << ','
            << bubble.nesting_level << ','
            << bubble.parent_id << ','
            << bubble.path_support << ','
            << bubble.min_inside_bp << ','
            << bubble.max_inside_bp << ','
            << bubble.long_path_support << ','
            << (bubble.inversion_signal ? 1 : 0) << ','
            << '"' << join_nodes(inside) << '"'
            << '\n';
    }
}

std::vector<Bubble> read_bubbles_csv(const std::string& input_path) {
    std::ifstream in(input_path);
    if (!in) {
        throw std::runtime_error("Failed to read bubbles CSV: " + input_path);
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
    const std::size_t idx_site_mode = required_column(index_by_name, "site_mode");
    const std::size_t idx_type = required_column(index_by_name, "type");
    const std::size_t idx_nesting = required_column(index_by_name, "nesting_level");
    const std::size_t idx_parent = required_column(index_by_name, "parent_id");
    const std::size_t idx_path_support = required_column(index_by_name, "path_support");
    const std::size_t idx_min_bp = required_column(index_by_name, "min_inside_bp");
    const std::size_t idx_max_bp = required_column(index_by_name, "max_inside_bp");
    const std::size_t idx_long_path = required_column(index_by_name, "long_path_support");
    const std::size_t idx_inversion = required_column(index_by_name, "inversion_signal");
    const std::size_t idx_inside_nodes = required_column(index_by_name, "inside_nodes");

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

        const std::vector<std::string> fields = split_csv_line(line);
        const auto require_field = [&](std::size_t idx, const std::string& name) -> const std::string& {
            if (idx >= fields.size()) {
                throw std::runtime_error(
                    "Malformed bubbles CSV line " + std::to_string(line_no) + ": missing field " + name);
            }
            return fields[idx];
        };

        Bubble bubble;
        bubble.id = parse_size_field(require_field(idx_id, "bubble_id"), "bubble_id");
        bubble.source = require_field(idx_source, "source");
        bubble.sink = require_field(idx_sink, "sink");
        bubble.site_mode = parse_site_mode_field(require_field(idx_site_mode, "site_mode"));
        bubble.type = parse_bubble_type_field(require_field(idx_type, "type"));
        bubble.nesting_level = parse_size_field(require_field(idx_nesting, "nesting_level"), "nesting_level");
        bubble.parent_id = parse_size_field(require_field(idx_parent, "parent_id"), "parent_id");
        bubble.path_support = parse_size_field(require_field(idx_path_support, "path_support"), "path_support");
        bubble.min_inside_bp = parse_size_field(require_field(idx_min_bp, "min_inside_bp"), "min_inside_bp");
        bubble.max_inside_bp = parse_size_field(require_field(idx_max_bp, "max_inside_bp"), "max_inside_bp");
        bubble.long_path_support =
            parse_size_field(require_field(idx_long_path, "long_path_support"), "long_path_support");
        const std::string inversion = require_field(idx_inversion, "inversion_signal");
        bubble.inversion_signal = (inversion == "1" || inversion == "true" || inversion == "TRUE");
        bubble.inside = split_semicolon(require_field(idx_inside_nodes, "inside_nodes"));
        std::sort(bubble.inside.begin(), bubble.inside.end());

        bubbles.push_back(std::move(bubble));
    }

    std::sort(bubbles.begin(), bubbles.end(), [](const Bubble& lhs, const Bubble& rhs) {
        return lhs.id < rhs.id;
    });

    return bubbles;
}

void write_bandage_node_colors_csv(
    const std::string& output_path,
    const std::vector<Bubble>& bubbles) {

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write Bandage color CSV: " + output_path);
    }

    const auto palette = level_palette();

    struct NodeColorInfo {
        std::size_t max_level = 0;
    };

    std::unordered_map<std::string, NodeColorInfo> by_node;

    for (const auto& bubble : bubbles) {
        auto nodes = nodes_in_bubble(bubble);
        for (const auto& node : nodes) {
            auto& info = by_node[node];
            info.max_level = std::max(info.max_level, bubble.nesting_level);
        }
    }

    std::vector<std::string> node_ids;
    node_ids.reserve(by_node.size());
    for (const auto& [node_id, _info] : by_node) {
        node_ids.push_back(node_id);
    }
    std::sort(node_ids.begin(), node_ids.end());

    out << "Name,Colour\n";
    for (const auto& node_id : node_ids) {
        auto& info = by_node[node_id];
        out << node_id << ',' << color_for_level(palette, info.max_level) << '\n';
    }
}

void write_snarl_debug_tsv(
    const std::string& output_path,
    const std::vector<SnarlDebugEntry>& entries) {

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write snarl debug TSV: " + output_path);
    }

    out << "candidate_id\tsource\tsink\tinside_node_count\tn_paths\tmin_inside_bp\tnested\taccepted\treason\n";

    for (const auto& e : entries) {
        out << e.candidate_id << '\t'
            << e.source << '\t'
            << e.sink << '\t'
            << e.inside_node_count << '\t'
            << e.n_paths << '\t'
            << e.min_inside_bp << '\t'
            << (e.nested ? 1 : 0) << '\t'
            << (e.accepted ? 1 : 0) << '\t'
            << e.reason << '\n';
    }
}

} // namespace panvar
