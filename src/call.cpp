#include "panvar/call.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <regex>
#include <sstream>
#include <string_view>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <zlib.h>

#include "call_internal.hpp"
#include "panvar/output.hpp"

namespace panvar {
namespace {

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

std::vector<std::size_t> split_semicolon_size_t(const std::string& value) {
    std::vector<std::size_t> out;
    if (value.empty()) {
        return out;
    }
    std::stringstream ss(value);
    std::string token;
    while (std::getline(ss, token, ';')) {
        if (token.empty()) {
            continue;
        }
        try {
            out.push_back(static_cast<std::size_t>(std::stoull(token)));
        } catch (const std::exception&) {
            throw std::runtime_error("Invalid allele id in member_alleles: " + token);
        }
    }
    return out;
}

struct CandidateInterval {
    std::size_t left = 0;
    std::size_t right = 0;
    std::size_t inside_count = 0;
    bool source_to_sink = true;
};

std::unordered_map<std::string, const PathRecord*> path_by_name(const Graph& graph) {
    std::unordered_map<std::string, const PathRecord*> out;
    out.reserve(graph.paths.size() * 2);
    for (const auto& path : graph.paths) {
        out[path.name] = &path;
    }
    return out;
}

std::vector<std::size_t> prefix_bp(
    const PathRecord& path,
    const std::unordered_map<std::string, Node>& nodes) {

    std::vector<std::size_t> pref(path.steps.size() + 1, 0);
    for (std::size_t i = 0; i < path.steps.size(); ++i) {
        const auto it = nodes.find(path.steps[i].node_id);
        const std::size_t len = (it == nodes.end()) ? 1 : std::max<std::size_t>(1, it->second.sequence.size());
        pref[i + 1] = pref[i] + len;
    }
    return pref;
}

std::string reverse_complement(const std::string& sequence) {
    std::string rc;
    rc.reserve(sequence.size());
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        const char c = *it;
        switch (c) {
            case 'A':
            case 'a':
                rc.push_back('T');
                break;
            case 'C':
            case 'c':
                rc.push_back('G');
                break;
            case 'G':
            case 'g':
                rc.push_back('C');
                break;
            case 'T':
            case 't':
                rc.push_back('A');
                break;
            default:
                rc.push_back('N');
                break;
        }
    }
    return rc;
}

std::string spell_sequence(
    const Graph& graph,
    const std::vector<PathStep>& steps,
    bool& complete) {

    complete = true;
    std::size_t total_len = 0;
    for (const auto& step : steps) {
        const auto node_it = graph.nodes.find(step.node_id);
        if (node_it == graph.nodes.end() || node_it->second.sequence.empty()) {
            complete = false;
            return "";
        }
        total_len += node_it->second.sequence.size();
    }

    std::string seq;
    seq.reserve(total_len);
    for (const auto& step : steps) {
        const auto& node_seq = graph.nodes.at(step.node_id).sequence;
        if (!step.reverse) {
            seq += node_seq;
        } else {
            seq += reverse_complement(node_seq);
        }
    }
    return seq;
}

std::vector<PathStep> canonical_steps_for_bubble(
    const PathRecord& path,
    const Bubble& bubble,
    const CandidateInterval& interval) {

    if (interval.left >= path.steps.size() || interval.right >= path.steps.size() || interval.left > interval.right) {
        return {};
    }

    std::vector<PathStep> out;
    out.reserve(interval.right - interval.left + 1);

    if (interval.source_to_sink) {
        for (std::size_t i = interval.left; i <= interval.right; ++i) {
            out.push_back(path.steps[i]);
        }
        return out;
    }

    for (std::size_t i = interval.right + 1; i > interval.left; --i) {
        const PathStep& s = path.steps[i - 1];
        out.push_back(PathStep{s.node_id, !s.reverse});
    }

    if (out.empty()) {
        return out;
    }

    if (out.front().node_id != bubble.source || out.back().node_id != bubble.sink) {
        std::reverse(out.begin(), out.end());
        for (auto& step : out) {
            step.reverse = !step.reverse;
        }
    }

    return out;
}

std::string build_walk_signature(const std::vector<PathStep>& steps) {
    std::string sig;
    sig.reserve(steps.size() * 8);
    for (std::size_t i = 0; i < steps.size(); ++i) {
        if (i > 0) {
            sig.push_back(',');
        }
        sig += steps[i].node_id;
        sig.push_back(steps[i].reverse ? '-' : '+');
    }
    return sig;
}

std::uint64_t hash_step_token(const PathStep& step) {
    std::uint64_t h = 1469598103934665603ULL;
    for (const unsigned char c : step.node_id) {
        h ^= static_cast<std::uint64_t>(c);
        h *= 1099511628211ULL;
    }
    h ^= (step.reverse ? 0xF0ULL : 0x0FULL);
    h *= 1099511628211ULL;
    return h;
}

std::vector<std::uint64_t> build_walk_tokens(const std::vector<PathStep>& steps) {
    std::vector<std::uint64_t> tokens;
    tokens.reserve(steps.size());
    for (const auto& step : steps) {
        tokens.push_back(hash_step_token(step));
    }
    return tokens;
}

double elapsed_seconds(const std::chrono::steady_clock::time_point& start) {
    const auto now = std::chrono::steady_clock::now();
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
    return static_cast<double>(ms) / 1000.0;
}

std::size_t parse_csv_size(const std::string& value, const std::string& field_name) {
    try {
        return static_cast<std::size_t>(std::stoull(value));
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid numeric field '" + field_name + "': " + value);
    }
}

bool parse_csv_bool(const std::string& value, const std::string& field_name) {
    if (value == "1" || value == "true" || value == "TRUE") {
        return true;
    }
    if (value == "0" || value == "false" || value == "FALSE") {
        return false;
    }
    throw std::runtime_error("Invalid boolean field '" + field_name + "': " + value);
}

struct PangeneGeneRecord {
    std::size_t start0 = 0;
    std::size_t end0 = 0;
    std::string gene_name;
    char orientation = '.';
};

using PangeneGeneIndex = std::unordered_map<std::string, std::vector<PangeneGeneRecord>>;
using GeneCopyCounts = std::unordered_map<std::string, std::size_t>;

struct CompiledGeneFilter {
    std::string pattern;
    bool use_regex = false;
    std::regex regex;
    std::string lowercase_pattern;
};

struct PangeneDeltaSummary {
    std::string cn_delta;
    std::string gain_genes;
    std::string loss_genes;
    std::size_t gain_copies = 0;
    std::size_t loss_copies = 0;
};

std::string ascii_lower(const std::string& input) {
    std::string out = input;
    std::transform(
        out.begin(),
        out.end(),
        out.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return out;
}

std::vector<std::string> split_tab_fields(std::string_view line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
        const std::size_t tab = line.find('\t', start);
        if (tab == std::string_view::npos) {
            fields.emplace_back(line.substr(start));
            break;
        }
        fields.emplace_back(line.substr(start, tab - start));
        start = tab + 1;
    }
    return fields;
}

bool has_gz_suffix(const std::string& path) {
    return path.size() >= 3 && path.substr(path.size() - 3) == ".gz";
}

PangeneGeneIndex load_pangene_bed(const std::string& bed_path) {
    PangeneGeneIndex index;

    auto handle_fields = [&](const std::vector<std::string>& fields, std::size_t line_no) {
        if (fields.size() < 4) {
            throw std::runtime_error(
                "Invalid pangene BED line " + std::to_string(line_no) +
                ": expected at least 4 columns");
        }
        std::size_t start0 = 0;
        std::size_t end0 = 0;
        try {
            start0 = static_cast<std::size_t>(std::stoull(fields[2]));
            end0 = static_cast<std::size_t>(std::stoull(fields[3]));
        } catch (const std::exception&) {
            const std::string col0 = ascii_lower(fields[0]);
            const std::string col2 = ascii_lower(fields[2]);
            const std::string col3 = ascii_lower(fields[3]);
            const bool header_like =
                (col0 == "molecule" || col0 == "path" || col0 == "seqname") &&
                (col2 == "start") &&
                (col3 == "end");
            if (header_like) {
                return;
            }
            throw std::runtime_error(
                "Invalid pangene BED coordinates on line " + std::to_string(line_no));
        }
        if (end0 < start0) {
            std::swap(start0, end0);
        }
        PangeneGeneRecord rec;
        rec.start0 = start0;
        rec.end0 = end0;
        rec.gene_name = fields[1];
        if (fields.size() >= 5 && !fields[4].empty()) {
            rec.orientation = fields[4][0];
        }
        index[fields[0]].push_back(std::move(rec));
    };

    if (has_gz_suffix(bed_path)) {
        gzFile fp = gzopen(bed_path.c_str(), "rb");
        if (fp == nullptr) {
            throw std::runtime_error("Failed to open pangene BED: " + bed_path);
        }
        constexpr int kBuf = 1 << 20;
        std::vector<char> buf(static_cast<std::size_t>(kBuf), '\0');
        std::size_t line_no = 0;
        while (gzgets(fp, buf.data(), kBuf) != nullptr) {
            ++line_no;
            std::string line(buf.data());
            while (!line.empty() &&
                   (line.back() == '\n' || line.back() == '\r')) {
                line.pop_back();
            }
            if (line.empty() || line[0] == '#') {
                continue;
            }
            handle_fields(split_tab_fields(line), line_no);
        }
        gzclose(fp);
    } else {
        std::ifstream in(bed_path);
        if (!in) {
            throw std::runtime_error("Failed to open pangene BED: " + bed_path);
        }
        std::string line;
        std::size_t line_no = 0;
        while (std::getline(in, line)) {
            ++line_no;
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            if (line.empty() || line[0] == '#') {
                continue;
            }
            handle_fields(split_tab_fields(line), line_no);
        }
    }

    for (auto& [path_name, genes] : index) {
        (void)path_name;
        std::sort(genes.begin(), genes.end(), [](const PangeneGeneRecord& lhs, const PangeneGeneRecord& rhs) {
            if (lhs.start0 != rhs.start0) {
                return lhs.start0 < rhs.start0;
            }
            if (lhs.end0 != rhs.end0) {
                return lhs.end0 < rhs.end0;
            }
            return lhs.gene_name < rhs.gene_name;
        });
    }
    return index;
}

std::vector<CompiledGeneFilter> compile_gene_filters(
    const std::vector<std::string>& patterns) {

    std::vector<CompiledGeneFilter> compiled;
    compiled.reserve(patterns.size());
    for (const auto& pattern : patterns) {
        if (pattern.empty()) {
            continue;
        }
        CompiledGeneFilter item;
        item.pattern = pattern;
        try {
            item.regex = std::regex(pattern, std::regex::icase);
            item.use_regex = true;
        } catch (const std::regex_error&) {
            item.use_regex = false;
            item.lowercase_pattern = ascii_lower(pattern);
        }
        compiled.push_back(std::move(item));
    }
    return compiled;
}

bool gene_name_matches_filters(
    const std::string& gene_name,
    const std::vector<CompiledGeneFilter>& filters) {

    if (filters.empty()) {
        return true;
    }
    const std::string lower_name = ascii_lower(gene_name);
    for (const auto& filter : filters) {
        if (filter.use_regex) {
            if (std::regex_search(gene_name, filter.regex)) {
                return true;
            }
            continue;
        }
        if (!filter.lowercase_pattern.empty() &&
            lower_name.find(filter.lowercase_pattern) != std::string::npos) {
            return true;
        }
    }
    return false;
}

GeneCopyCounts count_pangene_copies_for_interval(
    const PangeneGeneIndex& pangene_index,
    const std::string& path_name,
    std::size_t start_bp,
    std::size_t end_bp,
    const std::vector<CompiledGeneFilter>& filters) {

    GeneCopyCounts counts;
    const auto it = pangene_index.find(path_name);
    if (it == pangene_index.end()) {
        return counts;
    }
    if (end_bp < start_bp) {
        std::swap(start_bp, end_bp);
    }

    for (const auto& gene : it->second) {
        if (gene.end0 <= start_bp) {
            continue;
        }
        if (gene.start0 >= end_bp) {
            break;
        }
        if (!gene_name_matches_filters(gene.gene_name, filters)) {
            continue;
        }
        counts[gene.gene_name] += 1;
    }
    return counts;
}

PangeneDeltaSummary summarize_gene_copy_delta(
    const GeneCopyCounts& ref_counts,
    const GeneCopyCounts& alt_counts) {

    PangeneDeltaSummary out;
    std::unordered_set<std::string> genes;
    genes.reserve(ref_counts.size() + alt_counts.size() + 1);
    for (const auto& kv : ref_counts) {
        genes.insert(kv.first);
    }
    for (const auto& kv : alt_counts) {
        genes.insert(kv.first);
    }

    std::vector<std::string> sorted_genes(genes.begin(), genes.end());
    std::sort(sorted_genes.begin(), sorted_genes.end());

    std::vector<std::string> delta_items;
    std::vector<std::string> gain_items;
    std::vector<std::string> loss_items;
    for (const auto& gene : sorted_genes) {
        const std::size_t ref_cn = (ref_counts.count(gene) > 0) ? ref_counts.at(gene) : 0;
        const std::size_t alt_cn = (alt_counts.count(gene) > 0) ? alt_counts.at(gene) : 0;
        if (ref_cn == alt_cn) {
            continue;
        }
        delta_items.push_back(
            gene + ":" + std::to_string(ref_cn) + ">" + std::to_string(alt_cn));
        if (alt_cn > ref_cn) {
            const std::size_t diff = alt_cn - ref_cn;
            out.gain_copies += diff;
            gain_items.push_back(gene + "(+" + std::to_string(diff) + ")");
        } else {
            const std::size_t diff = ref_cn - alt_cn;
            out.loss_copies += diff;
            loss_items.push_back(gene + "(-" + std::to_string(diff) + ")");
        }
    }

    auto join_csv = [](const std::vector<std::string>& items) -> std::string {
        std::string out_csv;
        for (std::size_t i = 0; i < items.size(); ++i) {
            if (i > 0) {
                out_csv.push_back(',');
            }
            out_csv += items[i];
        }
        return out_csv;
    };
    out.cn_delta = join_csv(delta_items);
    out.gain_genes = join_csv(gain_items);
    out.loss_genes = join_csv(loss_items);
    return out;
}

std::size_t required_csv_column(
    const std::unordered_map<std::string, std::size_t>& index_by_name,
    const std::string& column_name) {
    const auto it = index_by_name.find(column_name);
    if (it == index_by_name.end()) {
        throw std::runtime_error("Missing required CSV column: " + column_name);
    }
    return it->second;
}

std::string endpoint_key(const std::string& source, const std::string& sink) {
    return source + '\t' + sink;
}

std::vector<PrecomputedClusterRow> read_precomputed_clusters_csv(
    const std::string& input_path) {

    std::ifstream in(input_path);
    if (!in) {
        throw std::runtime_error("Failed to read precomputed clusters CSV: " + input_path);
    }

    std::string header_line;
    if (!std::getline(in, header_line)) {
        throw std::runtime_error("Precomputed clusters CSV is empty: " + input_path);
    }
    if (!header_line.empty() && header_line.back() == '\r') {
        header_line.pop_back();
    }

    const auto header_fields = split_csv_line(header_line);
    std::unordered_map<std::string, std::size_t> header_idx;
    header_idx.reserve(header_fields.size() * 2);
    for (std::size_t i = 0; i < header_fields.size(); ++i) {
        header_idx[header_fields[i]] = i;
    }

    const std::size_t idx_bubble = required_csv_column(header_idx, "bubble_id");
    const std::size_t idx_source = required_csv_column(header_idx, "source");
    const std::size_t idx_sink = required_csv_column(header_idx, "sink");
    const std::size_t idx_cluster = required_csv_column(header_idx, "cluster_id");
    const std::size_t idx_rep_allele = required_csv_column(header_idx, "representative_allele_id");
    const std::size_t idx_support = required_csv_column(header_idx, "total_path_support");
    const std::size_t idx_members = required_csv_column(header_idx, "member_alleles");

    std::vector<PrecomputedClusterRow> rows;
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
        const auto fields = split_csv_line(line);
        auto field_or_throw = [&](std::size_t idx, const std::string& name) -> const std::string& {
            if (idx >= fields.size()) {
                throw std::runtime_error(
                    "Malformed precomputed clusters CSV line " + std::to_string(line_no) +
                    ": missing field " + name);
            }
            return fields[idx];
        };

        PrecomputedClusterRow row;
        row.bubble_id = parse_csv_size(field_or_throw(idx_bubble, "bubble_id"), "bubble_id");
        row.source = field_or_throw(idx_source, "source");
        row.sink = field_or_throw(idx_sink, "sink");
        row.cluster_id = parse_csv_size(field_or_throw(idx_cluster, "cluster_id"), "cluster_id");
        row.representative_allele_id =
            parse_csv_size(field_or_throw(idx_rep_allele, "representative_allele_id"), "representative_allele_id");
        row.total_path_support =
            parse_csv_size(field_or_throw(idx_support, "total_path_support"), "total_path_support");
        row.member_alleles = split_semicolon_size_t(field_or_throw(idx_members, "member_alleles"));
        rows.push_back(std::move(row));
    }

    return rows;
}

std::vector<PrecomputedAssignmentRow> read_precomputed_assignments_csv(
    const std::string& input_path) {

    std::ifstream in(input_path);
    if (!in) {
        throw std::runtime_error("Failed to read precomputed assignments CSV: " + input_path);
    }

    std::string header_line;
    if (!std::getline(in, header_line)) {
        throw std::runtime_error("Precomputed assignments CSV is empty: " + input_path);
    }
    if (!header_line.empty() && header_line.back() == '\r') {
        header_line.pop_back();
    }

    const auto header_fields = split_csv_line(header_line);
    std::unordered_map<std::string, std::size_t> header_idx;
    header_idx.reserve(header_fields.size() * 2);
    for (std::size_t i = 0; i < header_fields.size(); ++i) {
        header_idx[header_fields[i]] = i;
    }

    const std::size_t idx_bubble = required_csv_column(header_idx, "bubble_id");
    const std::size_t idx_source = required_csv_column(header_idx, "source");
    const std::size_t idx_sink = required_csv_column(header_idx, "sink");
    const std::size_t idx_path = required_csv_column(header_idx, "path_name");
    const std::size_t idx_cluster = required_csv_column(header_idx, "cluster_id");
    const std::size_t idx_allele = required_csv_column(header_idx, "allele_id");
    const std::size_t idx_allele_len = required_csv_column(header_idx, "allele_length");
    const std::size_t idx_interval_start = required_csv_column(header_idx, "interval_start");
    const std::size_t idx_interval_end = required_csv_column(header_idx, "interval_end");
    const std::size_t idx_source_to_sink = required_csv_column(header_idx, "source_to_sink");

    std::vector<PrecomputedAssignmentRow> rows;
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
        const auto fields = split_csv_line(line);
        auto field_or_throw = [&](std::size_t idx, const std::string& name) -> const std::string& {
            if (idx >= fields.size()) {
                throw std::runtime_error(
                    "Malformed precomputed assignments CSV line " + std::to_string(line_no) +
                    ": missing field " + name);
            }
            return fields[idx];
        };

        PrecomputedAssignmentRow row;
        row.bubble_id = parse_csv_size(field_or_throw(idx_bubble, "bubble_id"), "bubble_id");
        row.source = field_or_throw(idx_source, "source");
        row.sink = field_or_throw(idx_sink, "sink");
        row.path_name = field_or_throw(idx_path, "path_name");
        row.cluster_id = parse_csv_size(field_or_throw(idx_cluster, "cluster_id"), "cluster_id");
        row.allele_id = parse_csv_size(field_or_throw(idx_allele, "allele_id"), "allele_id");
        row.allele_length = parse_csv_size(field_or_throw(idx_allele_len, "allele_length"), "allele_length");
        row.interval_start = parse_csv_size(field_or_throw(idx_interval_start, "interval_start"), "interval_start");
        row.interval_end = parse_csv_size(field_or_throw(idx_interval_end, "interval_end"), "interval_end");
        row.source_to_sink = parse_csv_bool(field_or_throw(idx_source_to_sink, "source_to_sink"), "source_to_sink");
        rows.push_back(std::move(row));
    }

    return rows;
}

} // namespace

void call_variants_from_precomputed_grouped_impl(
    const Graph& graph,
    const AlleleCallOptions& options,
    const std::vector<Bubble>& bubbles,
    const PrecomputedClustersByBubble& clusters_by_bubble,
    const PrecomputedAssignmentsByBubble& assignments_by_bubble,
    AlleleCallSummary* summary_out) {

    if (options.reference_path.empty()) {
        throw std::runtime_error("--reference-path is required for module 'call'");
    }
    if (!options.write_region_vcf && !options.write_debug_reports) {
        throw std::runtime_error("No outputs requested: enable region VCF and/or --debug output");
    }

    const auto graph_paths = path_by_name(graph);
    const auto ref_it = graph_paths.find(options.reference_path);
    if (ref_it == graph_paths.end()) {
        throw std::runtime_error("Reference path not found in graph: " + options.reference_path);
    }
    const PathRecord* reference_path_record = ref_it->second;
    const std::vector<std::size_t> reference_prefix_bp = prefix_bp(*reference_path_record, graph.nodes);
    bool complete_ref_sequence = false;
    std::string reference_sequence = spell_sequence(graph, reference_path_record->steps, complete_ref_sequence);
    if (!complete_ref_sequence) {
        reference_sequence.clear();
    }
    const ParsedReferencePath reference_meta = ::panvar::parse_reference_path_label(options.reference_path);

    std::vector<std::string> vcf_sample_names;
    std::unordered_set<std::string> seen_samples;
    seen_samples.reserve(graph.paths.size() * 2);
    for (const auto& path : graph.paths) {
        if (seen_samples.insert(path.name).second) {
            vcf_sample_names.push_back(path.name);
        }
    }

    if (options.write_debug_reports) {
        if (options.debug_out_dir.empty()) {
            throw std::runtime_error("--debug-out-dir cannot be empty when call debug output is enabled");
        }
        std::filesystem::create_directories(options.debug_out_dir);
    }

    std::ofstream debug_summary_out;
    if (options.write_debug_reports) {
        const std::string debug_summary_path = options.debug_out_dir + "/debug_summary.tsv";
        debug_summary_out.open(debug_summary_path);
        if (!debug_summary_out) {
            throw std::runtime_error("Failed to write debug summary TSV: " + debug_summary_path);
        }
        debug_summary_out
            << "bubble_id\tsource\tsink\tstatus\tassignment_rows\tunique_alleles\tclusters\t"
            << "cluster_status_rows\tclusters_with_minimap_hit\tclusters_with_events\t"
            << "bubble_vcf_rows\tdotplots_written\tcluster_debug_reports\n";
    }

    auto write_empty_bubble_debug = [&](const Bubble& bubble, const std::string& status) {
        if (!options.write_debug_reports) {
            return;
        }
        const std::string bubble_debug_dir = options.debug_out_dir + "/bubble_" + std::to_string(bubble.id);
        std::filesystem::create_directories(bubble_debug_dir);

        const std::string bubble_status_path = bubble_debug_dir + "/bubble_status.tsv";
        std::ofstream bubble_out(bubble_status_path);
        if (!bubble_out) {
            throw std::runtime_error("Failed to write bubble debug status TSV: " + bubble_status_path);
        }
        bubble_out
            << "bubble_id\tstatus\thas_reference_assignment\treference_cluster_id\treference_interval_start\t"
            << "reference_interval_end\tunique_alleles\tclusters_total\tclusters_observed\tclusters_ok\t"
            << "clusters_with_minimap_hit\tclusters_with_events\tdotplots_written\tcluster_pairwise_vcfs\t"
            << "region_vcf_rows\n";
        bubble_out
            << bubble.id << '\t'
            << status << '\t'
            << 0 << '\t'
            << 0 << '\t'
            << 0 << '\t'
            << 0 << '\t'
            << 0 << '\t'
            << 0 << '\t'
            << 0 << '\t'
            << 0 << '\t'
            << 0 << '\t'
            << 0 << '\t'
            << 0 << '\t'
            << 0 << '\t'
            << 0 << '\n';

        const std::string cluster_status_path = bubble_debug_dir + "/cluster_status.tsv";
        std::ofstream cluster_out(cluster_status_path);
        if (!cluster_out) {
            throw std::runtime_error("Failed to write cluster debug status TSV: " + cluster_status_path);
        }
        cluster_out
            << "bubble_id\tcluster_id\tstatus\tis_reference_cluster\trepresentative_allele_id\t"
            << "representative_haplotype\treference_length_bp\tcluster_length_bp\tpaf_records\t"
            << "minimap_best_ok\torientation\tbest_norm_ed\tevent_count\tdotplot_written\t"
            << "pairwise_vcf_written\n";
    };

    GeneAnnotationIndex dotplot_gene_index;
    const GeneAnnotationIndex* dotplot_gene_ptr = nullptr;
    if (options.write_debug_reports && !options.dotplot_gtf_path.empty()) {
        dotplot_gene_index = ::panvar::load_gene_annotations_from_gtf(options.dotplot_gtf_path);
        dotplot_gene_ptr = &dotplot_gene_index;
    }

    const bool use_pangene = !options.pangene_bed_path.empty();
    PangeneGeneIndex pangene_index;
    std::vector<CompiledGeneFilter> pangene_gene_filters;
    std::unordered_map<std::string, std::vector<std::size_t>> prefix_bp_by_path;
    std::ofstream pangene_copy_out;
    if (use_pangene) {
        pangene_index = load_pangene_bed(options.pangene_bed_path);
        pangene_gene_filters = compile_gene_filters(options.pangene_gene_matches);
        prefix_bp_by_path.reserve(graph.paths.size() * 2);
        for (const auto& path : graph.paths) {
            prefix_bp_by_path[path.name] = prefix_bp(path, graph.nodes);
        }
        if (!options.pangene_copy_tsv_path.empty()) {
            pangene_copy_out.open(options.pangene_copy_tsv_path);
            if (!pangene_copy_out) {
                throw std::runtime_error(
                    "Failed to write pangene copy TSV: " + options.pangene_copy_tsv_path);
            }
            pangene_copy_out
                << "bubble_id\tcluster_id\tpath_name\tis_reference\t"
                << "step_start\tstep_end\tbp_start\tbp_end\tgene_name\tcopy_count\n";
        }
    }

    std::vector<RegionVcfBubble> region_vcf_bubbles;
    region_vcf_bubbles.reserve(bubbles.size());
    AlleleCallSummary summary;
    summary.bubbles_processed = bubbles.size();
    const auto run_start = std::chrono::steady_clock::now();

    for (std::size_t bubble_idx = 0; bubble_idx < bubbles.size(); ++bubble_idx) {
        const Bubble& bubble = bubbles[bubble_idx];
        const auto bubble_start = std::chrono::steady_clock::now();
        const auto assign_it = assignments_by_bubble.find(bubble.id);
        if (assign_it == assignments_by_bubble.end() || assign_it->second.empty()) {
            write_empty_bubble_debug(bubble, "skipped:no-precomputed-assignments");
            if (debug_summary_out) {
                debug_summary_out
                    << bubble.id << '\t'
                    << bubble.source << '\t'
                    << bubble.sink << '\t'
                    << "skipped:no-precomputed-assignments" << '\t'
                    << 0 << '\t'
                    << 0 << '\t'
                    << 0 << '\t'
                    << 0 << '\t'
                    << 0 << '\t'
                    << 0 << '\t'
                    << 0 << '\t'
                    << 0 << '\t'
                    << 0 << '\n';
            }
            continue;
        }
        const auto& bubble_assign_rows = assign_it->second;

        if (options.show_progress) {
            std::cerr
                << "[call] bubble " << (bubble_idx + 1) << "/" << bubbles.size()
                << " id=" << bubble.id
                << " precomputed_rows=" << bubble_assign_rows.size()
                << "\n";
        }

        std::unordered_map<std::size_t, std::size_t> unique_idx_by_allele_id;
        unique_idx_by_allele_id.reserve(bubble_assign_rows.size() * 2);
        std::vector<UniqueAllele> unique_alleles;
        unique_alleles.reserve(bubble_assign_rows.size());
        std::vector<PathAssignment> assignments;
        assignments.reserve(bubble_assign_rows.size());
        std::vector<std::size_t> assignment_cluster_ids;
        assignment_cluster_ids.reserve(bubble_assign_rows.size());

        for (const auto& row : bubble_assign_rows) {
            const auto path_it = graph_paths.find(row.path_name);
            if (path_it == graph_paths.end()) {
                continue;
            }
            const PathRecord* path = path_it->second;
            const auto known = unique_idx_by_allele_id.find(row.allele_id);

            std::size_t unique_idx = 0;
            if (known == unique_idx_by_allele_id.end()) {
                UniqueAllele allele;
                allele.allele_id = row.allele_id;
                allele.path_support = 0;

                CandidateInterval interval;
                interval.left = row.interval_start;
                interval.right = row.interval_end;
                interval.source_to_sink = row.source_to_sink;
                interval.inside_count = (interval.right > interval.left) ? (interval.right - interval.left - 1) : 0;

                allele.steps = canonical_steps_for_bubble(*path, bubble, interval);
                if (allele.steps.size() < 2 ||
                    allele.steps.front().node_id != bubble.source ||
                    allele.steps.back().node_id != bubble.sink) {
                    if (interval.left <= interval.right && interval.right < path->steps.size()) {
                        allele.steps.assign(
                            path->steps.begin() + static_cast<std::ptrdiff_t>(interval.left),
                            path->steps.begin() + static_cast<std::ptrdiff_t>(interval.right + 1));
                        if (!interval.source_to_sink) {
                            std::reverse(allele.steps.begin(), allele.steps.end());
                            for (auto& step : allele.steps) {
                                step.reverse = !step.reverse;
                            }
                        }
                    }
                }

                allele.signature = build_walk_signature(allele.steps);
                bool has_sequence = false;
                allele.sequence = spell_sequence(graph, allele.steps, has_sequence);
                allele.compare_steps = build_walk_tokens(allele.steps);
                if (has_sequence) {
                    allele.compare_token = allele.sequence;
                    allele.uses_sequence_similarity = true;
                    allele.sequence_length = allele.sequence.size();
                } else {
                    allele.compare_token = allele.signature;
                    allele.uses_sequence_similarity = false;
                    allele.sequence_length = row.allele_length > 0 ? row.allele_length : allele.steps.size();
                }

                unique_alleles.push_back(std::move(allele));
                unique_idx = unique_alleles.size() - 1;
                unique_idx_by_allele_id[row.allele_id] = unique_idx;
            } else {
                unique_idx = known->second;
            }

            unique_alleles[unique_idx].path_support += 1;

            PathAssignment assignment;
            assignment.unique_idx = unique_idx;
            assignment.path_name = row.path_name;
            assignment.interval_start = row.interval_start;
            assignment.interval_end = row.interval_end;
            assignment.source_to_sink = row.source_to_sink;
            assignments.push_back(std::move(assignment));
            assignment_cluster_ids.push_back(row.cluster_id);
        }

        if (assignments.empty() || unique_alleles.empty()) {
            write_empty_bubble_debug(bubble, "skipped:no-usable-assignments");
            if (debug_summary_out) {
                debug_summary_out
                    << bubble.id << '\t'
                    << bubble.source << '\t'
                    << bubble.sink << '\t'
                    << "skipped:no-usable-assignments" << '\t'
                    << bubble_assign_rows.size() << '\t'
                    << unique_alleles.size() << '\t'
                    << 0 << '\t'
                    << 0 << '\t'
                    << 0 << '\t'
                    << 0 << '\t'
                    << 0 << '\t'
                    << 0 << '\t'
                    << 0 << '\n';
            }
            continue;
        }

        summary.bubbles_with_assignments += 1;
        summary.unique_alleles += unique_alleles.size();
        summary.assignments += assignments.size();

        std::vector<std::size_t> cluster_of_unique(unique_alleles.size(), 0);
        std::unordered_map<std::size_t, std::vector<std::size_t>> member_by_cluster_from_assignments;
        std::unordered_map<std::size_t, std::size_t> support_by_cluster_from_assignments;
        member_by_cluster_from_assignments.reserve(unique_alleles.size() * 2);
        support_by_cluster_from_assignments.reserve(unique_alleles.size() * 2);

        for (std::size_t i = 0; i < assignments.size(); ++i) {
            const std::size_t cluster_id = assignment_cluster_ids[i];
            const std::size_t unique_idx = assignments[i].unique_idx;
            cluster_of_unique[unique_idx] = cluster_id;
            member_by_cluster_from_assignments[cluster_id].push_back(unique_idx);
            support_by_cluster_from_assignments[cluster_id] += 1;
        }
        for (auto& [cluster_id, members] : member_by_cluster_from_assignments) {
            (void)cluster_id;
            std::sort(members.begin(), members.end());
            members.erase(std::unique(members.begin(), members.end()), members.end());
        }

        std::vector<Cluster> clusters;
        const auto cluster_it = clusters_by_bubble.find(bubble.id);
        if (cluster_it != clusters_by_bubble.end()) {
            clusters.reserve(cluster_it->second.size());
            for (const auto& row : cluster_it->second) {
                Cluster cluster;
                cluster.cluster_id = row.cluster_id;
                cluster.total_path_support = row.total_path_support;

                if (!row.member_alleles.empty()) {
                    for (const std::size_t allele_id : row.member_alleles) {
                        const auto uid_it = unique_idx_by_allele_id.find(allele_id);
                        if (uid_it != unique_idx_by_allele_id.end()) {
                            cluster.member_unique_idxs.push_back(uid_it->second);
                        }
                    }
                    std::sort(cluster.member_unique_idxs.begin(), cluster.member_unique_idxs.end());
                    cluster.member_unique_idxs.erase(
                        std::unique(cluster.member_unique_idxs.begin(), cluster.member_unique_idxs.end()),
                        cluster.member_unique_idxs.end());
                } else {
                    const auto members_it = member_by_cluster_from_assignments.find(cluster.cluster_id);
                    if (members_it != member_by_cluster_from_assignments.end()) {
                        cluster.member_unique_idxs = members_it->second;
                    }
                }
                if (cluster.member_unique_idxs.empty()) {
                    continue;
                }

                const auto rep_it = unique_idx_by_allele_id.find(row.representative_allele_id);
                cluster.representative_idx = (rep_it == unique_idx_by_allele_id.end())
                                                 ? cluster.member_unique_idxs.front()
                                                 : rep_it->second;
                if (cluster.total_path_support == 0) {
                    const auto support_it = support_by_cluster_from_assignments.find(cluster.cluster_id);
                    if (support_it != support_by_cluster_from_assignments.end()) {
                        cluster.total_path_support = support_it->second;
                    } else {
                        for (const std::size_t member : cluster.member_unique_idxs) {
                            cluster.total_path_support += unique_alleles[member].path_support;
                        }
                    }
                }
                clusters.push_back(std::move(cluster));
            }
        }

        if (clusters.empty()) {
            clusters.reserve(member_by_cluster_from_assignments.size());
            for (const auto& [cluster_id, members] : member_by_cluster_from_assignments) {
                if (members.empty()) {
                    continue;
                }
                Cluster cluster;
                cluster.cluster_id = cluster_id;
                cluster.member_unique_idxs = members;
                cluster.representative_idx = members.front();
                cluster.total_path_support = support_by_cluster_from_assignments[cluster_id];
                clusters.push_back(std::move(cluster));
            }
        }

        std::unordered_set<std::size_t> existing_cluster_ids;
        existing_cluster_ids.reserve(clusters.size() * 2);
        for (auto& cluster : clusters) {
            existing_cluster_ids.insert(cluster.cluster_id);
            for (const auto unique_idx : cluster.member_unique_idxs) {
                cluster_of_unique[unique_idx] = cluster.cluster_id;
            }
            if (cluster.total_path_support == 0) {
                for (const auto unique_idx : cluster.member_unique_idxs) {
                    cluster.total_path_support += unique_alleles[unique_idx].path_support;
                }
            }
        }

        for (std::size_t unique_idx = 0; unique_idx < cluster_of_unique.size(); ++unique_idx) {
            const std::size_t cluster_id = cluster_of_unique[unique_idx];
            if (cluster_id == 0 || existing_cluster_ids.count(cluster_id) > 0) {
                continue;
            }
            Cluster cluster;
            cluster.cluster_id = cluster_id;
            cluster.member_unique_idxs = {unique_idx};
            cluster.representative_idx = unique_idx;
            cluster.total_path_support = unique_alleles[unique_idx].path_support;
            clusters.push_back(std::move(cluster));
            existing_cluster_ids.insert(cluster_id);
        }

        std::sort(clusters.begin(), clusters.end(), [](const Cluster& lhs, const Cluster& rhs) {
            return lhs.cluster_id < rhs.cluster_id;
        });

        summary.clusters += clusters.size();

        VariantBubbleReport report = ::panvar::write_variant_reports_for_bubble(
            bubble,
            unique_alleles,
            clusters,
            cluster_of_unique,
            assignments,
            options.reference_path,
            reference_sequence,
            options.min_sv_bp,
            options.write_debug_reports,
            options.debug_out_dir,
            dotplot_gene_ptr,
            options.dotplot_gene_matches,
            options.minimap_preset,
            options.minimap_best_n,
            options.minimap_emit_secondary,
            options.split_ins_use_geometric_svlen,
            options.classify_ins);

        summary.debug_reports_written += report.debug_reports_written;
        summary.dotplots_written += report.dotplots_written;

        if (debug_summary_out) {
            debug_summary_out
                << bubble.id << '\t'
                << bubble.source << '\t'
                << bubble.sink << '\t'
                << report.status << '\t'
                << assignments.size() << '\t'
                << unique_alleles.size() << '\t'
                << clusters.size() << '\t'
                << report.cluster_statuses.size() << '\t'
                << report.clusters_with_minimap_hit << '\t'
                << report.clusters_with_events << '\t'
                << report.vcf_rows.size() << '\t'
                << report.dotplots_written << '\t'
                << report.debug_reports_written << '\n';
        }

        if (use_pangene && report.has_reference_assignment && !report.vcf_rows.empty()) {
            struct SelectedAssignmentInterval {
                std::string path_name;
                std::size_t step_start = 0;
                std::size_t step_end = 0;
                bool source_to_sink = true;
            };

            auto cluster_path_key = [](std::size_t cluster_id, const std::string& path_name) {
                return std::to_string(cluster_id) + '\t' + path_name;
            };

            std::unordered_map<std::string, SelectedAssignmentInterval> best_by_cluster_path;
            std::unordered_map<std::size_t, SelectedAssignmentInterval> best_by_cluster;
            best_by_cluster_path.reserve(assignments.size() * 2);
            best_by_cluster.reserve(clusters.size() * 2);

            for (std::size_t i = 0; i < assignments.size(); ++i) {
                const PathAssignment& assignment = assignments[i];
                if (i >= assignment_cluster_ids.size()) {
                    continue;
                }
                const std::size_t cluster_id = assignment_cluster_ids[i];
                if (cluster_id == 0) {
                    continue;
                }
                SelectedAssignmentInterval candidate;
                candidate.path_name = assignment.path_name;
                candidate.step_start = assignment.interval_start;
                candidate.step_end = assignment.interval_end;
                candidate.source_to_sink = assignment.source_to_sink;

                const std::string key = cluster_path_key(cluster_id, assignment.path_name);
                auto by_path_it = best_by_cluster_path.find(key);
                if (by_path_it == best_by_cluster_path.end() ||
                    (!by_path_it->second.source_to_sink && candidate.source_to_sink)) {
                    best_by_cluster_path[key] = candidate;
                }

                auto by_cluster_it = best_by_cluster.find(cluster_id);
                if (by_cluster_it == best_by_cluster.end() ||
                    (!by_cluster_it->second.source_to_sink && candidate.source_to_sink) ||
                    (by_cluster_it->second.source_to_sink == candidate.source_to_sink &&
                     candidate.path_name < by_cluster_it->second.path_name)) {
                    best_by_cluster[cluster_id] = std::move(candidate);
                }
            }

            auto counts_for_selected = [&](const SelectedAssignmentInterval& selection) -> GeneCopyCounts {
                auto pref_it = prefix_bp_by_path.find(selection.path_name);
                if (pref_it == prefix_bp_by_path.end()) {
                    return {};
                }
                const auto& pref = pref_it->second;
                const std::size_t lo = std::min(selection.step_start, selection.step_end);
                const std::size_t hi = std::max(selection.step_start, selection.step_end);
                if (lo >= pref.size() || hi + 1 >= pref.size()) {
                    return {};
                }
                const std::size_t bp_start = pref[lo];
                const std::size_t bp_end = pref[hi + 1];
                return count_pangene_copies_for_interval(
                    pangene_index,
                    selection.path_name,
                    bp_start,
                    bp_end,
                    pangene_gene_filters);
            };

            std::unordered_map<std::string, GeneCopyCounts> counts_by_cluster_path;
            counts_by_cluster_path.reserve(best_by_cluster_path.size() * 2);
            for (const auto& kv : best_by_cluster_path) {
                counts_by_cluster_path[kv.first] = counts_for_selected(kv.second);
            }
            std::unordered_map<std::size_t, GeneCopyCounts> counts_by_cluster;
            counts_by_cluster.reserve(best_by_cluster.size() * 2);
            for (const auto& kv : best_by_cluster) {
                counts_by_cluster[kv.first] = counts_for_selected(kv.second);
            }

            const std::string ref_key = cluster_path_key(report.reference_cluster_id, options.reference_path);
            GeneCopyCounts reference_counts;
            const auto ref_by_path_it = counts_by_cluster_path.find(ref_key);
            if (ref_by_path_it != counts_by_cluster_path.end()) {
                reference_counts = ref_by_path_it->second;
            } else {
                const auto ref_by_cluster_it = counts_by_cluster.find(report.reference_cluster_id);
                if (ref_by_cluster_it != counts_by_cluster.end()) {
                    reference_counts = ref_by_cluster_it->second;
                }
            }

            if (pangene_copy_out) {
                for (const auto& kv : best_by_cluster) {
                    const std::size_t cluster_id = kv.first;
                    const auto& selection = kv.second;
                    auto pref_it = prefix_bp_by_path.find(selection.path_name);
                    if (pref_it == prefix_bp_by_path.end()) {
                        continue;
                    }
                    const auto& pref = pref_it->second;
                    const std::size_t lo = std::min(selection.step_start, selection.step_end);
                    const std::size_t hi = std::max(selection.step_start, selection.step_end);
                    if (lo >= pref.size() || hi + 1 >= pref.size()) {
                        continue;
                    }
                    const std::size_t bp_start = pref[lo];
                    const std::size_t bp_end = pref[hi + 1];
                    const auto counts_it = counts_by_cluster.find(cluster_id);
                    if (counts_it == counts_by_cluster.end() || counts_it->second.empty()) {
                        continue;
                    }
                    for (const auto& gene_kv : counts_it->second) {
                        pangene_copy_out
                            << bubble.id << '\t'
                            << cluster_id << '\t'
                            << selection.path_name << '\t'
                            << (cluster_id == report.reference_cluster_id ? 1 : 0) << '\t'
                            << lo << '\t'
                            << hi << '\t'
                            << bp_start << '\t'
                            << bp_end << '\t'
                            << gene_kv.first << '\t'
                            << gene_kv.second << '\n';
                    }
                }
            }

            for (auto& row : report.vcf_rows) {
                GeneCopyCounts alt_counts;
                const std::string key =
                    cluster_path_key(row.cluster_id, row.representative_haplotype);
                const auto by_path_it = counts_by_cluster_path.find(key);
                if (by_path_it != counts_by_cluster_path.end()) {
                    alt_counts = by_path_it->second;
                } else {
                    const auto by_cluster_it = counts_by_cluster.find(row.cluster_id);
                    if (by_cluster_it != counts_by_cluster.end()) {
                        alt_counts = by_cluster_it->second;
                    }
                }

                const PangeneDeltaSummary delta =
                    summarize_gene_copy_delta(reference_counts, alt_counts);
                row.pangene_cn_delta = delta.cn_delta;
                row.pangene_gain_genes = delta.gain_genes;
                row.pangene_loss_genes = delta.loss_genes;
                row.pangene_gain_copies = delta.gain_copies;
                row.pangene_loss_copies = delta.loss_copies;
                if (options.pangene_tune_ins &&
                    row.event_type == "INS" &&
                    row.pangene_gain_copies > 0 &&
                    (row.event_subtype == "." || row.event_subtype == "NOVEL")) {
                    row.event_subtype = "DUP_PANGENE";
                }
            }
        }

        if (options.write_region_vcf &&
            report.has_reference_assignment && !report.vcf_rows.empty()) {
            RegionVcfBubble item;
            item.bubble_id = bubble.id;
            item.source = bubble.source;
            item.sink = bubble.sink;
            item.reference_cluster_id = report.reference_cluster_id;
            item.reference_interval_start = report.reference_interval_start;
            item.reference_interval_end = report.reference_interval_end;
            item.cluster_by_path = std::move(report.cluster_by_path);
            item.allele_length_by_path = std::move(report.allele_length_by_path);
            item.variant_rows = std::move(report.vcf_rows);
            region_vcf_bubbles.push_back(std::move(item));
        }

        if (options.show_progress) {
            std::cerr
                << "[call] bubble " << bubble.id
                << " done: unique=" << unique_alleles.size()
                << ", clusters=" << clusters.size()
                << ", paths=" << assignments.size()
                << ", elapsed=" << std::fixed << std::setprecision(2)
                << elapsed_seconds(bubble_start) << "s"
                << ", total=" << elapsed_seconds(run_start) << "s\n";
        }
    }

    if (options.write_region_vcf) {
        if (options.region_vcf_path.empty()) {
            throw std::runtime_error("--vcf-out cannot be empty when region-level VCF export is enabled");
        }
        if (reference_prefix_bp.empty()) {
            throw std::runtime_error("Reference path prefix lengths are unavailable for region-level VCF export");
        }
        const std::size_t region_start_1based =
            reference_meta.has_interval ? reference_meta.region_start_1based : 1;
        ::panvar::write_region_level_vcf(
            options.region_vcf_path,
            options.reference_path,
            reference_meta,
            reference_sequence,
            reference_prefix_bp,
            region_vcf_bubbles,
            vcf_sample_names,
            region_start_1based,
            options.vcf_merge_window_bp,
            options.vcf_merge_mode,
            options.vcf_merge_lenient_window_bp,
            options.vcf_merge_lenient_min_ref_jaccard,
            options.vcf_merge_min_seq_similarity,
            options.vcf_merge_max_seq_edit_fraction,
            &summary.region_vcf_records);
    }

    if (summary_out != nullptr) {
        *summary_out = summary;
    }
}

void call_variants_from_precomputed_alleles(
    const Graph& graph,
    const AlleleCallOptions& options,
    const std::string& clusters_csv_in_path,
    const std::string& assignments_csv_in_path,
    AlleleCallSummary* summary_out) {

    if (options.bubbles_csv_in.empty()) {
        throw std::runtime_error(
            "Module 'call' with precomputed alleles requires module-1 bubbles "
            "(--bubble-prefix-in <prefix> or --bubbles-csv-in <path>).");
    }
    if (options.reference_path.empty()) {
        throw std::runtime_error("--reference-path is required for module 'call'");
    }
    if (!options.write_region_vcf && !options.write_debug_reports) {
        throw std::runtime_error("No outputs requested: enable region VCF and/or --debug output");
    }
    validate_call_merge_options(options);

    const std::vector<Bubble> bubbles = read_bubbles_csv(options.bubbles_csv_in);
    const std::vector<PrecomputedClusterRow> pre_clusters = read_precomputed_clusters_csv(clusters_csv_in_path);
    const std::vector<PrecomputedAssignmentRow> pre_assignments =
        read_precomputed_assignments_csv(assignments_csv_in_path);

    std::unordered_map<std::size_t, const Bubble*> bubble_by_id;
    bubble_by_id.reserve(bubbles.size() * 2);
    std::unordered_map<std::string, std::size_t> unique_bubble_id_by_endpoint;
    unique_bubble_id_by_endpoint.reserve(bubbles.size() * 2);
    std::unordered_set<std::string> ambiguous_endpoints;
    ambiguous_endpoints.reserve(8);
    for (const auto& bubble : bubbles) {
        bubble_by_id[bubble.id] = &bubble;
        const std::string key = endpoint_key(bubble.source, bubble.sink);
        const auto it = unique_bubble_id_by_endpoint.find(key);
        if (it == unique_bubble_id_by_endpoint.end()) {
            unique_bubble_id_by_endpoint.emplace(key, bubble.id);
        } else if (it->second != bubble.id) {
            ambiguous_endpoints.insert(key);
        }
    }
    for (const auto& key : ambiguous_endpoints) {
        unique_bubble_id_by_endpoint.erase(key);
    }

    auto remap_bubble_id = [&](
                              std::size_t input_bubble_id,
                              const std::string& source,
                              const std::string& sink,
                              const char* row_type) -> std::size_t {
        const std::string key = endpoint_key(source, sink);
        const auto endpoint_it = unique_bubble_id_by_endpoint.find(key);
        const auto id_it = bubble_by_id.find(input_bubble_id);

        if (id_it != bubble_by_id.end()) {
            const Bubble* by_id = id_it->second;
            if (by_id->source == source && by_id->sink == sink) {
                return input_bubble_id;
            }
        }

        if (endpoint_it != unique_bubble_id_by_endpoint.end()) {
            return endpoint_it->second;
        }

        std::ostringstream oss;
        oss
            << "Precomputed " << row_type
            << " row could not be mapped to bubbles CSV (bubble_id=" << input_bubble_id
            << ", source=" << source
            << ", sink=" << sink << "). "
            << "Make sure module-1 input (--bubble-prefix-in/--bubbles-csv-in) "
            << "matches module-2 CSVs.";
        throw std::runtime_error(oss.str());
    };

    std::size_t remapped_cluster_rows = 0;
    std::size_t remapped_assignment_rows = 0;
    PrecomputedClustersByBubble clusters_by_bubble;
    clusters_by_bubble.reserve(bubbles.size() * 2);
    for (const auto& input_row : pre_clusters) {
        PrecomputedClusterRow row = input_row;
        const std::size_t mapped_bubble_id =
            remap_bubble_id(row.bubble_id, row.source, row.sink, "cluster");
        if (mapped_bubble_id != row.bubble_id) {
            ++remapped_cluster_rows;
        }
        row.bubble_id = mapped_bubble_id;
        const Bubble* canonical = bubble_by_id.at(mapped_bubble_id);
        row.source = canonical->source;
        row.sink = canonical->sink;
        clusters_by_bubble[row.bubble_id].push_back(std::move(row));
    }

    PrecomputedAssignmentsByBubble assignments_by_bubble;
    assignments_by_bubble.reserve(bubbles.size() * 2);
    for (const auto& input_row : pre_assignments) {
        PrecomputedAssignmentRow row = input_row;
        const std::size_t mapped_bubble_id =
            remap_bubble_id(row.bubble_id, row.source, row.sink, "assignment");
        if (mapped_bubble_id != row.bubble_id) {
            ++remapped_assignment_rows;
        }
        row.bubble_id = mapped_bubble_id;
        const Bubble* canonical = bubble_by_id.at(mapped_bubble_id);
        row.source = canonical->source;
        row.sink = canonical->sink;
        assignments_by_bubble[row.bubble_id].push_back(std::move(row));
    }

    if (options.show_progress &&
        (remapped_cluster_rows > 0 || remapped_assignment_rows > 0)) {
        std::cerr
            << "[call] remapped precomputed rows to current bubbles by source/sink: "
            << "clusters=" << remapped_cluster_rows
            << ", assignments=" << remapped_assignment_rows
            << '\n';
    }

    call_variants_from_precomputed_grouped_impl(
        graph,
        options,
        bubbles,
        clusters_by_bubble,
        assignments_by_bubble,
        summary_out);
}

} // namespace panvar
