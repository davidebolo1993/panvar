#include "panvar/call.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

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
    const std::size_t idx_nesting = required_csv_column(header_idx, "nesting_level");
    const std::size_t idx_cluster = required_csv_column(header_idx, "cluster_id");
    const std::size_t idx_rep_allele = required_csv_column(header_idx, "representative_allele_id");
    const std::size_t idx_rep_len = required_csv_column(header_idx, "representative_length");
    const std::size_t idx_rep_mode = required_csv_column(header_idx, "representative_mode");
    const std::size_t idx_member_count = required_csv_column(header_idx, "member_allele_count");
    const std::size_t idx_support = required_csv_column(header_idx, "total_path_support");
    const std::size_t idx_members = required_csv_column(header_idx, "member_alleles");
    const std::size_t idx_signature = required_csv_column(header_idx, "representative_signature");

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
        row.nesting_level = parse_csv_size(field_or_throw(idx_nesting, "nesting_level"), "nesting_level");
        row.cluster_id = parse_csv_size(field_or_throw(idx_cluster, "cluster_id"), "cluster_id");
        row.representative_allele_id =
            parse_csv_size(field_or_throw(idx_rep_allele, "representative_allele_id"), "representative_allele_id");
        row.representative_length =
            parse_csv_size(field_or_throw(idx_rep_len, "representative_length"), "representative_length");
        row.representative_mode = field_or_throw(idx_rep_mode, "representative_mode");
        row.member_allele_count =
            parse_csv_size(field_or_throw(idx_member_count, "member_allele_count"), "member_allele_count");
        row.total_path_support =
            parse_csv_size(field_or_throw(idx_support, "total_path_support"), "total_path_support");
        row.member_alleles = split_semicolon_size_t(field_or_throw(idx_members, "member_alleles"));
        row.representative_signature = field_or_throw(idx_signature, "representative_signature");
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
    const std::size_t idx_path_type = required_csv_column(header_idx, "path_type");
    const std::size_t idx_cluster = required_csv_column(header_idx, "cluster_id");
    const std::size_t idx_allele = required_csv_column(header_idx, "allele_id");
    const std::size_t idx_allele_len = required_csv_column(header_idx, "allele_length");
    const std::size_t idx_interval_start = required_csv_column(header_idx, "interval_start");
    const std::size_t idx_interval_end = required_csv_column(header_idx, "interval_end");
    const std::size_t idx_interval_steps = required_csv_column(header_idx, "interval_step_count");
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
        const std::string type_field = field_or_throw(idx_path_type, "path_type");
        row.path_type = type_field.empty() ? 'P' : type_field[0];
        row.cluster_id = parse_csv_size(field_or_throw(idx_cluster, "cluster_id"), "cluster_id");
        row.allele_id = parse_csv_size(field_or_throw(idx_allele, "allele_id"), "allele_id");
        row.allele_length = parse_csv_size(field_or_throw(idx_allele_len, "allele_length"), "allele_length");
        row.interval_start = parse_csv_size(field_or_throw(idx_interval_start, "interval_start"), "interval_start");
        row.interval_end = parse_csv_size(field_or_throw(idx_interval_end, "interval_end"), "interval_end");
        row.interval_step_count =
            parse_csv_size(field_or_throw(idx_interval_steps, "interval_step_count"), "interval_step_count");
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

    GeneAnnotationIndex dotplot_gene_index;
    const GeneAnnotationIndex* dotplot_gene_ptr = nullptr;
    if (options.write_debug_reports && !options.dotplot_gtf_path.empty()) {
        dotplot_gene_index = ::panvar::load_gene_annotations_from_gtf(options.dotplot_gtf_path);
        dotplot_gene_ptr = &dotplot_gene_index;
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
            assignment.path_type = row.path_type;
            assignment.interval_start = row.interval_start;
            assignment.interval_end = row.interval_end;
            assignment.interval_step_count = row.interval_step_count;
            assignment.source_to_sink = row.source_to_sink;
            assignments.push_back(std::move(assignment));
            assignment_cluster_ids.push_back(row.cluster_id);
        }

        if (assignments.empty() || unique_alleles.empty()) {
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
            "Module 'call' with precomputed alleles still requires --bubbles-csv-in from module 1.");
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
            << "Make sure --bubbles-csv-in matches the module-2 CSVs.";
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
