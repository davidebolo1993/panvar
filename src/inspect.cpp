#include "panvar/inspect.hpp"

#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/output.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <zlib.h>

namespace panvar {
namespace {

void print_inspect_help() {
    std::cout
        << "Usage:\n"
        << "  panvar inspect -i <graph.gfa> --bubbles-csv <bubble.bubbles.csv> [--bubble-id <N>] [options]\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Input GFA file (required)\n"
        << "      --bubbles-csv <path>         Module-1 bubbles CSV (required if no prefix)\n"
        << "      --bubble-prefix-in <prefix>  Module-1 output prefix from 'panvar bubble'\n"
        << "                                   (auto-uses <prefix>.bubbles.csv)\n"
        << "      --bubble-id <N>              Bubble ID to inspect (default: inspect all bubbles)\n"
        << "  -o, --out-prefix <prefix>        Output prefix (default: inspect)\n"
        << "      --fasta-out <path>           Explicit FASTA.GZ path (requires --bubble-id)\n"
        << "      --table-out <path>           Explicit node-count TSV path (requires --bubble-id)\n"
        << "      --edge-table-out <path>      Explicit edge-count TSV path (requires --bubble-id)\n"
        << "  -h, --help                       Show this help\n";
}

struct InspectNodeCount {
    std::size_t total = 0;
    std::size_t forward = 0;
    std::size_t reverse = 0;
};

std::string tsv_sanitize(std::string value) {
    for (char& c : value) {
        if (c == '\t' || c == '\n' || c == '\r') {
            c = ' ';
        }
    }
    return value;
}

void gz_write_string(gzFile file, const std::string& text, const std::string& output_path) {
    if (text.empty()) {
        return;
    }
    const int written = gzwrite(file, text.data(), static_cast<unsigned int>(text.size()));
    if (written != static_cast<int>(text.size())) {
        throw std::runtime_error("Failed to write compressed FASTA: " + output_path);
    }
}

void gz_write_fasta_record(
    gzFile file,
    const std::string& output_path,
    const std::string& name,
    const std::string& sequence) {

    gz_write_string(file, ">" + name + "\n", output_path);
    constexpr std::size_t kLineWidth = 80;
    for (std::size_t i = 0; i < sequence.size(); i += kLineWidth) {
        gz_write_string(file, sequence.substr(i, kLineWidth) + "\n", output_path);
    }
}

struct InspectBubbleResult {
    std::size_t paths_written = 0;
    std::string fasta_out_path;
    std::string table_out_path;
    std::string edge_table_out_path;
};

// Orientation-aware edge key for a step pair, e.g. "12+>13-". Adjacency-aware so a
// tandem self-loop (the same edge traversed repeatedly) shows up as a high count.
std::string edge_key(const PathStep& a, const PathStep& b) {
    return a.node_id + (a.reverse ? '-' : '+') + ">" + b.node_id + (b.reverse ? '-' : '+');
}

// One path's traversal of a bubble: node counts and edge counts.
struct InspectPathRow {
    std::string name;
    std::size_t sequence_length = 0;
    std::unordered_map<std::string, InspectNodeCount> node_counts;
    std::unordered_map<std::string, std::size_t> edge_counts;
};

InspectBubbleResult write_inspect_outputs_for_bubble(
    const Graph& graph,
    const Bubble& bubble,
    const std::string& fasta_out_path,
    const std::string& table_out_path,
    const std::string& edge_table_out_path) {

    cli::ensure_parent_dir_for_file(fasta_out_path);
    cli::ensure_parent_dir_for_file(table_out_path);
    cli::ensure_parent_dir_for_file(edge_table_out_path);

    gzFile fasta = gzopen(fasta_out_path.c_str(), "wb");
    if (fasta == nullptr) {
        throw std::runtime_error("Failed to write compressed FASTA: " + fasta_out_path);
    }

    std::unordered_set<std::string> inside_nodes;
    inside_nodes.reserve(bubble.inside.size() * 2);
    inside_nodes.insert(bubble.inside.begin(), bubble.inside.end());

    // Pass 1: collect per-path node + edge counts and stream the FASTA. Edge columns
    // are not known up front, so the matrices are written after the union is built.
    std::vector<InspectPathRow> rows;
    std::vector<std::string> edge_order;          // edge columns, first-seen order
    std::unordered_set<std::string> edge_seen;
    try {
        for (const auto& path : graph.paths) {
            const BubblePathIndex index = build_bubble_path_index(path);
            const auto interval = find_best_bubble_path_interval(index, bubble);
            if (!interval.has_value()) {
                continue;
            }

            const std::vector<PathStep> steps = canonical_bubble_path_steps(path, bubble, *interval);
            if (steps.empty()) {
                continue;
            }

            const std::string sequence = spell_path_steps_sequence(graph, steps);
            InspectPathRow row;
            row.name = tsv_sanitize(path.name);
            row.sequence_length = sequence.size();
            row.node_counts.reserve(bubble.inside.size());
            for (const auto& step : steps) {
                if (inside_nodes.find(step.node_id) == inside_nodes.end()) {
                    continue;
                }
                auto& c = row.node_counts[step.node_id];
                c.total += 1;
                if (step.reverse) {
                    c.reverse += 1;
                } else {
                    c.forward += 1;
                }
            }
            for (std::size_t i = 1; i < steps.size(); ++i) {
                const std::string key = edge_key(steps[i - 1], steps[i]);
                if (edge_seen.insert(key).second) {
                    edge_order.push_back(key);
                }
                row.edge_counts[key] += 1;
            }

            const std::string fasta_name =
                row.name +
                " bubble=" + std::to_string(bubble.id) +
                " source=" + bubble.source +
                " sink=" + bubble.sink +
                " length_bp=" + std::to_string(row.sequence_length) +
                " source_to_sink=" + (interval->source_to_sink ? "1" : "0") +
                " interval=" + std::to_string(interval->left) + "-" + std::to_string(interval->right);
            gz_write_fasta_record(fasta, fasta_out_path, fasta_name, sequence);

            rows.push_back(std::move(row));
        }
    } catch (...) {
        gzclose(fasta);
        throw;
    }
    if (gzclose(fasta) != Z_OK) {
        throw std::runtime_error("Failed to close compressed FASTA: " + fasta_out_path);
    }

    // Stable, deterministic edge column order.
    std::sort(edge_order.begin(), edge_order.end());

    std::ofstream table_out(table_out_path);
    if (!table_out) {
        throw std::runtime_error("Failed to write inspect table: " + table_out_path);
    }
    std::ofstream edge_out(edge_table_out_path);
    if (!edge_out) {
        throw std::runtime_error("Failed to write inspect edge table: " + edge_table_out_path);
    }

    table_out << "path_name\tpath_length_bp";
    for (const auto& node_id : bubble.inside) {
        table_out << "\tnode." << tsv_sanitize(node_id);
    }
    table_out << "\n";
    edge_out << "path_name\tpath_length_bp";
    for (const auto& key : edge_order) {
        edge_out << "\tedge." << key;
    }
    edge_out << "\n";

    InspectBubbleResult result;
    result.fasta_out_path = fasta_out_path;
    result.table_out_path = table_out_path;
    result.edge_table_out_path = edge_table_out_path;
    for (const InspectPathRow& row : rows) {
        table_out << row.name << '\t' << row.sequence_length;
        for (const auto& node_id : bubble.inside) {
            const auto it = row.node_counts.find(node_id);
            if (it == row.node_counts.end()) {
                table_out << "\t0:0:0";
            } else {
                table_out << '\t' << it->second.total << ':' << it->second.forward << ':'
                          << it->second.reverse;
            }
        }
        table_out << '\n';

        edge_out << row.name << '\t' << row.sequence_length;
        for (const auto& key : edge_order) {
            const auto it = row.edge_counts.find(key);
            edge_out << '\t' << (it == row.edge_counts.end() ? 0 : it->second);
        }
        edge_out << '\n';
        ++result.paths_written;
    }
    return result;
}

} // namespace

int run_inspect_command(const std::vector<std::string>& args) {
    std::string gfa_path;
    std::string bubbles_csv_path;
    std::string bubble_prefix_in;
    std::string out_prefix = "inspect";
    std::string fasta_out_path;
    std::string table_out_path;
    std::string edge_table_out_path;
    std::size_t bubble_id = 0;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) {
                throw std::runtime_error("Missing value after " + flag);
            }
            return args[++i];
        };

        if (arg == "-h" || arg == "--help") {
            print_inspect_help();
            return 0;
        }
        if (arg == "-i" || arg == "--gfa") {
            gfa_path = require_value(arg);
            continue;
        }
        if (arg == "--bubbles-csv") {
            bubbles_csv_path = require_value(arg);
            continue;
        }
        if (arg == "--bubble-prefix-in") {
            bubble_prefix_in = require_value(arg);
            continue;
        }
        if (arg == "--bubble-id") {
            bubble_id = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "-o" || arg == "--out-prefix") {
            out_prefix = require_value(arg);
            continue;
        }
        if (arg == "--fasta-out") {
            fasta_out_path = require_value(arg);
            continue;
        }
        if (arg == "--table-out") {
            table_out_path = require_value(arg);
            continue;
        }
        if (arg == "--edge-table-out") {
            edge_table_out_path = require_value(arg);
            continue;
        }
        throw std::runtime_error("Unknown option: " + arg);
    }

    if (gfa_path.empty()) {
        throw std::runtime_error("Missing required input: --gfa <path>");
    }
    if (!bubble_prefix_in.empty()) {
        const std::string derived = bubble_prefix_in + ".bubbles.csv";
        if (bubbles_csv_path.empty()) {
            bubbles_csv_path = derived;
        } else if (bubbles_csv_path != derived) {
            throw std::runtime_error(
                "Conflicting bubble inputs: --bubble-prefix-in resolves to '" +
                derived + "' but --bubbles-csv is '" + bubbles_csv_path + "'");
        }
    }
    if (bubbles_csv_path.empty()) {
        throw std::runtime_error("Missing required input: --bubbles-csv <path> or --bubble-prefix-in <prefix>");
    }
    if (bubble_id == 0 && (!fasta_out_path.empty() || !table_out_path.empty() || !edge_table_out_path.empty())) {
        throw std::runtime_error("--fasta-out/--table-out/--edge-table-out require --bubble-id; use --out-prefix when inspecting all bubbles");
    }

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(gfa_path, parse_options);
    if (graph.paths.empty()) {
        throw std::runtime_error("Input GFA has no P/W paths; inspect requires paths");
    }

    const std::vector<Bubble> bubbles = read_bubbles_csv(bubbles_csv_path);
    if (bubbles.empty()) {
        throw std::runtime_error("Bubbles CSV has no bubbles: " + bubbles_csv_path);
    }

    std::vector<const Bubble*> selected_bubbles;
    if (bubble_id == 0) {
        selected_bubbles.reserve(bubbles.size());
        for (const auto& bubble : bubbles) {
            selected_bubbles.push_back(&bubble);
        }
    } else {
        const auto bubble_it = std::find_if(
            bubbles.begin(),
            bubbles.end(),
            [&](const Bubble& b) { return b.id == bubble_id; });
        if (bubble_it == bubbles.end()) {
            throw std::runtime_error("Bubble ID not found in bubbles CSV: " + std::to_string(bubble_id));
        }
        selected_bubbles.push_back(&(*bubble_it));
    }

    std::size_t total_paths_written = 0;
    std::cout
        << "Input graph: " << gfa_path << "\n"
        << "Bubble source: " << bubbles_csv_path << "\n"
        << "Bubbles inspected: " << selected_bubbles.size() << "\n";

    for (const Bubble* bubble_ptr : selected_bubbles) {
        const Bubble& bubble = *bubble_ptr;
        std::string bubble_fasta_out_path = fasta_out_path;
        std::string bubble_table_out_path = table_out_path;
        std::string bubble_edge_table_out_path = edge_table_out_path;
        if (bubble_fasta_out_path.empty()) {
            bubble_fasta_out_path = out_prefix + ".bubble_" + std::to_string(bubble.id) + ".paths.fa.gz";
        }
        if (bubble_table_out_path.empty()) {
            bubble_table_out_path = out_prefix + ".bubble_" + std::to_string(bubble.id) + ".node_counts.tsv";
        }
        if (bubble_edge_table_out_path.empty()) {
            bubble_edge_table_out_path = out_prefix + ".bubble_" + std::to_string(bubble.id) + ".edge_counts.tsv";
        }

        const InspectBubbleResult result = write_inspect_outputs_for_bubble(
            graph,
            bubble,
            bubble_fasta_out_path,
            bubble_table_out_path,
            bubble_edge_table_out_path);
        total_paths_written += result.paths_written;

        std::cout
            << "Bubble ID: " << bubble.id << "\n"
            << "Source/sink: " << bubble.source << "/" << bubble.sink << "\n"
            << "Inside nodes: " << bubble.inside.size() << "\n"
            << "Paths written: " << result.paths_written << "\n"
            << "Wrote: " << result.fasta_out_path << "\n"
            << "Wrote: " << result.table_out_path << "\n"
            << "Wrote: " << result.edge_table_out_path << "\n";
    }

    std::cout << "Total paths written: " << total_paths_written << "\n";

    return 0;
}

} // namespace panvar
