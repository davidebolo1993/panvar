// Developer probe for the mosaic genotyper: a seeded paired-read simulator, a leave-one-out graph
// filter, and a truth evaluator. It is not part of the product surface and is never installed.
//
// Everything it measures is computed against a truth the product never sees. The evaluator rebuilds
// the model from the same GFA the run used and refuses unless that model reproduces the run's model
// fingerprint; truth paths are then projected from the full graph onto that model's block frame.

#include <sys/resource.h>
#if defined(__APPLE__)
#include <mach/mach.h>
#include <malloc/malloc.h>
#endif
#include "edlib.h"
#include "panvar/bubbles.hpp"
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_command.hpp"
#include "panvar/genotype_locus.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/genotype_repair.hpp"
#include "panvar/gfa.hpp"
#include "panvar/gfa_io.hpp"
#include "panvar/graph_sort.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/gz_reader.hpp"
#include "panvar/integrated_snarls.hpp"
#include "panvar/ls_hmm.hpp"
#include "panvar/marker_model.hpp"
#include "panvar/mosaic_beam.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <map>
#include <mutex>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace {

struct Args {
    std::map<std::string, std::string> values;

    std::string required(const std::string& key) const {
        const auto found = values.find(key);
        if (found == values.end()) throw std::runtime_error("missing required option --" + key);
        return found->second;
    }
    std::string get(const std::string& key, const std::string& fallback) const {
        const auto found = values.find(key);
        return found == values.end() ? fallback : found->second;
    }
};

Args parse_args(int argc, char** argv, int from) {
    Args args;
    for (int i = from; i < argc; ++i) {
        const std::string flag = argv[i];
        if (flag.rfind("--", 0) != 0 || i + 1 >= argc) {
            throw std::runtime_error("expected --key value pairs, got: " + flag);
        }
        args.values[flag.substr(2)] = argv[++i];
    }
    return args;
}

// fnv1a64_fmix64_top53_v1: FNV-1a over the bytes, finished with the murmur3 avalanche because FNV-1a alone
// leaves the high bits -- the ones every caller here reads -- poorly mixed for short sequential names.
std::uint64_t salted_hash(const std::string& text) {
    std::uint64_t hash = 1469598103934665603ULL;
    for (const unsigned char byte : text) {
        hash ^= byte;
        hash *= 1099511628211ULL;
    }
    hash ^= hash >> 33;
    hash *= 0xff51afd7ed558ccdULL;
    hash ^= hash >> 33;
    hash *= 0xc4ceb9fe1a85ec53ULL;
    hash ^= hash >> 33;
    return hash;
}

// A donor's own simulation stream: the master seed and the salt together, so two donors never draw the same
// reads and a salt-free run keeps the declared seed exactly.
std::uint64_t salted_seed(std::uint64_t seed, const std::string& salt) {
    return salt.empty() ? seed : salted_hash(salt + ":" + std::to_string(seed));
}

std::string revcomp(const std::string& bases) {
    std::string out(bases.size(), 'N');
    for (std::size_t i = 0; i < bases.size(); ++i) {
        const char base = bases[bases.size() - 1 - i];
        out[i] = base == 'A' ? 'T' : base == 'C' ? 'G' : base == 'G' ? 'C' : base == 'T' ? 'A' : 'N';
    }
    return out;
}

std::vector<const panvar::PathRecord*> donor_paths(const panvar::Graph& graph, const std::string& donor) {
    std::vector<const panvar::PathRecord*> out;
    for (const panvar::PathRecord& path : graph.paths) {
        if (path.name.rfind(donor + "#", 0) == 0) out.push_back(&path);
    }
    std::sort(out.begin(), out.end(),
              [](const panvar::PathRecord* lhs, const panvar::PathRecord* rhs) {
                  return lhs->name < rhs->name;
              });
    if (out.size() != 2) {
        throw std::runtime_error("donor " + donor + " must have exactly two haplotype paths, found " +
                                 std::to_string(out.size()));
    }
    return out;
}

struct SortedGraph {
    std::string resolved;
    panvar::GfaModel gfa;
    panvar::Graph graph;
    // input id -> sorted id, kept so a second graph can be renumbered into this one's frame.
    std::unordered_map<std::string, std::string> id_remap;
};

SortedGraph load_sorted(const std::string& gfa_path, const std::string& reference,
                        const std::string& stage_path) {
    SortedGraph out;
    out.gfa = panvar::read_gfa_model(gfa_path);
    panvar::GraphSortOptions sort_options;
    sort_options.reference_path = reference;
    sort_options.flip = true;
    const panvar::GraphSortResult sort_result = panvar::sort_graph_reference(out.gfa, sort_options);
    out.resolved = sort_result.resolved_reference;
    out.id_remap = sort_result.id_remap;
    panvar::write_gfa_model(stage_path, out.gfa);
    panvar::ParseGfaOptions parse_options;
    out.graph = panvar::parse_gfa(stage_path, parse_options);
    return out;
}

// Truth is a different graph from the model -- under leave-one-out it still carries the donor's two
// paths, which is exactly what makes the scheme leave-one-out -- but every truth projection runs
// through the MODEL's bubbles and block report, so the two graphs must number nodes identically.
//
// They do not do so by construction. sort_graph_reference builds its topological order from arcs
// taken off every path walk, so dropping a donor's two P lines removes arcs, changes in-degrees, and
// can change the order Kahn emits -- and with it the 1..N renumbering. On the C4 development panel
// this shifted 2926 of 4740 node ids for HG02818 while leaving HG00733 untouched: whether a donor
// moves the frame depends on whether a removed arc frees a node early enough to beat the numeric-id
// tie-break, so it cannot be predicted per donor. Leave-zero-out never exposed it because model and
// truth were then the same file.
//
// Everything else about the two sorts agrees -- same segments, same links, same reference path, so
// the same flip set and the same sequences -- which leaves the frames differing only by a
// permutation of ids. Recover the model's numbering by composing the two remaps:
// truth sorted id -> input id -> model sorted id.
SortedGraph load_truth_in_model_frame(const std::string& gfa_path, const std::string& reference,
                                      const SortedGraph& model_sorted, const std::string& stage_path) {
    SortedGraph out = load_sorted(gfa_path, reference, stage_path);
    std::unordered_map<std::string, std::string> to_model;
    to_model.reserve(out.id_remap.size() * 2);
    for (const auto& [input_id, truth_id] : out.id_remap) {
        const auto found = model_sorted.id_remap.find(input_id);
        if (found == model_sorted.id_remap.end()) {
            throw std::runtime_error("truth graph node " + input_id +
                                     " is absent from the model graph, so truth cannot be placed in "
                                     "the model's node frame");
        }
        to_model.emplace(truth_id, found->second);
    }
    const auto rm = [&](const std::string& id) -> std::string {
        const auto it = to_model.find(id);
        return it != to_model.end() ? it->second : id;
    };

    panvar::GfaModel& gfa = out.gfa;
    std::unordered_map<std::string, std::string> new_seq;
    new_seq.reserve(gfa.seq.size() * 2);
    for (const auto& [id, sequence] : gfa.seq) new_seq[rm(id)] = sequence;
    gfa.seq = std::move(new_seq);
    for (panvar::GfaEdge& e : gfa.edges) {
        e.from = rm(e.from);
        e.to = rm(e.to);
    }
    for (panvar::GfaPath& p : gfa.paths) {
        for (panvar::PathStep& s : p.steps) s.node_id = rm(s.node_id);
    }
    // node_order must stay the sorted order, which after renumbering is ascending id.
    gfa.node_order.clear();
    gfa.node_order.reserve(gfa.seq.size());
    for (const auto& [id, sequence] : gfa.seq) {
        (void)sequence;
        gfa.node_order.push_back(id);
    }
    std::sort(gfa.node_order.begin(), gfa.node_order.end(),
              [](const std::string& a, const std::string& b) {
                  return std::stoull(a) < std::stoull(b);
              });
    for (auto& [input_id, truth_id] : out.id_remap) {
        (void)input_id;
        truth_id = rm(truth_id);
    }

    // Re-stage so the written graph and the parsed one both carry the model's frame.
    panvar::write_gfa_model(stage_path, gfa);
    panvar::ParseGfaOptions parse_options;
    out.graph = panvar::parse_gfa(stage_path, parse_options);
    return out;
}

// Mirrors the command's structural stage exactly. The fingerprint check in `evaluate` is what makes
// this duplication safe: any divergence from the product refuses the evaluation.
struct LoadedModel {
    SortedGraph sorted;
    std::vector<panvar::Bubble> bubbles;
    panvar::LocusInputResult built;
    std::optional<panvar::PreparedGenotypeModel> model;
    std::string fingerprint;
};

// Mirrors read_curation_manifest in src/genotype_command.cpp: one template name per line, '#' opening
// a comment only when it OPENS the line, because template names contain '#' themselves. Duplicated
// rather than exported because the production parser is internal to the command; the rule that matters
// -- comments only at line start -- is gated on the production side by test_curated_panel. If the two
// ever disagree, a run and its evaluation would disagree about which panel was scored.
std::vector<std::string> read_panel_curation(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("panel curation manifest could not be opened: " + path);
    std::vector<std::string> names;
    std::string line;
    while (std::getline(file, line)) {
        const std::size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') continue;
        const std::size_t last = line.find_last_not_of(" \t\r\n");
        names.push_back(line.substr(first, last - first + 1));
    }
    if (names.empty()) throw std::runtime_error("panel curation manifest names no templates: " + path);
    return names;
}

LoadedModel load_model(const std::string& gfa_path, const std::string& reference,
                       const std::string& stage_path, std::size_t min_variant_bp,
                       const std::vector<std::string>& curated = {}) {
    LoadedModel out;
    out.sorted = load_sorted(gfa_path, reference, stage_path);
    panvar::validate_graph_paths(out.sorted.graph, "genotype_probe", true, true);
    panvar::BubbleCallOptions bubble_options;
    bubble_options.reference_path = out.sorted.resolved;
    bubble_options.quiet = true;
    bubble_options.min_variant_bp = min_variant_bp;
    bubble_options.snarl_pairs_override =
        panvar::find_top_level_snarls_cactus(panvar::snarl_input_from_model(out.sorted.gfa));
    bubble_options.snarl_source_supplied = true;
    out.bubbles = panvar::call_bubbles(out.sorted.graph, bubble_options);
    panvar::LocusInputOptions locus_options;
    locus_options.locus_id = out.sorted.resolved;
    locus_options.reference_path_name = out.sorted.resolved;
    // A curated run scored a curated model, so rebuilding it here without the same curation produces a
    // different fingerprint and every fingerprint check downstream refuses. Curation belongs to the
    // MODEL graph only -- truth is loaded separately and keeps every path, which is what lets a donor
    // be evaluated against full truth rather than against a curated view of it.
    locus_options.excluded_templates = curated;
    out.built = panvar::build_locus_input(out.sorted.graph, out.bubbles, locus_options);
    out.model.emplace(panvar::prepare_genotype_model(out.built.input));
    out.fingerprint = panvar::model_fingerprint(*out.model);
    return out;
}

std::map<std::string, std::string> read_key_values(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot read " + path);
    std::map<std::string, std::string> values;
    std::string line;
    std::getline(input, line);
    while (std::getline(input, line)) {
        const std::size_t tab = line.find('\t');
        if (tab != std::string::npos) values[line.substr(0, tab)] = line.substr(tab + 1);
    }
    return values;
}

std::vector<std::map<std::string, std::string>> read_table(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot read " + path);
    std::vector<std::map<std::string, std::string>> rows;
    std::string line;
    if (!std::getline(input, line)) return rows;
    std::vector<std::string> header;
    std::stringstream head(line);
    for (std::string cell; std::getline(head, cell, '\t');) header.push_back(cell);
    while (std::getline(input, line)) {
        std::map<std::string, std::string> row;
        std::stringstream fields(line);
        std::size_t column = 0;
        for (std::string cell; std::getline(fields, cell, '\t') && column < header.size(); ++column) {
            row[header[column]] = cell;
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

// Process resources and atomic replacement live in the library, shared with panvar genotype.
using panvar::FootprintSample;
using panvar::atomic_write;
using panvar::current_rss_bytes;
using panvar::peak_rss_gib;
using panvar::physical_footprint;
using panvar::release_free_pages;

std::vector<std::size_t> parse_alleles(const std::string& joined) {
    std::vector<std::size_t> out;
    std::stringstream stream(joined);
    for (std::string cell; std::getline(stream, cell, ',');) out.push_back(std::stoul(cell));
    return out;
}

// Global edit distance. Negative bound means unbounded; returns -1 when the distance exceeds it.
long edit_distance(const std::string& a, const std::string& b, long bound = -1) {
    if (a == b) return 0;
    if (a.empty() || b.empty()) {
        const long distance = static_cast<long>(std::max(a.size(), b.size()));
        return bound >= 0 && distance > bound ? -1 : distance;
    }
    EdlibAlignResult result = edlibAlign(
        a.data(), static_cast<int>(a.size()), b.data(), static_cast<int>(b.size()),
        edlibNewAlignConfig(static_cast<int>(bound), EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0));
    const long distance = result.status == EDLIB_STATUS_OK ? result.editDistance : -2;
    edlibFreeAlignResult(result);
    if (distance == -2) throw std::runtime_error("edlib failed");
    return distance;
}

// The availability floor for one truth sequence: the closest allele the model offers at this block.
long closest_allele_distance(const std::string& truth, const panvar::PreparedGenotypeBlock& block,
                             bool& exactly_available, std::size_t* closest = nullptr) {
    std::vector<std::size_t> alleles(block.alleles.size());
    for (std::size_t a = 0; a < alleles.size(); ++a) alleles[a] = a;
    // stable: equal length gaps keep allele order, so the closest allele reported is deterministic
    std::stable_sort(alleles.begin(), alleles.end(), [&](std::size_t lhs, std::size_t rhs) {
        const long dl = std::labs(static_cast<long>(block.alleles[lhs].sequence.size()) - static_cast<long>(truth.size()));
        const long dr = std::labs(static_cast<long>(block.alleles[rhs].sequence.size()) - static_cast<long>(truth.size()));
        return dl < dr;
    });
    long best = -1;
    exactly_available = false;
    for (const std::size_t index : alleles) {
        const std::string* allele = &block.alleles[index].sequence;
        const long length_gap =
            std::labs(static_cast<long>(allele->size()) - static_cast<long>(truth.size()));
        if (best >= 0 && length_gap >= best) break;   // length difference bounds edit distance below
        const long distance = edit_distance(truth, *allele, best < 0 ? -1 : best - 1);
        if (distance >= 0 && (best < 0 || distance < best)) {
            best = distance;
            if (closest != nullptr) *closest = index;
        }
        if (best == 0) {
            exactly_available = true;
            break;
        }
    }
    return best;
}

// ------------------------------------------------------------------------------------------------
// simulate
// ------------------------------------------------------------------------------------------------

int simulate(const Args& args) {
    const std::string out_dir = args.required("out");
    const std::string donor = args.required("donor");
    const double depth = std::stod(args.get("depth", "30"));
    const std::size_t read_length = std::stoul(args.get("read-length", "100"));
    const std::uint64_t declared_seed = std::stoull(args.get("seed", "7"));
    // Each donor gets its own stream when a salt is given; without one the declared seed is used unchanged.
    const std::string seed_salt = args.get("seed-salt", "");
    const std::uint64_t seed = salted_seed(declared_seed, seed_salt);

    panvar::LibraryModel library;
    library.insert_mean_bp = std::stod(args.get("insert-mean", "350"));
    library.insert_sd_bp = std::stod(args.get("insert-sd", "40"));
    library.support_sigmas = std::stod(args.get("sigmas", "4"));
    library.substitution_error = std::stod(args.get("error", "0.002"));
    library.outlier_weight = std::stod(args.get("outlier", "0.001"));
    // Signal templates per base per chromosome copy, summed over both orientations: per-copy coverage is
    // depth/2, and each template contributes two mates of read_length bases. The declared lambda is that
    // rate divided by (1-eta), because the likelihood's signal intensity is (1-eta) * lambda.
    const double signal_rate = (depth / 2.0) / (2.0 * static_cast<double>(read_length));
    library.fragment_start_intensity = signal_rate / (1.0 - library.outlier_weight);
    library.provenance = "simulation";
    panvar::validate_library_model(library);
    const panvar::DiscreteInsertDistribution insert =
        panvar::build_library_insert_distribution(library);

    const SortedGraph sorted =
        load_sorted(args.required("gfa"), args.required("reference"), out_dir + "/truth.sorted.gfa");
    const std::vector<const panvar::PathRecord*> paths = donor_paths(sorted.graph, donor);

    std::mt19937_64 rng(seed);
    std::ofstream r1(out_dir + "/reads_1.fq");
    std::ofstream r2(out_dir + "/reads_2.fq");
    std::ofstream truth(out_dir + "/truth.fa");
    std::ofstream manifest(out_dir + "/simulation.tsv");
    manifest << "key\tvalue\n" << "donor\t" << donor << "\nseed\t" << declared_seed << "\nseed_salt\t"
             << (seed_salt.empty() ? "none" : seed_salt) << "\neffective_seed\t" << seed << "\ndepth\t" << depth
             << "\nread_length\t" << read_length << "\nlambda\t" << library.fragment_start_intensity
             << "\nsimulator\tgenotype_probe simulate\n";

    std::size_t written = 0;
    std::size_t skipped_ambiguous = 0;
    for (std::size_t p = 0; p < paths.size(); ++p) {
        const std::string haplotype = panvar::spell_path_steps_sequence(sorted.graph, paths[p]->steps);
        truth << '>' << paths[p]->name << '\n' << haplotype << '\n';
        // The model's generative process: a Poisson number of templates with mean lambda * E_h, each
        // placed with probability proportional to P(i) * (L - i + 1) over inserts that can produce
        // two mates of this length, read in either orientation with probability 1/2.
        std::vector<std::size_t> lengths;
        std::vector<double> weights;
        for (std::size_t i = std::max(insert.support.min_bp, read_length);
             i <= insert.support.max_bp && i <= haplotype.size(); ++i) {
            lengths.push_back(i);
            weights.push_back(std::exp(insert.log_probability(i)) *
                              static_cast<double>(haplotype.size() - i + 1));
        }
        const double exposure = panvar::haplotype_exposure(haplotype.size(), insert, read_length);
        std::poisson_distribution<long> count_draw(signal_rate * exposure);
        std::discrete_distribution<std::size_t> length_draw(weights.begin(), weights.end());
        std::bernoulli_distribution orientation_draw(0.5);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        const long templates = count_draw(rng);
        manifest << "haplotype_" << p << "\t" << paths[p]->name << "\nhaplotype_" << p << "_bp\t"
                 << haplotype.size() << "\nhaplotype_" << p << "_expected_templates\t"
                 << signal_rate * exposure << "\nhaplotype_" << p
                 << "_templates\t" << templates << '\n';

        const auto mutate = [&](std::string& read) {
            for (char& base : read) {
                if (unit(rng) >= library.substitution_error) continue;
                const char* others = base == 'A' ? "CGT" : base == 'C' ? "AGT" : base == 'G' ? "ACT" : "ACG";
                base = others[std::uniform_int_distribution<int>(0, 2)(rng)];
            }
        };
        for (long t = 0; t < templates; ++t) {
            const std::size_t length = lengths[length_draw(rng)];
            const std::size_t start =
                std::uniform_int_distribution<std::size_t>(0, haplotype.size() - length)(rng);
            const std::string left = haplotype.substr(start, read_length);
            const std::string right = revcomp(haplotype.substr(start + length - read_length, read_length));
            const bool first_forward = orientation_draw(rng);
            std::string mate1 = first_forward ? left : right;
            std::string mate2 = first_forward ? right : left;
            mutate(mate1);
            mutate(mate2);
            if (mate1.find_first_not_of("ACGT") != std::string::npos ||
                mate2.find_first_not_of("ACGT") != std::string::npos) {
                ++skipped_ambiguous;
                continue;
            }
            ++written;
            const std::string quality(read_length, 'I');
            r1 << "@probe_" << written << "/1\n" << mate1 << "\n+\n" << quality << '\n';
            r2 << "@probe_" << written << "/2\n" << mate2 << "\n+\n" << quality << '\n';
        }
    }
    manifest << "fragments_written\t" << written << "\nfragments_skipped_non_acgt\t"
             << skipped_ambiguous << '\n';

    std::ofstream lib(out_dir + "/library.txt");
    lib << "# declared simulation parameters, genotype_probe simulate seed " << declared_seed
        << (seed_salt.empty() ? std::string() : " salted with " + seed_salt) << '\n'
        << "insert_mean_bp " << library.insert_mean_bp << "\ninsert_sd_bp " << library.insert_sd_bp
        << "\nsupport_sigmas " << library.support_sigmas << "\nsubstitution_error "
        << library.substitution_error << "\noutlier_weight " << library.outlier_weight
        << "\nfragment_start_intensity " << library.fragment_start_intensity << '\n';
    std::cout << "simulated " << written << " fragments for " << donor << " (" << skipped_ambiguous
              << " skipped for non-ACGT bases)\n";
    return 0;
}

// ------------------------------------------------------------------------------------------------
// loo-gfa
// ------------------------------------------------------------------------------------------------

int loo_gfa(const Args& args) {
    const std::string donor = args.required("donor");
    panvar::GzLineReader reader(args.required("gfa"));
    if (!reader.ok()) throw std::runtime_error("cannot open GFA");
    std::ofstream out(args.required("out"));
    std::size_t removed = 0;
    std::string line;
    while (reader.getline(line)) {
        if (line.rfind("P\t", 0) == 0) {
            const std::size_t end = line.find('\t', 2);
            if (line.substr(2, end - 2).rfind(donor + "#", 0) == 0) {
                ++removed;
                continue;
            }
        } else if (line.rfind("W\t", 0) == 0) {
            const std::size_t end = line.find('\t', 2);
            if (line.substr(2, end - 2) == donor) {
                ++removed;
                continue;
            }
        }
        out << line << '\n';
    }
    if (removed != 2) {
        throw std::runtime_error("leave-one-out must remove exactly the donor's two paths, removed " +
                                 std::to_string(removed));
    }
    std::cout << "removed " << removed << " paths of " << donor << '\n';
    return 0;
}

// ------------------------------------------------------------------------------------------------
// evaluate
// ------------------------------------------------------------------------------------------------

int evaluate(const Args& args) {
    const std::string run = args.required("run-prefix");
    const std::string out = args.required("out");
    const std::string scheme = args.required("scheme");
    const std::string donor = args.required("donor");
    const std::string reference = args.required("reference");
    const std::map<std::string, std::string> status = read_key_values(run + ".mosaic_status.tsv");

    std::ofstream summary(out + ".summary.tsv");
    summary << "key\tvalue\nscheme\t" << scheme << "\ndonor\t" << donor << "\nrun\t" << run << '\n';
    for (const char* key : {"status", "note", "max_mosaic_candidates", "mismatch_band", "proposal_state_limit",
                            "mosaic_candidates", "mosaic_candidates_found", "proposal_states_pruned",
                            "proposal_search_exhaustive", "mosaic_target_reached", "panel_diplotypes",
                            "panel_baseline_complete", "candidate_scope", "phase_support_within_candidates",
                            "phase_audit_status",
                            "candidate_set_conditional", "sequence_candidates", "phase_support",
                            "fragment_contenders", "blocks_genotype_read_supported", "numeric_separation_gap",
                            "fragment_term_gap_to_runner_up", "winner_candidate_set_probability",
                            "seconds_proposal", "seconds_total", "seconds_rescoring", "model_fingerprint",
                            "marker_depth_rate_per_copy", "marker_background"}) {
        const auto found = status.find(key);
        summary << key << '\t' << (found == status.end() ? "NA" : found->second) << '\n';
    }

    const std::string curation_path = args.get("panel-curation", "");
    const LoadedModel loaded =
        load_model(args.required("model-gfa"), reference, out + ".model.sorted.gfa", 50,
                   curation_path.empty() ? std::vector<std::string>{}
                                         : read_panel_curation(curation_path));
    if (status.count("model_fingerprint") == 0 || status.at("model_fingerprint") != loaded.fingerprint) {
        throw std::runtime_error("the rebuilt model does not reproduce the run's model fingerprint" +
                                 std::string(curation_path.empty()
                                                 ? "; if the run was curated, pass the same "
                                                   "--panel-curation manifest"
                                                 : " even with the curation manifest given"));
    }
    summary << "panel_curation\t" << (curation_path.empty() ? "none" : curation_path) << '\n';
    const panvar::PreparedGenotypeModel& model = *loaded.model;

    // Truth comes from the full graph, placed in the model's node frame -- left to sort itself it does
    // not land in the same frame once the donor's paths are gone. The check below then confirms the
    // frame rather than discovering it is broken.
    const SortedGraph truth_graph = load_truth_in_model_frame(
        args.required("truth-gfa"), reference, loaded.sorted, out + ".truth.sorted.gfa");
    for (const auto& item : loaded.sorted.gfa.seq) {
        const auto found = truth_graph.gfa.seq.find(item.first);
        if (found == truth_graph.gfa.seq.end() || found->second != item.second) {
            throw std::runtime_error("model and truth graphs do not share a node frame at node " +
                                     item.first);
        }
    }
    const std::vector<const panvar::PathRecord*> truth_paths = donor_paths(truth_graph.graph, donor);
    std::vector<panvar::PathProjection> truth;
    for (const panvar::PathRecord* path : truth_paths) {
        const auto projection = panvar::project_path_onto_locus(
            truth_graph.graph, loaded.bubbles, loaded.built.report, path->name);
        if (!projection.has_value()) {
            throw std::runtime_error("truth path " + path->name + " does not span the block chain");
        }
        truth.push_back(*projection);
    }

    // In a leave-zero-out run the donor is a panel template, so its projection must be exactly its own
    // row of the model. This is the check that the projection frame is the model frame.
    std::size_t lzo_rows_checked = 0;
    for (std::size_t t = 0; t < truth_paths.size(); ++t) {
        for (std::size_t row = 0; row < model.panel().template_count(); ++row) {
            if (model.panel().template_names[row] != truth_paths[t]->name) continue;
            for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                const std::string& allele =
                    model.blocks()[b].alleles[model.panel().template_alleles[row][b]].sequence;
                if (!truth[t].blocks[b].has_value() || *truth[t].blocks[b] != allele) {
                    throw std::runtime_error("projection of panel template " + truth_paths[t]->name +
                                             " disagrees with its model row at block " +
                                             model.blocks()[b].block_id);
                }
            }
            ++lzo_rows_checked;
        }
    }
    summary << "truth_rows_checked_against_model\t" << lzo_rows_checked << '\n';

    // Truth availability: can the model spell each truth haplotype at all?
    std::vector<bool> spellable(2, true);
    std::size_t blocks_with_both_available = 0;

    std::ofstream table(out + ".blocks.tsv");
    table << "scheme\tdonor\tmax_mosaic_candidates\tblock_id\tkind\ttruth1_bp\ttruth2_bp\tcalled1_bp\t"
             "called2_bp\tunordered_distance\tavailability_floor\texcess\ttruth1_available\t"
             "truth2_available\tgenotype_support\tcalled_allele1\tcalled_allele2\ttruth1_closest_allele\t"
             "truth2_closest_allele\ttruth1_closest_distance\ttruth2_closest_distance\n";

    // A genotype-equivalent call emits one unordered genotype shared by its whole top set, and the
    // distances below are unordered per block, so any representative of the set measures it exactly.
    // DIAGNOSTIC FALLBACK ONLY: --fallback-alleles <hap1>|<hap2> evaluates a declared diplotype -- the
    // initial union's top candidate of a run that emitted no genotype -- instead of the run's call. The
    // summary labels it, and it is never an official call: callers keep it out of called-only accuracy.
    const std::string fallback = args.get("fallback-alleles", "");
    const bool use_fallback = !fallback.empty();
    const bool complete = use_fallback ||
                          (status.count("status") &&
                           (status.at("status") == "COMPLETE" || status.at("status") == "CONDITIONAL_GENOTYPE_EQUIVALENT"));
    std::vector<std::size_t> hap1;
    std::vector<std::size_t> hap2;
    std::map<std::string, std::string> genotype_support;
    const std::string phase_support = status.count("phase_support") ? status.at("phase_support") : "NA";
    // Phase-aware accuracy only where the reads chose the phase: within the candidates, or audited.
    const bool phase_claimed = phase_support == "READS_WITHIN_CANDIDATES" || phase_support == "READS_AUDITED";
    if (use_fallback) {
        const std::size_t bar = fallback.find('|');
        if (bar == std::string::npos) throw std::runtime_error("--fallback-alleles wants <hap1>|<hap2>");
        hap1 = parse_alleles(fallback.substr(0, bar));
        hap2 = parse_alleles(fallback.substr(bar + 1));
        if (hap1.size() != model.blocks().size() || hap2.size() != model.blocks().size()) {
            throw std::runtime_error("--fallback-alleles needs one allele per block on both haplotypes");
        }
    } else if (complete) {
        for (const auto& row : read_table(run + ".block_genotypes.tsv")) {
            genotype_support[row.at("block_id")] = row.at("genotype_support");
        }
        std::string representative;
        for (const auto& row : read_table(run + ".mosaic_scores.tsv")) {
            if (row.at("in_winner_set") == "true") representative = row.at("representative_beam_candidate");
        }
        for (const auto& row : read_table(run + ".mosaic_candidates.tsv")) {
            if (row.at("candidate_id") == representative) {
                hap1 = parse_alleles(row.at("hap1_alleles"));
                hap2 = parse_alleles(row.at("hap2_alleles"));
            }
        }
        if (hap1.size() != model.blocks().size() || hap2.size() != model.blocks().size()) {
            throw std::runtime_error("could not recover the winner's allele paths from the run");
        }
    }

    long unordered_total = 0;
    long floor_total = 0;
    long sv_distance_total = 0;
    long context_distance_total = 0;
    std::size_t read_supported_blocks = 0;
    std::size_t read_supported_exact = 0;
    long read_supported_distance = 0;
    long read_supported_floor = 0;
    long order_a = 0;
    long order_b = 0;
    std::size_t sv_exact = 0;
    std::size_t sv_blocks = 0;
    std::size_t all_exact = 0;
    std::size_t undefined_blocks = 0;
    std::vector<int> local_best;   // +1 when direct pairing is strictly better, -1 swapped, 0 tie
    for (std::size_t b = 0; b < model.blocks().size(); ++b) {
        const panvar::PreparedGenotypeBlock& block = model.blocks()[b];
        const bool variant = block.kind == panvar::GenotypeBlockKind::Variant;
        if (variant) ++sv_blocks;
        const auto& t1 = truth[0].blocks[b];
        const auto& t2 = truth[1].blocks[b];
        if (!t1.has_value() || !t2.has_value()) {
            ++undefined_blocks;
            spellable[0] = spellable[0] && t1.has_value();
            spellable[1] = spellable[1] && t2.has_value();
            table << scheme << '\t' << donor << '\t' << status.at("max_mosaic_candidates") << '\t'
                  << block.block_id << '\t' << (variant ? "SV" : "context")
                  << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
            continue;
        }
        bool available1 = false;
        bool available2 = false;
        // The floor pair: the allele closest to each truth haplotype (first by allele order among equals).
        std::size_t closest1 = 0;
        std::size_t closest2 = 0;
        const long floor1 = closest_allele_distance(*t1, block, available1, &closest1);
        const long floor2 = closest_allele_distance(*t2, block, available2, &closest2);
        spellable[0] = spellable[0] && available1;
        spellable[1] = spellable[1] && available2;
        if (available1 && available2) ++blocks_with_both_available;
        const long floor = floor1 + floor2;
        floor_total += floor;

        std::string called1_bp = "NA";
        std::string called2_bp = "NA";
        std::string distance_cell = "NA";
        std::string excess_cell = "NA";
        if (complete) {
            const std::string& c1 = block.alleles[hap1[b]].sequence;
            const std::string& c2 = block.alleles[hap2[b]].sequence;
            const long direct = edit_distance(*t1, c1) + edit_distance(*t2, c2);
            const long swapped = edit_distance(*t1, c2) + edit_distance(*t2, c1);
            const long unordered = std::min(direct, swapped);
            unordered_total += unordered;
            (variant ? sv_distance_total : context_distance_total) += unordered;
            const auto support = genotype_support.find(block.block_id);
            if (support != genotype_support.end() && support->second == "READS") {
                ++read_supported_blocks;
                if (unordered == 0) ++read_supported_exact;
                read_supported_distance += unordered;
                read_supported_floor += floor;
            }
            order_a += direct;
            order_b += swapped;
            local_best.push_back(direct < swapped ? 1 : swapped < direct ? -1 : 0);
            if (unordered == 0) {
                ++all_exact;
                if (variant) ++sv_exact;
            }
            called1_bp = std::to_string(c1.size());
            called2_bp = std::to_string(c2.size());
            distance_cell = std::to_string(unordered);
            excess_cell = std::to_string(unordered - floor);
        }
        table << scheme << '\t' << donor << '\t' << status.at("max_mosaic_candidates") << '\t' << block.block_id
              << '\t' << (variant ? "SV" : "context") << '\t' << t1->size() << '\t' << t2->size() << '\t'
              << called1_bp << '\t' << called2_bp << '\t' << distance_cell << '\t' << floor << '\t'
              << excess_cell << '\t' << (available1 ? "true" : "false") << '\t'
              << (available2 ? "true" : "false") << '\t'
              << (genotype_support.count(block.block_id) ? genotype_support.at(block.block_id) : "NA") << '\t'
              << (complete ? std::to_string(hap1[b]) : std::string("NA")) << '\t'
              << (complete ? std::to_string(hap2[b]) : std::string("NA")) << '\t' << closest1 << '\t' << closest2
              << '\t' << floor1 << '\t' << floor2 << '\n';
    }

    // The only no-call: no genotype result exists.
    summary << "run_no_call\t" << (complete ? "false" : "true") << '\n';
    summary << "evaluated_diplotype\t"
            << (use_fallback ? args.get("evaluated-label", "diagnostic_fallback_not_a_call") : std::string("the_run_call"))
            << '\n';
    summary << "blocks\t" << model.blocks().size() << "\nsv_blocks\t" << sv_blocks
            << "\nundefined_truth_blocks\t" << undefined_blocks
            << "\nblocks_with_both_truth_alleles_available\t" << blocks_with_both_available
            << "\ntruth_hap1_spellable\t" << (spellable[0] ? "true" : "false")
            << "\ntruth_hap2_spellable\t" << (spellable[1] ? "true" : "false")
            << "\navailability_floor_total\t" << floor_total << '\n';
    if (complete) {
        const int global = order_a <= order_b ? 1 : -1;
        std::size_t switch_blocks = 0;
        for (int local : local_best) {
            if (local != 0 && local != global) ++switch_blocks;
        }
        summary << "unordered_distance_total\t" << unordered_total
                << "\nexcess_over_floor_total\t" << unordered_total - floor_total
                // Phase-aware accuracy only where the reads chose the phase.
                << "\nphased_distance_total\t"
                << (phase_claimed ? std::to_string(std::min(order_a, order_b)) : std::string("NA"))
                << "\nphase_switch_blocks\t"
                << (phase_claimed ? std::to_string(switch_blocks) : std::string("NA"))
                << "\nblocks_exact\t" << all_exact << "\nsv_blocks_exact\t" << sv_exact
                << "\nsv_unordered_distance_total\t" << sv_distance_total
                << "\ncontext_unordered_distance_total\t" << context_distance_total
                // Accuracy over every emitted block genotype is above; this is the read-supported subset.
                << "\nread_supported_blocks\t" << read_supported_blocks
                << "\nread_supported_blocks_exact\t" << read_supported_exact
                << "\nread_supported_unordered_distance_total\t" << read_supported_distance
                << "\nread_supported_excess_over_floor_total\t" << read_supported_distance - read_supported_floor
                // Neither of these is a no-call: a genotype was emitted.
                << "\nblock_unsupported\t"
                << std::count_if(genotype_support.begin(), genotype_support.end(),
                                 [](const auto& item) { return item.second != "READS"; })
                << "\nphase_unsupported\t" << (phase_claimed ? "false" : "true") << '\n';

        // The candidate-set diagnostics below read the run's candidate and score tables, which a run that
        // emitted no genotype never writes: a fallback evaluation has no candidate set to describe.
        if (use_fallback) {
            summary << "candidate_set_diagnostics\tnot_applicable_to_a_fallback\n";
        } else {
        // Truth in candidate set: some retained beam candidate spells both truth haplotypes.
        std::string truth_proposal_rank = "absent";
        std::string truth_sequence_candidate = "absent";
        std::string truth_source = "absent";
        for (const auto& row : read_table(run + ".mosaic_candidates.tsv")) {
            const std::vector<std::size_t> a = parse_alleles(row.at("hap1_alleles"));
            const std::vector<std::size_t> c = parse_alleles(row.at("hap2_alleles"));
            const auto spells = [&](const std::vector<std::size_t>& path, std::size_t which) {
                for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                    if (!truth[which].blocks[b].has_value() ||
                        model.blocks()[b].alleles[path[b]].sequence != *truth[which].blocks[b]) {
                        return false;
                    }
                }
                return true;
            };
            if ((spells(a, 0) && spells(c, 1)) || (spells(a, 1) && spells(c, 0))) {
                truth_proposal_rank = row.at("unique_rank");
                truth_sequence_candidate = row.at("sequence_candidate_id");
                truth_source = row.count("source") ? row.at("source") : "mosaic";
                break;
            }
        }
        std::string truth_score_rank = "absent";
        if (truth_sequence_candidate != "absent") {
            std::size_t rank = 0;
            for (const auto& row : read_table(run + ".mosaic_scores.tsv")) {
                ++rank;
                if (row.at("sequence_candidate_id") == truth_sequence_candidate) {
                    truth_score_rank = std::to_string(rank);
                }
            }
        }
        // Named for what it actually measures. candidate_prefix() truncates the beam to K BEFORE
        // mosaic_candidates.tsv is written, and the loop above searches only that file, so an "absent"
        // verdict means absent from the K rows written for this run. It cannot distinguish a truth
        // ranked beyond K in the final beam, one pruned at an earlier layer, and one never generated.
        // The former name, truth_in_candidates, invited exactly that conflation: a round-14 conclusion
        // was drawn from it and had to be withdrawn.
        summary << "truth_in_ranked_prefix\t" << (truth_proposal_rank != "absent" ? "true" : "false")
                << "\ntruth_ranked_prefix_limit\t"
                << (status.count("max_mosaic_candidates") ? status.at("max_mosaic_candidates")
                                                          : std::string("NA"))
                << "\ntruth_candidate_source\t" << truth_source
                << "\ntruth_ranked_prefix_rank\t" << truth_proposal_rank
                << "\ntruth_fragment_score_rank\t" << truth_score_rank << '\n';

        // Availability diagnostics -- evaluation only: they read the truth and never feed a call. The
        // candidate closest to the truth, whether it attains the availability floor, and how far the
        // likelihood placed it below the winner. Among several closest candidates, the highest-scoring.
        std::map<std::string, double> score_of;
        std::map<std::string, double> fragment_of;
        std::string winner_sequence;
        for (const auto& row : read_table(run + ".mosaic_scores.tsv")) {
            score_of[row.at("sequence_candidate_id")] = std::stod(row.at("score"));
            fragment_of[row.at("sequence_candidate_id")] = std::stod(row.at("fragment_term"));
            if (row.at("in_winner_set") == "true") winner_sequence = row.at("sequence_candidate_id");
        }
        std::map<std::tuple<std::size_t, std::size_t, std::size_t>, long> distance_cache;
        const auto truth_distance = [&](std::size_t which, std::size_t b, std::size_t allele) {
            const auto key = std::make_tuple(which, b, allele);
            const auto found = distance_cache.find(key);
            if (found != distance_cache.end()) return found->second;
            const long value = edit_distance(*truth[which].blocks[b], model.blocks()[b].alleles[allele].sequence);
            distance_cache.emplace(key, value);
            return value;
        };
        long minimum = -1;
        std::string closest = "NA";
        std::string closest_rank = "NA";
        for (const auto& row : read_table(run + ".mosaic_candidates.tsv")) {   // in proposal rank order
            const std::vector<std::size_t> a = parse_alleles(row.at("hap1_alleles"));
            const std::vector<std::size_t> c = parse_alleles(row.at("hap2_alleles"));
            long total = 0;
            for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                if (!truth[0].blocks[b].has_value() || !truth[1].blocks[b].has_value()) continue;
                total += std::min(truth_distance(0, b, a[b]) + truth_distance(1, b, c[b]),
                                  truth_distance(0, b, c[b]) + truth_distance(1, b, a[b]));
            }
            const std::string id = row.at("sequence_candidate_id");
            if (minimum < 0 || total < minimum ||
                (total == minimum && score_of.at(id) > score_of.at(closest))) {
                minimum = total;
                closest = id;
                closest_rank = row.at("unique_rank");
            }
        }
        summary << "availability_floor\t" << floor_total << "\nminimum_candidate_distance\t" << minimum
                << "\nfloor_attained_in_candidates\t" << (minimum == floor_total ? "true" : "false")
                << "\nclosest_candidate_proposal_rank\t" << closest_rank
                << "\nclosest_candidate_sequence_id\t" << closest << std::setprecision(12)
                << "\nclosest_candidate_fragment_gap_to_winner\t"
                << fragment_of.at(winner_sequence) - fragment_of.at(closest)
                << "\nclosest_candidate_total_score_gap_to_winner\t"
                << score_of.at(winner_sequence) - score_of.at(closest) << std::setprecision(6)
                << "\nwinner_excess_over_floor\t" << unordered_total - floor_total << '\n';
        }
    } else {
        summary << "accuracy\texcluded: run status is not COMPLETE\n";
    }
    std::cout << "evaluated " << scheme << " " << donor << " K="
              << (status.count("max_mosaic_candidates") ? status.at("max_mosaic_candidates") : std::string("NA")) << '\n';
    return 0;
}


// ------------------------------------------------------------------------------------------------
// compare: a fixed pair of diplotypes, decomposed fragment by fragment
// ------------------------------------------------------------------------------------------------

std::vector<panvar::ReadPair> read_pairs(const std::string& first, const std::string& second) {
    const auto load = [](const std::string& file) {
        panvar::GzLineReader reader(file);
        if (!reader.ok()) throw std::runtime_error("cannot read " + file);
        std::vector<panvar::SequencingRead> reads;
        std::string header, bases, plus, quality;
        while (reader.getline(header)) {
            if (header.empty()) continue;
            reader.getline(bases);
            reader.getline(plus);
            reader.getline(quality);
            std::string name = header.substr(1);
            const std::size_t space = name.find_first_of(" \t");
            if (space != std::string::npos) name.resize(space);
            reads.push_back({name, bases, quality});
        }
        return reads;
    };
    const auto a = load(first);
    const auto b = load(second);
    if (a.size() != b.size()) throw std::runtime_error("paired read files differ in record count");
    std::vector<panvar::ReadPair> out;
    for (std::size_t i = 0; i < a.size(); ++i) out.push_back({a[i], b[i]});
    return out;
}

// ------------------------------------------------------------------------------------------------
// subsample: nested lower-coverage read sets from one simulation
// ------------------------------------------------------------------------------------------------

// Keeps a paired fragment when the unit value of a salted hash of its name falls below
// coverage/source_coverage. The decision depends only on the salt and the name, so the retained sets are nested
// within a donor -- every fragment kept at a coverage is kept at every higher one -- while two donors given
// different salts get different masks. Resimulating each coverage would not nest, and would confound coverage
// with a fresh draw.
int subsample(const Args& args) {
    const std::string out_dir = args.required("out");
    const double source_coverage = std::stod(args.required("source-coverage"));
    const double coverage = std::stod(args.required("coverage"));
    if (!(source_coverage > 0.0) || !(coverage > 0.0) || coverage > source_coverage) {
        throw std::runtime_error("coverage must be positive and at most the source coverage");
    }
    const double fraction = coverage / source_coverage;
    const std::string name_salt = args.get("name-salt", "");
    const std::vector<panvar::ReadPair> fragments =
        read_pairs(args.required("reads1"), args.required("reads2"));

    std::ofstream r1(out_dir + "/reads_1.fq");
    std::ofstream r2(out_dir + "/reads_2.fq");
    std::size_t kept = 0;
    for (const panvar::ReadPair& fragment : fragments) {
        // Both mates carry the same fragment name up to its /1 or /2 suffix.
        std::string name = fragment.first.name;
        const std::size_t slash = name.rfind('/');
        if (slash != std::string::npos) name.resize(slash);
        // The top 53 bits as a unit value: exact in double, and monotone in the hash.
        const double unit =
            static_cast<double>(salted_hash(name_salt + "/" + name) >> 11) * (1.0 / 9007199254740992.0);
        if (unit >= fraction) continue;
        ++kept;
        r1 << '@' << fragment.first.name << '\n' << fragment.first.bases << "\n+\n" << fragment.first.qualities << '\n';
        r2 << '@' << fragment.second.name << '\n' << fragment.second.bases << "\n+\n" << fragment.second.qualities << '\n';
    }

    // The declared library follows the coverage: the fragment-start intensity is proportional to it.
    std::map<std::string, std::string> library;
    std::vector<std::string> order;
    {
        std::ifstream input(args.required("library"));
        if (!input) throw std::runtime_error("cannot read " + args.required("library"));
        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream fields(line);
            std::string key, value;
            if (fields >> key >> value) {
                if (!library.count(key)) order.push_back(key);
                library[key] = value;
            }
        }
    }
    if (!library.count("fragment_start_intensity")) {
        throw std::runtime_error("the source library must declare fragment_start_intensity");
    }
    const double source_lambda = std::stod(library.at("fragment_start_intensity"));
    {
        std::ofstream lib(out_dir + "/library.txt");
        lib << "# declared parameters, genotype_probe subsample to " << coverage << "x of " << source_coverage
            << "x\n"
            << std::setprecision(12);
        for (const std::string& key : order) {
            if (key == "fragment_start_intensity") continue;
            lib << key << ' ' << library.at(key) << '\n';
        }
        lib << "fragment_start_intensity " << source_lambda * fraction << '\n';
    }
    {
        std::ofstream manifest(out_dir + "/subsample.tsv");
        manifest << std::setprecision(12) << "key\tvalue\nsource_coverage\t" << source_coverage << "\ncoverage\t"
                 << coverage << "\nrequested_fraction\t" << fraction << "\nsource_fragments\t" << fragments.size()
                 << "\nkept_fragments\t" << kept << "\nobserved_fraction\t"
                 << (fragments.empty() ? 0.0 : static_cast<double>(kept) / static_cast<double>(fragments.size()))
                 << "\nfragment_start_intensity\t" << source_lambda * fraction
                 << "\nname_salt\t" << (name_salt.empty() ? "none" : name_salt)
                 << "\nselector\tfnv1a64_fmix64_top53_v1 over <salt>/<fragment name>, nested across coverages\n";
    }
    std::cout << "subsampled " << kept << " of " << fragments.size() << " fragments to " << coverage << "x\n";
    return 0;
}


panvar::LibraryModel read_library_file(const std::string& file) {
    std::ifstream input(file);
    if (!input) throw std::runtime_error("cannot read " + file);
    std::map<std::string, double> values;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::string key;
        double value = 0.0;
        if (fields >> key >> value) values[key] = value;
    }
    panvar::LibraryModel model;
    model.insert_mean_bp = values.at("insert_mean_bp");
    model.insert_sd_bp = values.at("insert_sd_bp");
    model.support_sigmas = values.at("support_sigmas");
    model.substitution_error = values.at("substitution_error");
    model.outlier_weight = values.at("outlier_weight");
    model.fragment_start_intensity = values.at("fragment_start_intensity");
    model.provenance = "declared: " + file;
    panvar::validate_library_model(model);
    return model;
}

std::pair<std::vector<std::size_t>, std::vector<std::size_t>> winner_paths(const std::string& run) {
    std::string representative;
    for (const auto& row : read_table(run + ".mosaic_scores.tsv")) {
        if (row.at("in_winner_set") == "true") representative = row.at("representative_beam_candidate");
    }
    for (const auto& row : read_table(run + ".mosaic_candidates.tsv")) {
        if (row.at("candidate_id") == representative) {
            return {parse_alleles(row.at("hap1_alleles")), parse_alleles(row.at("hap2_alleles"))};
        }
    }
    throw std::runtime_error("no winner in " + run);
}

// The model's generative process on arbitrary haplotype sequences, in memory.
std::vector<panvar::ReadPair> simulate_from(const std::vector<std::string>& haplotypes,
                                            const panvar::LibraryModel& library, double depth,
                                            std::size_t read_length, std::uint64_t seed) {
    const panvar::DiscreteInsertDistribution insert = panvar::build_library_insert_distribution(library);
    const double signal_rate = (depth / 2.0) / (2.0 * static_cast<double>(read_length));
    std::mt19937_64 rng(seed);
    std::vector<panvar::ReadPair> out;
    std::size_t serial = 0;
    for (const std::string& haplotype : haplotypes) {
        std::vector<std::size_t> lengths;
        std::vector<double> weights;
        for (std::size_t i = std::max(insert.support.min_bp, read_length);
             i <= insert.support.max_bp && i <= haplotype.size(); ++i) {
            lengths.push_back(i);
            weights.push_back(std::exp(insert.log_probability(i)) * static_cast<double>(haplotype.size() - i + 1));
        }
        const double exposure = panvar::haplotype_exposure(haplotype.size(), insert, read_length);
        std::poisson_distribution<long> count(signal_rate * exposure);
        std::discrete_distribution<std::size_t> length_draw(weights.begin(), weights.end());
        std::bernoulli_distribution orientation(0.5);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        const long templates = count(rng);
        for (long t = 0; t < templates; ++t) {
            const std::size_t length = lengths[length_draw(rng)];
            const std::size_t start = std::uniform_int_distribution<std::size_t>(0, haplotype.size() - length)(rng);
            const std::string left = haplotype.substr(start, read_length);
            const std::string right = revcomp(haplotype.substr(start + length - read_length, read_length));
            const bool first_forward = orientation(rng);
            std::string m1 = first_forward ? left : right;
            std::string m2 = first_forward ? right : left;
            for (std::string* mate : {&m1, &m2}) {
                for (char& base : *mate) {
                    if (unit(rng) >= library.substitution_error) continue;
                    const char* others = base == 'A' ? "CGT" : base == 'C' ? "AGT" : base == 'G' ? "ACT" : "ACG";
                    base = others[std::uniform_int_distribution<int>(0, 2)(rng)];
                }
            }
            if (m1.find_first_not_of("ACGT") != std::string::npos || m2.find_first_not_of("ACGT") != std::string::npos) {
                continue;
            }
            ++serial;
            const std::string q(read_length, 'I');
            out.push_back({{"cal_" + std::to_string(serial) + "/1", m1, q}, {"cal_" + std::to_string(serial) + "/2", m2, q}});
        }
    }
    return out;
}

// Every accepted origin found by scanning every start position -- no seed index. The Hamming count and
// FR geometry are the shared primitives; the search is independent of the pigeonhole index.
// Every field a placement carries must be filled here, not left defaulted. log_placement_mass no
// longer recomputes an emission from a mismatch count, so a producer that forgets the emission gets
// zeros and a SILENTLY wrong mass rather than a compile error -- and this scan exists precisely to
// cross-check the accelerated scorer, so a silent zero would make the two paths agree for the wrong
// reason.
std::vector<panvar::FragmentPlacement> scan_origins(const std::string& haplotype, const panvar::ReadPair& fragment,
                                                    const panvar::InsertSupport& support,
                                                    const panvar::LibraryModel& library, std::size_t band) {
    struct Hit {
        std::size_t start = 0;
        panvar::GappedFit fit;
    };
    // The forward mate is anchored at its START and the reverse mate at its RIGHT ENDPOINT, because
    // that is how the template was sequenced and how the production path scores it. Anchoring both
    // at a start would make this cross-check disagree with production wherever an alignment is
    // gapped -- silently, since both would still be proper distributions.
    const auto starts = [&](const std::string& read) {
        std::vector<Hit> out;
        if (haplotype.empty()) return out;
        for (std::size_t p = 0; p < haplotype.size(); ++p) {
            const auto fit = panvar::gapped_emission(haplotype, read, p, library, band);
            if (fit.has_value() && std::isfinite(fit->log_emission)) out.push_back({p, *fit});
        }
        return out;
    };
    const auto ends = [&](const std::string& read) {
        std::vector<Hit> out;
        if (haplotype.empty()) return out;
        for (std::size_t p = 0; p <= haplotype.size(); ++p) {
            const auto fit = panvar::gapped_emission_ending_at(haplotype, read, p, library, band);
            if (fit.has_value() && std::isfinite(fit->log_emission)) out.push_back({p, *fit});
        }
        return out;
    };
    const std::string rc1 = revcomp(fragment.first.bases);
    const std::string rc2 = revcomp(fragment.second.bases);
    std::vector<panvar::FragmentPlacement> out;
    for (int arm = 0; arm < 2; ++arm) {
        const std::string& forward = arm == 0 ? fragment.first.bases : fragment.second.bases;
        const std::string& reverse = arm == 0 ? rc2 : rc1;
        const auto f_starts = starts(forward);
        const auto r_ends = ends(reverse);
        // (left endpoint, insert) -> the summed mass of every latent alignment pair reaching those
        // two template ends. This scan exists to cross-check the production collapse, so it must
        // collapse the same way: alternative alignments of one template are ONE origin, and pairing
        // them without collapsing would report alignment ambiguity as extra physical copies.
        struct Accumulated {
            double log_latent_mass = -std::numeric_limits<double>::infinity();
            std::size_t latent_alignments = 0;
            double best_pair = -std::numeric_limits<double>::infinity();
            std::size_t reverse_start = 0;
            std::size_t forward_consumed_bp = 0;
            std::size_t reverse_consumed_bp = 0;
            std::size_t forward_mismatches = 0;
            std::size_t reverse_mismatches = 0;
        };
        std::map<std::pair<std::size_t, std::size_t>, Accumulated> origins;
        for (const auto& fs : f_starts) {
            for (const auto& rs : r_ends) {
                for (const auto& forward_span : fs.fit.components) {
                    for (const auto& reverse_span : rs.fit.components) {
                        const auto insert = panvar::fr_insert_from_endpoints(
                            fs.start, forward_span.consumed_bp,
                            rs.start, reverse_span.consumed_bp, support);
                        if (!insert.has_value()) continue;
                        Accumulated& origin = origins[{fs.start, *insert}];
                        const double pair = forward_span.log_emission + reverse_span.log_emission;
                        origin.log_latent_mass = panvar::log_add_exp(origin.log_latent_mass, pair);
                        ++origin.latent_alignments;
                        if (pair <= origin.best_pair) continue;
                        origin.best_pair = pair;
                        origin.reverse_start = rs.start - reverse_span.consumed_bp;
                        origin.forward_consumed_bp = forward_span.consumed_bp;
                        origin.reverse_consumed_bp = reverse_span.consumed_bp;
                        origin.forward_mismatches = forward_span.edits;
                        origin.reverse_mismatches = reverse_span.edits;
                    }
                }
            }
        }
        for (const auto& item : origins) {
            panvar::FragmentPlacement placement;
            placement.forward_start = item.first.first;
            placement.insert_bp = item.first.second;
            placement.log_latent_mass = item.second.log_latent_mass;
            placement.latent_alignments = item.second.latent_alignments;
            placement.reverse_start = item.second.reverse_start;
            placement.forward_mismatches = item.second.forward_mismatches;
            placement.reverse_mismatches = item.second.reverse_mismatches;
            placement.forward_consumed_bp = item.second.forward_consumed_bp;
            placement.reverse_consumed_bp = item.second.reverse_consumed_bp;
            placement.orientation = arm == 0 ? panvar::MateOrientation::FirstForward : panvar::MateOrientation::SecondForward;
            out.push_back(placement);
        }
    }
    return out;
}

// A physical origin is identified by its two template ends and its orientation. The reverse start and
// the mismatch counts describe whichever latent alignment carried the most mass, which is a property
// of the explanation rather than of the place, so they are deliberately NOT part of the identity: two
// implementations that agree on the origins may legitimately pick different best alignments.
auto placement_key(const panvar::FragmentPlacement& p) {
    return std::make_tuple(static_cast<int>(p.orientation), p.forward_start, p.insert_bp);
}

struct EditRun {
    char op = 'M';
    std::size_t query_pos = 0;
    std::size_t target_pos = 0;
    std::size_t length = 0;
};

// True edlib global traceback, returned as runs of non-matching operations with positions in both
// sequences. 'I': bases in the query absent from the target; 'D': bases in the target absent from the
// query; 'X': substitutions.
std::vector<EditRun> edit_runs(const std::string& query, const std::string& target, long& distance) {
    EdlibAlignResult result = edlibAlign(query.data(), static_cast<int>(query.size()), target.data(),
                                         static_cast<int>(target.size()),
                                         edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
    if (result.status != EDLIB_STATUS_OK) throw std::runtime_error("edlib traceback failed");
    distance = result.editDistance;
    std::vector<EditRun> runs;
    std::size_t q = 0;
    std::size_t t = 0;
    for (int i = 0; i < result.alignmentLength; ++i) {
        const unsigned char op = result.alignment[i];
        char code = op == EDLIB_EDOP_MATCH ? 'M' : op == EDLIB_EDOP_INSERT ? 'I' : op == EDLIB_EDOP_DELETE ? 'D' : 'X';
        if (code != 'M') {
            if (!runs.empty() && runs.back().op == code &&
                runs.back().query_pos + (code == 'D' ? 0 : runs.back().length) == q &&
                runs.back().target_pos + (code == 'I' ? 0 : runs.back().length) == t) {
                ++runs.back().length;
            } else {
                runs.push_back({code, q, t, 1});
            }
        }
        if (op != EDLIB_EDOP_DELETE) ++q;
        if (op != EDLIB_EDOP_INSERT) ++t;
    }
    edlibFreeAlignResult(result);
    return runs;
}

int compare(const Args& args) {
    const std::string out = args.required("out");
    const std::string reference = args.required("reference");
    const std::string block_id = args.required("block");
    const std::size_t band = std::stoul(args.get("band", "6"));
    // The staged sorted graphs are an implementation detail of loading, not a result, and they are
    // 3.3 MB each. tied-consensus removes its own; compare did not, so every invocation left them
    // beside its output, and pointing --out at a results directory filled it with exactly the
    // redundant material such a directory exists to exclude.
    //
    // A guard rather than a call at the end of the function: every argument refusal below throws
    // AFTER the model graph has been staged, so a success-path cleanup still leaks one file per
    // rejected invocation -- measured, two leaked graphs from the two refusal cases alone. The
    // error_code overload is deliberate: a destructor that throws while an exception is unwinding
    // calls std::terminate, and unwinding is the case this guard exists for.
    struct StagedGraphs {
        std::string prefix;
        ~StagedGraphs() {
            std::error_code ignored;
            std::filesystem::remove(prefix + ".model.sorted.gfa", ignored);
            std::filesystem::remove(prefix + ".truth.sorted.gfa", ignored);
        }
    } staged_graphs{out};
    const LoadedModel loaded = load_model(args.required("model-gfa"), reference, out + ".model.sorted.gfa", 50);
    const panvar::PreparedGenotypeModel& model = *loaded.model;
    // --run-b is optional: with --alleles-b the second pair is the first one with a single block
    // replaced, which is the only way to attribute a score difference to one block.
    const std::string run_b_path = args.get("run-b", "");
    for (const std::string& run : {std::string("run-a"), std::string("run-b")}) {
        const std::string path = run == "run-a" ? args.required("run-a") : run_b_path;
        if (path.empty()) continue;
        if (read_key_values(path + ".mosaic_status.tsv").at("model_fingerprint") != loaded.fingerprint) {
            throw std::runtime_error(run + " was not scored on the rebuilt model");
        }
    }
    auto a = winner_paths(args.required("run-a"));
    auto b = run_b_path.empty() ? a : winner_paths(run_b_path);
    std::size_t block = model.blocks().size();
    for (std::size_t i = 0; i < model.blocks().size(); ++i) {
        if (model.blocks()[i].block_id == block_id) block = i;
    }
    if (block == model.blocks().size()) throw std::runtime_error("unknown block " + block_id);

    // Bounded single-block override. Comparing two runs' winners does NOT isolate a block: the
    // HG04036 10x-versus-30x comparison differed in seven of eleven, so attributing its 2126-nat
    // fragment-term gap to the sv:7 substitution alone was never established by that artifact.
    // Replacing only this block's two alleles on a common base pair leaves every other block
    // identical, so differing_blocks contains exactly one entry and the attribution is sound.
    const auto override_block = [&](const char* option,
                                    std::pair<std::vector<std::size_t>, std::vector<std::size_t>>& pair) {
        const std::string spec = args.get(option, "");
        if (spec.empty()) return std::string("none");
        const std::size_t comma = spec.find(',');
        if (comma == std::string::npos) {
            throw std::runtime_error(std::string(option) + " expects <allele1>,<allele2>");
        }
        const std::size_t first = std::stoul(spec.substr(0, comma));
        const std::size_t second = std::stoul(spec.substr(comma + 1));
        const std::size_t alleles = model.blocks()[block].alleles.size();
        if (first >= alleles || second >= alleles) {
            throw std::runtime_error(std::string(option) + " names an allele outside block " + block_id +
                                     ", which has " + std::to_string(alleles) + " alleles");
        }
        pair.first[block] = first;
        pair.second[block] = second;
        return spec;
    };
    const std::string overrode_a = override_block("alleles-a", a);
    const std::string overrode_b = override_block("alleles-b", b);

    std::ofstream summary(out + ".summary.tsv");
    summary << std::setprecision(12) << "key\tvalue\n";
    std::vector<std::string> differing;
    for (std::size_t i = 0; i < model.blocks().size(); ++i) {
        if (std::minmax(a.first[i], a.second[i]) != std::minmax(b.first[i], b.second[i])) {
            differing.push_back(model.blocks()[i].block_id);
        }
    }
    summary << "differing_blocks\t";
    for (std::size_t i = 0; i < differing.size(); ++i) summary << (i ? "," : "") << differing[i];
    summary << '\n';

    const std::vector<panvar::MosaicPath> paths = {panvar::MosaicPath{a.first}, panvar::MosaicPath{a.second},
                                                   panvar::MosaicPath{b.first}, panvar::MosaicPath{b.second}};
    std::vector<panvar::HaplotypeSequence> spelled;
    std::vector<std::string> sequences;
    for (const auto& path : paths) {
        spelled.push_back(panvar::spell_mosaic(model, path));
        sequences.push_back(spelled.back().bases);
    }
    const panvar::InternedHaplotypes all = panvar::intern_haplotypes(paths, sequences);
    const panvar::HaplotypeCandidate cand_a{all.index_of[0], all.index_of[1]};
    const panvar::HaplotypeCandidate cand_b{all.index_of[2], all.index_of[3]};

    // The prior is part of the score, and it differs between candidates because the Li-Stephens
    // switching path differs. A comparison reporting only fragment and exposure is not the model's
    // preference, which is the same mistake as reporting a fragment term alone. The pair prior is
    // taken from the same function the command uses, so the log(2) for two distinct sequences is
    // applied identically rather than reconstructed here.
    //
    // A haploid prior of -inf means no allele path spells that haplotype: blockwise spellability
    // asked exactly the way the model asks it, rather than as a separate reimplemented check.
    const double switch_probability = std::stod(args.get("switch-probability", "0.01"));
    const std::vector<double> switches_for_prior(model.blocks().size() - 1, switch_probability);
    const double prior_a =
        panvar::log_unordered_sequence_pair_prior(model, sequences[0], sequences[1], switches_for_prior);
    const double prior_b =
        panvar::log_unordered_sequence_pair_prior(model, sequences[2], sequences[3], switches_for_prior);
    double haploid_prior[4];
    for (std::size_t i = 0; i < 4; ++i) {
        haploid_prior[i] = panvar::log_haploid_sequence_prior(model, sequences[i], switches_for_prior);
    }
    const panvar::LibraryModel library = read_library_file(args.required("library"));
    const std::vector<panvar::ReadPair> fragments = read_pairs(args.required("reads1"), args.required("reads2"));
    panvar::FragmentScoringOptions options;
    options.mismatch_band = band;
    options.minimum_piece_bp = 12;

    // ---- 2. candidate-set independence ------------------------------------------------------------------
    const panvar::InternedHaplotypes only_a = panvar::intern_haplotypes({paths[0], paths[1]}, {sequences[0], sequences[1]});
    const panvar::InternedHaplotypes only_b = panvar::intern_haplotypes({paths[2], paths[3]}, {sequences[2], sequences[3]});
    const auto alone_a = panvar::score_diplotypes(only_a, {{only_a.index_of[0], only_a.index_of[1]}}, fragments, library, options);
    const auto alone_b = panvar::score_diplotypes(only_b, {{only_b.index_of[0], only_b.index_of[1]}}, fragments, library, options);
    panvar::FragmentScoringOptions detailed = options;
    detailed.record_origin_summary = true;
    const auto together = panvar::score_diplotypes(all, {cand_a, cand_b}, fragments, library, detailed);
    const auto same = [](const panvar::DiplotypeScore& x, const panvar::DiplotypeScore& y) {
        return x.fragment_term == y.fragment_term && x.exposure_term == y.exposure_term &&
               x.fragments_with_mass == y.fragments_with_mass;
    };
    summary << "unique_haplotypes\t" << all.sequences.size() << "\ncandidate_set_independent\t"
            << (same(alone_a.scores[0], together.scores[0]) && same(alone_b.scores[0], together.scores[1]) ? "true" : "false")
            << "\nfragment_term_a\t" << together.scores[0].fragment_term << "\nfragment_term_b\t"
            << together.scores[1].fragment_term << "\nfragments_with_mass_a\t" << together.scores[0].fragments_with_mass
            << "\nfragments_with_mass_b\t" << together.scores[1].fragments_with_mass << "\nexposure_term_a\t"
            << together.scores[0].exposure_term << "\nexposure_term_b\t" << together.scores[1].exposure_term << '\n';
    // The COMPLETE likelihood, not its observed-fragment summand. The exposure term is
    // candidate-dependent whenever two diplotypes differ in length, so a comparison reported as a
    // fragment-term gap alone can invert the model's actual preference -- which is exactly how the
    // HG04036 calibration was misread as a scorer inconsistency.
    summary << "log_likelihood_a\t" << together.scores[0].log_likelihood << "\nlog_likelihood_b\t"
            << together.scores[1].log_likelihood << "\ntotal_gap_b_minus_a\t"
            << (together.scores[1].log_likelihood - together.scores[0].log_likelihood)
            << "\nexposure_gap_b_minus_a\t"
            << (together.scores[1].exposure_term - together.scores[0].exposure_term)
            << "\nalleles_a_override\t" << overrode_a << "\nalleles_b_override\t" << overrode_b
            << "\nrun_b\t" << (run_b_path.empty() ? "none" : run_b_path) << '\n';
    // DiplotypeScore::log_likelihood carries fragment + exposure; the prior is added at certification.
    // The COMPLETE score is therefore log_likelihood + pair prior, and that is what decides a ranking.
    summary << std::setprecision(12) << "switch_probability\t" << switch_probability
            << "\nprior_a\t" << prior_a << "\nprior_b\t" << prior_b
            << "\nprior_gap_b_minus_a\t" << (prior_b - prior_a)
            << "\ncomplete_score_a\t" << (together.scores[0].log_likelihood + prior_a)
            << "\ncomplete_score_b\t" << (together.scores[1].log_likelihood + prior_b)
            << "\ncomplete_score_gap_b_minus_a\t"
            << ((together.scores[1].log_likelihood + prior_b) -
                (together.scores[0].log_likelihood + prior_a))
            << "\nb_outscores_a\t"
            << (((together.scores[1].log_likelihood + prior_b) >
                 (together.scores[0].log_likelihood + prior_a))
                    ? "true"
                    : "false")
            << '\n';
    // Blockwise spellability, asked as the model asks it: -inf means no allele path spells that
    // haplotype, so the panel cannot represent it at all.
    summary << "hap_a1_prior\t" << haploid_prior[0] << "\nhap_a2_prior\t" << haploid_prior[1]
            << "\nhap_b1_prior\t" << haploid_prior[2] << "\nhap_b2_prior\t" << haploid_prior[3]
            << "\na_spellable\t"
            << ((std::isfinite(haploid_prior[0]) && std::isfinite(haploid_prior[1])) ? "true" : "false")
            << "\nb_spellable\t"
            << ((std::isfinite(haploid_prior[2]) && std::isfinite(haploid_prior[3])) ? "true" : "false")
            << '\n';

    // ---- 3. per-fragment decomposition ------------------------------------------------------------------
    const std::size_t u = together.unique_haplotypes;
    struct Row {
        std::size_t fragment;
        double mass_a, mass_b, background, contribution_a, contribution_b, delta;
    };
    std::vector<Row> rows;
    double delta_sum = 0.0;
    double magnitude = 0.0;
    for (std::size_t f = 0; f < fragments.size(); ++f) {
        Row row;
        row.fragment = f;
        row.mass_a = panvar::log_add_exp(together.log_mass[f * u + cand_a.first], together.log_mass[f * u + cand_a.second]);
        row.mass_b = panvar::log_add_exp(together.log_mass[f * u + cand_b.first], together.log_mass[f * u + cand_b.second]);
        row.background = panvar::log_background_probability(fragments[f]);
        row.contribution_a = panvar::log_fragment_contribution(row.mass_a, row.background, library);
        row.contribution_b = panvar::log_fragment_contribution(row.mass_b, row.background, library);
        row.delta = row.contribution_b - row.contribution_a;
        delta_sum += row.delta;
        magnitude += std::abs(row.contribution_a) + std::abs(row.contribution_b);
        rows.push_back(row);
    }
    const double gap = together.scores[1].fragment_term - together.scores[0].fragment_term;
    const double tolerance = static_cast<double>(fragments.size() + 16) * std::numeric_limits<double>::epsilon() * magnitude;
    summary << "fragment_term_gap_b_minus_a\t" << gap << "\nsum_of_fragment_deltas\t" << delta_sum
            << "\ndecomposition_residual\t" << delta_sum - gap << "\ndecomposition_tolerance\t" << tolerance
            << "\ndecomposition_exact\t" << (std::abs(delta_sum - gap) <= tolerance ? "true" : "false") << '\n';

    std::sort(rows.begin(), rows.end(), [](const Row& x, const Row& y) { return x.delta > y.delta; });
    const auto best_of = [&](std::size_t f, const panvar::HaplotypeCandidate& c) {
        const auto& s1 = together.origin_summary[f * u + c.first];
        const auto& s2 = together.origin_summary[f * u + c.second];
        return s1.best_log_term >= s2.best_log_term ? std::make_pair(c.first, s1) : std::make_pair(c.second, s2);
    };
    const auto describe = [](const panvar::FragmentPlacement& p) {
        return std::string(p.orientation == panvar::MateOrientation::FirstForward ? "F1" : "F2") + ":" +
               std::to_string(p.forward_start) + "-" + std::to_string(p.reverse_start) + ":ins" +
               std::to_string(p.insert_bp) + ":mm" + std::to_string(p.forward_mismatches) + "," +
               std::to_string(p.reverse_mismatches);
    };
    {
        std::ofstream table(out + ".fragments.tsv");
        table << std::setprecision(10)
              << "fragment\tread\tmass_a\tmass_b\tbackground\tcontribution_a\tcontribution_b\tdelta\t"
                 "origins_a1\torigins_a2\torigins_b1\torigins_b2\tbest_a\tbest_b\tsignal_a\tsignal_b\n";
        for (const Row& row : rows) {
            const auto ba = best_of(row.fragment, cand_a);
            const auto bb = best_of(row.fragment, cand_b);
            table << row.fragment << '\t' << fragments[row.fragment].first.name << '\t' << row.mass_a << '\t'
                  << row.mass_b << '\t' << row.background << '\t' << row.contribution_a << '\t' << row.contribution_b
                  << '\t' << row.delta << '\t' << together.origin_summary[row.fragment * u + cand_a.first].origins << '\t'
                  << together.origin_summary[row.fragment * u + cand_a.second].origins << '\t'
                  << together.origin_summary[row.fragment * u + cand_b.first].origins << '\t'
                  << together.origin_summary[row.fragment * u + cand_b.second].origins << '\t'
                  << (ba.second.origins ? "h" + std::to_string(ba.first) + ":" + describe(ba.second.best) : "none") << '\t'
                  << (bb.second.origins ? "h" + std::to_string(bb.first) + ":" + describe(bb.second.best) : "none") << '\t'
                  << (std::isfinite(row.mass_a) ? "true" : "false") << '\t'
                  << (std::isfinite(row.mass_b) ? "true" : "false") << '\n';
        }
    }

    // ---- 4. signal-only fragments and the 90% set --------------------------------------------------------
    double only_b_delta = 0.0;
    double only_a_delta = 0.0;
    std::size_t only_b_count = 0;
    std::size_t only_a_count = 0;
    for (const Row& row : rows) {
        if (std::isfinite(row.mass_b) && !std::isfinite(row.mass_a)) { ++only_b_count; only_b_delta += row.delta; }
        if (std::isfinite(row.mass_a) && !std::isfinite(row.mass_b)) { ++only_a_count; only_a_delta += row.delta; }
    }
    std::size_t ninety = 0;
    double cumulative = 0.0;
    while (ninety < rows.size() && cumulative < 0.9 * gap) cumulative += rows[ninety++].delta;
    summary << "fragments_signal_only_under_b\t" << only_b_count << "\ndelta_from_signal_only_under_b\t" << only_b_delta
            << "\nshare_of_gap_signal_only_under_b\t" << only_b_delta / gap << "\nfragments_signal_only_under_a\t"
            << only_a_count << "\ndelta_from_signal_only_under_a\t" << only_a_delta << "\nsmallest_set_explaining_90pct\t"
            << ninety << "\nits_delta\t" << cumulative << '\n';
    std::size_t positive = 0;
    std::size_t negative = 0;
    for (const Row& row : rows) {
        if (row.delta > 1e-9) ++positive;
        if (row.delta < -1e-9) ++negative;
    }
    summary << "fragments_favouring_b\t" << positive << "\nfragments_favouring_a\t" << negative << '\n';

    // ---- 5. true edlib traceback at the block ------------------------------------------------------------
    std::size_t allele_a = 0;
    std::size_t allele_b = 0;
    std::size_t hap_a = 0;   // index into sequences of A's haplotype carrying the differing allele
    std::size_t hap_b = 0;
    {
        std::multiset<std::size_t> in_a = {a.first[block], a.second[block]};
        std::multiset<std::size_t> in_b = {b.first[block], b.second[block]};
        for (std::size_t x : {a.first[block], a.second[block]}) {
            if (in_b.count(x) < in_a.count(x)) allele_a = x;
        }
        for (std::size_t x : {b.first[block], b.second[block]}) {
            if (in_a.count(x) < in_b.count(x)) allele_b = x;
        }
        hap_a = a.first[block] == allele_a ? 0 : 1;
        hap_b = b.first[block] == allele_b ? 2 : 3;
    }
    const std::string& seq_a = model.blocks()[block].alleles[allele_a].sequence;
    const std::string& seq_b = model.blocks()[block].alleles[allele_b].sequence;
    const SortedGraph truth_graph = load_truth_in_model_frame(
        args.required("truth-gfa"), reference, loaded.sorted, out + ".truth.sorted.gfa");
    std::vector<std::string> truth_blocks;
    for (const panvar::PathRecord* path : donor_paths(truth_graph.graph, args.required("donor"))) {
        const auto projection = panvar::project_path_onto_locus(truth_graph.graph, loaded.bubbles, loaded.built.report, path->name);
        if (!projection.has_value() || !projection->blocks[block].has_value()) throw std::runtime_error("truth has no frame at the block");
        truth_blocks.push_back(*projection->blocks[block]);
    }

    // Every allele in this block against each truth haplotype, so the availability-floor alleles are
    // read off the data instead of guessed. Without this a --pair comparison has to be told which
    // alleles to use, and the obvious way to obtain them -- reading a run's winner -- can only ever
    // name alleles the caller already chose.
    {
        std::ofstream table(out + ".alleles.tsv");
        table << "allele\tbp";
        for (std::size_t t = 0; t < truth_blocks.size(); ++t) table << "\tdistance_to_truth" << t + 1;
        table << '\n';
        std::vector<long> floor_distance(truth_blocks.size(), -1);
        std::vector<std::size_t> floor_allele(truth_blocks.size(), 0);
        for (std::size_t i = 0; i < model.blocks()[block].alleles.size(); ++i) {
            const std::string& sequence = model.blocks()[block].alleles[i].sequence;
            table << i << '\t' << sequence.size();
            for (std::size_t t = 0; t < truth_blocks.size(); ++t) {
                const long d = edit_distance(truth_blocks[t], sequence);
                table << '\t' << d;
                if (floor_distance[t] < 0 || d < floor_distance[t]) {
                    floor_distance[t] = d;
                    floor_allele[t] = i;
                }
            }
            table << '\n';
        }
        summary << "block_alleles\t" << model.blocks()[block].alleles.size() << '\n';
        for (std::size_t t = 0; t < truth_blocks.size(); ++t) {
            summary << "floor_allele_truth" << t + 1 << '\t' << floor_allele[t]
                    << "\nfloor_distance_truth" << t + 1 << '\t' << floor_distance[t] << '\n';
        }
    }
    std::size_t truth_index = 0;
    long best_distance = -1;
    for (std::size_t t = 0; t < truth_blocks.size(); ++t) {
        const long d = edit_distance(truth_blocks[t], seq_a);
        if (best_distance < 0 || d < best_distance) { best_distance = d; truth_index = t; }
    }
    const std::string& truth = truth_blocks[truth_index];
    long d_truth_a = 0;
    long d_truth_b = 0;
    long d_a_b = 0;
    const auto runs_truth_a = edit_runs(truth, seq_a, d_truth_a);
    const auto runs_truth_b = edit_runs(truth, seq_b, d_truth_b);
    const auto runs_a_b = edit_runs(seq_a, seq_b, d_a_b);
    summary << "allele_a\t" << allele_a << "\nallele_b\t" << allele_b << "\nallele_a_bp\t" << seq_a.size()
            << "\nallele_b_bp\t" << seq_b.size() << "\ntruth_bp\t" << truth.size() << "\nedlib_truth_to_a\t" << d_truth_a
            << "\nedlib_truth_to_b\t" << d_truth_b << "\nedlib_a_to_b\t" << d_a_b << '\n';
    {
        std::ofstream table(out + ".alignments.tsv");
        table << "pair\top\tquery_pos\ttarget_pos\tlength\n";
        for (const auto& item : {std::make_pair("truth_to_a", &runs_truth_a), std::make_pair("truth_to_b", &runs_truth_b),
                                 std::make_pair("a_to_b", &runs_a_b)}) {
            for (const EditRun& run : *item.second) {
                table << item.first << '\t' << run.op << '\t' << run.query_pos << '\t' << run.target_pos << '\t' << run.length << '\n';
            }
        }
    }

    // ---- map the dominant fragments onto the differing positions -------------------------------------------
    const std::size_t dominant_count = std::min<std::size_t>(ninety, 200);
    const auto block_begin_a = spelled[hap_a].block_intervals[block].begin;
    const auto block_begin_b = spelled[hap_b].block_intervals[block].begin;
    const auto a_to_b_positions = [&](bool in_a_coordinates) {
        std::vector<std::pair<std::size_t, std::size_t>> spans;   // [begin, end) of each differing run
        for (const EditRun& run : runs_a_b) {
            const std::size_t begin = in_a_coordinates ? run.query_pos : run.target_pos;
            const std::size_t length = in_a_coordinates ? (run.op == 'D' ? 0 : run.length) : (run.op == 'I' ? 0 : run.length);
            spans.emplace_back(begin, begin + std::max<std::size_t>(length, 1));
        }
        return spans;
    };
    const auto spans_a = a_to_b_positions(true);
    const auto spans_b = a_to_b_positions(false);
    const auto overlap = [](const std::vector<std::pair<std::size_t, std::size_t>>& spans, long begin, long end) {
        std::size_t hits = 0;
        for (const auto& s : spans) {
            if (static_cast<long>(s.first) < end && static_cast<long>(s.second) > begin) ++hits;
        }
        return hits;
    };
    {
        std::ofstream table(out + ".dominant.tsv");
        table << "rank\tfragment\tread\tdelta\tsignal_a\tsignal_b\tfootprints_on_a_block\tdiffs_a\tfootprints_on_b_block\tdiffs_b\n";
        for (std::size_t i = 0; i < dominant_count; ++i) {
            const Row& row = rows[i];
            const auto& sa = together.origin_summary[row.fragment * u + all.index_of[hap_a]];
            const auto& sb = together.origin_summary[row.fragment * u + all.index_of[hap_b]];
            const std::size_t len1 = fragments[row.fragment].first.bases.size();
            const auto footprint = [&](const panvar::OriginSummary& s, std::size_t begin,
                                       const std::vector<std::pair<std::size_t, std::size_t>>& spans, std::size_t& diffs) {
                if (s.origins == 0) return std::string("none");
                const long f0 = static_cast<long>(s.best.forward_start) - static_cast<long>(begin);
                const long r0 = static_cast<long>(s.best.reverse_start) - static_cast<long>(begin);
                diffs = overlap(spans, f0, f0 + static_cast<long>(len1)) + overlap(spans, r0, r0 + static_cast<long>(len1));
                return std::to_string(f0) + "," + std::to_string(r0);
            };
            std::size_t diffs_a = 0;
            std::size_t diffs_b = 0;
            const std::string fa = footprint(sa, block_begin_a, spans_a, diffs_a);
            const std::string fb = footprint(sb, block_begin_b, spans_b, diffs_b);
            table << i + 1 << '\t' << row.fragment << '\t' << fragments[row.fragment].first.name << '\t' << row.delta << '\t'
                  << (std::isfinite(row.mass_a) ? "true" : "false") << '\t' << (std::isfinite(row.mass_b) ? "true" : "false")
                  << '\t' << fa << '\t' << diffs_a << '\t' << fb << '\t' << diffs_b << '\n';
        }
    }

    // ---- 6. accelerated origins against an index-free scan, for the dominant fragments -----------------------
    const std::size_t checked = std::min<std::size_t>(dominant_count, 40);
    std::vector<panvar::ReadPair> subset;
    std::vector<std::size_t> subset_ids;
    for (std::size_t i = 0; i < checked; ++i) {
        subset.push_back(fragments[rows[i].fragment]);
        subset_ids.push_back(i);
    }
    panvar::FragmentScoringOptions recording = options;
    recording.record_origins_for_fragments = subset_ids;
    const auto recorded = panvar::score_diplotypes(all, {cand_a, cand_b}, subset, library, recording);
    const panvar::DiscreteInsertDistribution insert = panvar::build_library_insert_distribution(library);
    std::size_t identity_mismatch = 0;
    std::size_t mass_mismatch = 0;
    std::size_t origins_compared = 0;
    for (std::size_t i = 0; i < checked; ++i) {
        const panvar::InsertSupport support = panvar::fragment_insert_support(subset[i], library);
        for (std::size_t h = 0; h < all.sequences.size(); ++h) {
            auto scanned = scan_origins(all.sequences[h], subset[i], support, library, band);
            const auto found = recorded.recorded_origins.find({i, h});
            auto produced = found == recorded.recorded_origins.end() ? std::vector<panvar::FragmentPlacement>{} : found->second;
            const auto order = [](const panvar::FragmentPlacement& x, const panvar::FragmentPlacement& y) {
                return placement_key(x) < placement_key(y);
            };
            std::sort(scanned.begin(), scanned.end(), order);
            std::sort(produced.begin(), produced.end(), order);
            bool identical = scanned.size() == produced.size();
            for (std::size_t k = 0; identical && k < scanned.size(); ++k) {
                identical = placement_key(scanned[k]) == placement_key(produced[k]);
            }
            if (!identical) ++identity_mismatch;
            origins_compared += scanned.size();
            const double scanned_mass = panvar::log_placement_mass(scanned, subset[i], insert, library);
            const double production_mass = recorded.mass_of(i, h);
            const double full_mass = together.mass_of(rows[i].fragment, h);
            const bool masses_agree =
                (!std::isfinite(scanned_mass) && !std::isfinite(production_mass) && !std::isfinite(full_mass)) ||
                (std::abs(scanned_mass - production_mass) <= 1e-9 && std::abs(production_mass - full_mass) <= 1e-9);
            if (!masses_agree) ++mass_mismatch;
        }
    }
    summary << "origin_check_fragments\t" << checked << "\norigin_check_origins\t" << origins_compared
            << "\norigin_identity_mismatches\t" << identity_mismatch << "\norigin_mass_mismatches\t" << mass_mismatch << '\n';

    // ---- 7. calibration: reads from each allele's diplotype ---------------------------------------------------
    const double depth = std::stod(args.get("depth", "30"));
    for (const auto& item : {std::make_pair(std::string("from_a"), std::vector<std::string>{sequences[0], sequences[1]}),
                             std::make_pair(std::string("from_b"), std::vector<std::string>{sequences[2], sequences[3]})}) {
        const auto reads = simulate_from(item.second, library, depth, fragments.front().first.bases.size(),
                                         item.first == "from_a" ? 811 : 812);
        const auto scored = panvar::score_diplotypes(all, {cand_a, cand_b}, reads, library, options);
        const double fragment_gap =
            scored.scores[1].fragment_term - scored.scores[0].fragment_term;
        const double exposure_gap =
            scored.scores[1].exposure_term - scored.scores[0].exposure_term;
        const double total_gap =
            scored.scores[1].log_likelihood - scored.scores[0].log_likelihood;
        // Generator preference is a statement about the declared point-process likelihood, not
        // merely its observed-fragment summand.  The exposure term is candidate-dependent whenever
        // the two diplotypes have different lengths; omitting it made reads generated from the
        // shorter HG04036 candidate appear to prefer the longer candidate even though the complete
        // likelihood favoured their own generator by roughly 1,800 nats at 30x.
        const bool expected = item.first == "from_a" ? total_gap < 0.0 : total_gap > 0.0;
        summary << "calibration_" << item.first << "_fragments\t" << reads.size() << "\ncalibration_" << item.first
                << "_fragment_gap_b_minus_a\t" << fragment_gap << "\ncalibration_" << item.first
                << "_exposure_gap_b_minus_a\t" << exposure_gap << "\ncalibration_" << item.first
                << "_total_gap_b_minus_a\t" << total_gap << "\ncalibration_" << item.first << "_prefers_generator\t"
                << (expected ? "true" : "false") << '\n';
    }

    // ---- 8. mismatch band, diagnostic only -----------------------------------------------------------------------
    for (const std::size_t trial : {3, 4, 5, 6, 7}) {
        panvar::FragmentScoringOptions banded = options;
        banded.mismatch_band = trial;
        const auto scored = panvar::score_diplotypes(all, {cand_a, cand_b}, fragments, library, banded);
        summary << "band_" << trial << "_gap_b_minus_a\t" << scored.scores[1].fragment_term - scored.scores[0].fragment_term
                << "\nband_" << trial << "_fragments_with_mass_a\t" << scored.scores[0].fragments_with_mass
                << "\nband_" << trial << "_fragments_with_mass_b\t" << scored.scores[1].fragments_with_mass << '\n';
    }
    // The staged graphs are removed by StagedGraphs above, on every exit path rather than this one.
    std::cout << "compared " << block_id << ": gap " << gap << ", decomposition residual " << delta_sum - gap << '\n';
    return 0;
}

} // namespace

// ------------------------------------------------------------------------------------------------
// preflight: what one locus would cost, before any donor is simulated or run
// ------------------------------------------------------------------------------------------------

// Reports a locus's structure and the work a leave-zero-out run would do on it, with every predicted count
// placed beside the declared cap it would hit. The runtime estimate is an extrapolation from the C4 coverage
// experiment, whose rescoring time was 4.32 s + 1.793 ms per fragment at 8,523 candidates over ~130 unique
// haplotypes (twelve runs, residuals under 3.8 s). One locus cannot separate the per-haplotype mass term from
// the per-candidate-pair term, so the estimate is bracketed: the low end attributes C4's slope entirely to the
// pair term, the high end entirely to the haplotype-mass term. Neither is a measurement of this locus.
// Every panel template respelled from the chain and compared against the graph's own path bytes, without
// stopping at the first disagreement. The round-trip gate throws on one template and says nothing about the
// rest or about where the bases went; this reports the shape of every failure: which end, how much, and
// whether the prefix agrees byte for byte (sequence lost from the tiling) or diverges inside (a wrong walk).
int respell(const Args& args) {
    const std::string reference = args.required("reference");
    const std::size_t min_variant_bp = std::stoul(args.get("min-variant-bp", "50"));
    const std::string stage =
        (std::filesystem::temp_directory_path() / ("genotype_probe_respell_" + std::to_string(::getpid()) + ".sorted.gfa")).string();
    LoadedModel loaded = load_model(args.required("gfa"), reference, stage, min_variant_bp);
    std::filesystem::remove(stage);

    const panvar::PreparedGenotypeModel& model = *loaded.model;
    const panvar::PanelAlleleMatrix& panel = model.panel();
    std::map<std::string, const std::string*> authoritative_by_name;
    for (const panvar::AuthoritativePanelSequence& sequence : loaded.built.authoritative) {
        authoritative_by_name.emplace(sequence.template_name, &sequence.source_frame_bases);
    }

    std::ostringstream out;
    out << "template\torientation\tspelled_bp\tauthoritative_bp\tdelta\tfirst_diff_source\tshape\n";
    std::size_t failing = 0;
    std::size_t first_failing = panel.template_count();
    for (std::size_t t = 0; t < panel.template_count(); ++t) {
        const std::string& name = panel.template_names[t];
        const auto found = authoritative_by_name.find(name);
        if (found == authoritative_by_name.end()) continue;
        const panvar::HaplotypeSequence spelled = panvar::spell_panel_template(model, t);
        const panvar::ChainOrientation orientation = model.template_source_orientations()[t];
        const std::string source_frame = panvar::sequence_in_source_frame(spelled.bases, orientation);
        const std::string& authoritative = *found->second;
        if (source_frame == authoritative) continue;
        if (first_failing == panel.template_count()) first_failing = t;
        ++failing;

        const std::size_t common = std::min(source_frame.size(), authoritative.size());
        std::size_t at = 0;
        while (at < common && source_frame[at] == authoritative[at]) ++at;
        // A pure truncation means every spelled byte is right and the chain simply stops: bases the blocks
        // never claimed. An internal difference means a block walked somewhere the path did not.
        const std::string shape = at == common
            ? (source_frame.size() < authoritative.size() ? "truncated_in_source_frame"
                                                          : "overlong_in_source_frame")
            : "internal_difference";
        out << name << '\t'
            << (orientation == panvar::ChainOrientation::Forward ? "forward" : "reverse") << '\t'
            << source_frame.size() << '\t' << authoritative.size() << '\t'
            << static_cast<long long>(source_frame.size()) - static_cast<long long>(authoritative.size())
            << '\t' << at << '\t' << shape << '\n';
    }
    out << "#templates\t" << panel.template_count() << "\n#failing\t" << failing << '\n';

    // The per-block breakdown of one failure. In chain frame a reverse template's source-frame tail is its
    // head, so naming the chain block that abuts the loss is what points at the primitive responsible.
    if (first_failing < panel.template_count()) {
        const panvar::HaplotypeSequence spelled = panvar::spell_panel_template(model, first_failing);
        out << "#blocks_of\t" << panel.template_names[first_failing] << '\n'
            << "block\tkind\tbp\tchain_begin\tchain_end\n";
        for (std::size_t bi = 0; bi < model.blocks().size(); ++bi) {
            const panvar::SequenceInterval interval = spelled.block_intervals[bi];
            out << model.blocks()[bi].block_id << '\t'
                << (model.blocks()[bi].kind == panvar::GenotypeBlockKind::Variant ? "variant" : "context")
                << '\t' << interval.size() << '\t' << interval.begin << '\t' << interval.end << '\n';
        }
    }

    std::cout << out.str();
    if (!args.get("out", "").empty()) {
        std::ofstream file(args.get("out", ""));
        file << out.str();
    }
    return 0;
}

// Per (block, template) coverage audit. Built from the locus input alone, never from a prepared
// model, so a locus whose model is refused is still auditable. "Bypass" is not one state: this
// separates a bracketed deletion from absent coverage at a terminal anchor.
int coverage_audit(const Args& args) {
    const std::string reference = args.required("reference");
    const std::size_t min_variant_bp = std::stoul(args.get("min-variant-bp", "50"));
    const std::string stage =
        (std::filesystem::temp_directory_path() / ("genotype_probe_coverage_" + std::to_string(::getpid()) + ".sorted.gfa")).string();
    const SortedGraph sorted = load_sorted(args.required("gfa"), reference, stage);
    std::filesystem::remove(stage);

    panvar::BubbleCallOptions bubble_options;
    bubble_options.reference_path = sorted.resolved;
    bubble_options.quiet = true;
    bubble_options.min_variant_bp = min_variant_bp;
    bubble_options.snarl_pairs_override =
        panvar::find_top_level_snarls_cactus(panvar::snarl_input_from_model(sorted.gfa));
    bubble_options.snarl_source_supplied = true;
    panvar::LocusInputOptions locus_options;
    locus_options.locus_id = sorted.resolved;
    locus_options.reference_path_name = sorted.resolved;
    const panvar::LocusInputResult built = panvar::build_locus_input(
        sorted.graph, panvar::call_bubbles(sorted.graph, bubble_options), locus_options);
    const panvar::LocusInputReport& report = built.report;
    const std::size_t framed = report.forward_templates + report.reverse_templates;

    std::map<std::string, std::map<std::string, std::size_t>> by_block;
    std::ostringstream out;
    out << "block\ttemplate\tclassification\n";
    for (const panvar::TemplateBlockCoverage& row : report.non_traversing_coverage) {
        const char* label = panvar::block_coverage_label(row.coverage, row.uncovered_on_left);
        out << row.block_id << '\t' << row.template_name << '\t' << label << '\n';
        ++by_block[row.block_id][label];
    }

    // The bounding nodes are what a template must actually reach for a block to be covered at all,
    // so they are what settles whether a block is uncovered or merely bypassed.
    out << "#block\tkind\tsource_node\tsink_node\ttraversing\tempty_between_boundaries\tbridged_bypass\tuncovered_left\tuncovered_right\n";
    for (const panvar::LocusBlockReport& block : report.blocks) {
        const auto found = by_block.find(block.block_id);
        const auto count = [&](const char* label) -> std::size_t {
            if (found == by_block.end()) return 0;
            const auto at = found->second.find(label);
            return at == found->second.end() ? 0 : at->second;
        };
        const std::size_t empty = count("empty_between_boundaries");
        const std::size_t bridged = count("bridged_bypass");
        const std::size_t left = count("uncovered_left");
        const std::size_t right = count("uncovered_right");
        out << block.block_id << '\t'
            << (block.kind == panvar::GenotypeBlockKind::Variant ? "variant" : "context") << '\t'
            << (block.source_node.empty() ? "-" : block.source_node) << '\t'
            << (block.sink_node.empty() ? "-" : block.sink_node) << '\t'
            << framed - empty - bridged - left - right << '\t' << empty << '\t' << bridged << '\t'
            << left << '\t' << right << '\n';
    }
    out << "#templates_framed\t" << framed << "\n#templates_refused\t" << report.refused_templates.size()
        << "\n#blocks\t" << report.blocks.size() << "\n#representation_dependent_blocks\t";
    for (std::size_t i = 0; i < report.representation_dependent_blocks.size(); ++i) {
        out << (i ? "," : "") << report.representation_dependent_blocks[i];
    }
    if (report.representation_dependent_blocks.empty()) out << "none";
    out << '\n';

    std::cout << out.str();
    if (!args.get("out", "").empty()) {
        std::ofstream file(args.get("out", ""));
        file << out.str();
    }
    return 0;
}

int preflight(const Args& args) {
    const std::string reference = args.required("reference");
    const double depth = std::stod(args.get("depth", "30"));
    const std::size_t read_length = std::stoul(args.get("read-length", "100"));
    const std::size_t mosaics = std::stoul(args.get("max-mosaic-candidates", "8"));
    const std::size_t audit = std::stoul(args.get("phase-audit-candidates", "0"));
    const std::size_t state_limit = std::stoul(args.get("state-limit", "1024"));
    const std::size_t min_variant_bp = std::stoul(args.get("min-variant-bp", "50"));
    const std::uint64_t max_mass_cells = std::stoull(args.get("max-mass-cells", "200000000"));
    const std::uint64_t max_panel_pair_evaluations = std::stoull(args.get("max-panel-pair-evaluations", "20000000000"));
    const std::uint64_t max_beam_extensions = std::stoull(args.get("max-beam-extensions", "2000000000"));
    const std::string stage =
        (std::filesystem::temp_directory_path() / ("genotype_probe_preflight_" + std::to_string(::getpid()) + ".sorted.gfa")).string();
    // Curation is part of the MODEL, not a reporting detail: curating two templates out changes the
    // panel, the diplotype count and the fingerprint. Without it, preflight could only ever describe the
    // uncurated graph, which is why the curated staged model could not be gated before a call existed.
    const std::string curation_path = args.get("panel-curation", "");
    const std::vector<std::string> curated =
        curation_path.empty() ? std::vector<std::string>{} : read_panel_curation(curation_path);
    LoadedModel loaded;
    std::string model_refusal;
    try {
        loaded = load_model(args.required("gfa"), reference, stage, min_variant_bp, curated);
    } catch (const std::exception& error) {
        model_refusal = error.what();
    }
    std::filesystem::remove(stage);
    if (!model_refusal.empty()) {
        // The locus cannot be modelled as it stands. Report what the structural pass did see -- partial frames
        // are the usual cause -- rather than only the exception.
        const SortedGraph sorted = load_sorted(args.required("gfa"), reference, stage + ".refused");
        std::filesystem::remove(stage + ".refused");
        panvar::BubbleCallOptions bubble_options;
        bubble_options.reference_path = sorted.resolved;
        bubble_options.quiet = true;
        bubble_options.min_variant_bp = min_variant_bp;
        bubble_options.snarl_pairs_override =
            panvar::find_top_level_snarls_cactus(panvar::snarl_input_from_model(sorted.gfa));
        bubble_options.snarl_source_supplied = true;
        panvar::LocusInputOptions locus_options;
        locus_options.locus_id = sorted.resolved;
        locus_options.reference_path_name = sorted.resolved;
        const panvar::LocusInputResult built =
            panvar::build_locus_input(sorted.graph, panvar::call_bubbles(sorted.graph, bubble_options), locus_options);
        std::ostringstream refused;
        refused << "key\tvalue\nlocus\t" << sorted.resolved << "\nreference\t" << reference
                << "\nmodel_refused\t" << model_refusal << "\ntemplates_in_graph\t"
                << built.report.forward_templates + built.report.reverse_templates
                << "\npartial_frames_refused\t" << built.report.refused_templates.size() << "\nblocks\t"
                << built.report.variant_blocks + built.report.context_blocks << "\nvariant_blocks\t"
                << built.report.variant_blocks << "\ncontext_blocks\t" << built.report.context_blocks
                << "\ncanonical_alleles\t" << built.report.total_canonical_alleles
                << "\nrepresentation_dependent_blocks\t" << built.report.representation_dependent_blocks.size()
                // A locus whose model cannot be built runs in no arm. Saying so explicitly beats omitting
                // the field, which a reader can only interpret as "not measured".
                << "\nruns_as_is\tfalse\nmodel_runnable\tfalse\npanel_arm_runnable\tfalse"
                   "\nmosaic_arm_runnable\tfalse\nphase_audit_runnable\tfalse\n";
        for (std::size_t i = 0; i < built.report.refused_templates.size() && i < 10; ++i) {
            refused << "refused_template_" << i << '\t' << built.report.refused_templates[i] << '\n';
        }
        std::cout << refused.str();
        if (!args.get("out", "").empty()) {
            std::ofstream file(args.get("out", ""));
            file << refused.str();
        }
        return 0;
    }
    const panvar::PreparedGenotypeModel& model = *loaded.model;
    const panvar::LocusInputReport& report = loaded.built.report;
    const panvar::PanelAlleleMatrix& panel = model.panel();

    // The two gates a real run applies after preparing the model. A preflight that skipped them called a locus
    // affordable that cannot build its model at all: cyp2d6 passes preparation and then fails the round trip.
    std::string round_trip = "ok";
    try {
        panvar::verify_complete_panel_round_trip(model, loaded.built.authoritative);
    } catch (const std::exception& error) {
        round_trip = error.what();
    }
    std::string representation = "ok";
    try {
        panvar::require_representation_independent(report);
    } catch (const std::exception& error) {
        representation = error.what();
    }

    // Every panel haplotype, spelled once: distinct sequences are what the scorer actually holds.
    std::set<std::string> distinct;
    std::size_t total_bp = 0;
    std::size_t min_bp = std::numeric_limits<std::size_t>::max();
    std::size_t max_bp = 0;
    for (std::size_t row = 0; row < panel.template_count(); ++row) {
        panvar::MosaicPath path;
        path.alleles.assign(panel.template_alleles[row].begin(), panel.template_alleles[row].end());
        const std::string bases = panvar::spell_mosaic(model, path).bases;
        distinct.insert(bases);
        total_bp += bases.size();
        min_bp = std::min(min_bp, bases.size());
        max_bp = std::max(max_bp, bases.size());
    }
    const double mean_bp = panel.template_count() == 0 ? 0.0
                                                       : static_cast<double>(total_bp) / static_cast<double>(panel.template_count());
    const std::size_t panel_pairs = panvar::panel_diplotypes(panel).size();
    const panvar::PreparedMarkerPanel markers =
        panvar::build_marker_panel(model, std::stoul(args.get("kmer", "31")), std::stoul(args.get("syncmer", "15")));
    std::size_t blocks_without_markers = 0;
    for (const panvar::BlockMarkerInventory& block : markers.blocks) {
        if (block.marker_slots.empty()) ++blocks_without_markers;
    }
    std::size_t bypass_blocks = 0;
    for (const panvar::LocusBlockReport& block : report.blocks) {
        if (block.bypass_templates > 0) ++bypass_blocks;
    }

    // One fragment covers two mates of read_length on one chromosome copy, and both copies are sequenced.
    const double fragments = mean_bp * depth / (2.0 * static_cast<double>(read_length));
    // Production scores the LARGER of the primary mosaic set and the phase audit, so whenever the audit
    // exceeds K it is the audit that sizes the scoring pass. Estimating from K alone understated the
    // haplotypes, the mass matrix and the runtime of every run that carries an audit.
    const std::size_t scored_mosaics = std::max(mosaics, audit);
    // The audit enlarges the CANDIDATE set, not the haplotype pool: candidates are diplotypes drawn
    // from haplotypes the panel and the beam already supply, so unique haplotypes do not scale with it.
    // Measured on c4 at 30x, the scored haplotypes are identical at audit 0 and audit 256 (130 in the
    // panel arm, 140 at K = 8), while the run takes 1.6x as long -- so the audit's cost is its own pass
    // and is priced below, not folded into the mass matrix.
    const double unique_haplotypes = static_cast<double>(distinct.size() + 2 * mosaics);
    const double mass_cells = unique_haplotypes * fragments;
    // The primary set only: the audit's candidates are priced by their own term below, and counting
    // them here as well double-counts them. Measured, the two costs are separable and additive.
    const double pair_evaluations = static_cast<double>(panel_pairs + mosaics) * fragments;

    // The same integer, overflow-checked planner production uses, not a second formula in floating point.
    std::uint64_t beam_extensions = 0;
    std::string beam_overflow = "no";
    try {
        beam_extensions =
            panvar::extension_upper_bound(panel.template_count(), model.blocks().size(), state_limit);
    } catch (const std::exception& error) {
        beam_overflow = error.what();
        beam_extensions = std::numeric_limits<std::uint64_t>::max();
    }
    // C4's measured slope, and the two ways of attributing it. Re-fitted from runs that hold everything
    // but one variable fixed. The previous constants came from the round 8 cohort, where EVERY run
    // carried the 256-candidate phase audit, so the audit's cost sat inside the base slope; adding an
    // audit term on top of that fit double-counted it and predicted 93-99 s for a 46.5 s run.
    // Base, phase audit off, one donor at three coverages: 7.571 s at 2,073 fragments, 13.790 s at
    // 10,540 and 29.065 s at 31,358 -- linear to within 5 ms.
    const double c4_fragments = 31358.0;
    const double c4_slope_seconds = 0.00073387 * c4_fragments;   // per-run cost above the fixed term
    const double by_pairs = c4_slope_seconds * (pair_evaluations / (8523.0 * c4_fragments));
    const double by_mass = c4_slope_seconds * (mass_cells * mean_bp) / (130.0 * c4_fragments * 226000.0);
    // Re-derived with c4's own prior cost (below) removed, so the two do not double-count: c4's
    // audit-off points become 6.427 / 12.646 / 27.921 s, which refit to this intercept and the slope above.
    const double fixed_seconds = 4.91;   // model, markers, proposal, spelling and the rescoring intercept
    // The exact-path prior is computed once per candidate diplotype and does NOT scale with coverage, so
    // neither term above can see it. On c4 it is 1.1 s and invisible; on a 108,000-diplotype panel it is
    // a third to a half of the whole run, which is why constants fitted on c4 alone under-predicted acot
    // by 2.3x. Measured at 2x on three loci, panel-only: acot 23.8 s over 2.07M candidate-blocks, gstm1
    // 27.6 s over 2.17M, ankrd36c 28.2 s over 2.27M -- 11.5, 12.7 and 12.5 us each, flat within
    // single-run noise, and confirmed by sampling the running process (log_unordered_diplotype_prior).
    const double prior_seconds =
        12.2e-6 * static_cast<double>(panel_pairs + mosaics) * static_cast<double>(model.blocks().size());

    // Peak memory is dominated by the spelled panel, not by the mass matrix: every template is held as
    // bases, so it scales with templates x length and not with coverage at all. Pricing only the mass
    // matrix was ~180x low on these loci. Measured at 2x, panel-only: 1,539 MB over 111.6 Mbp of panel
    // (acot), 1,331 over 105.4 (gstm1), 2,584 over 225.6 (ankrd36c) -- 13.8, 12.6 and 11.5 bytes per
    // panel base, flat within 9% across a 2.1x spread in template length. The mass matrix is added on
    // top of it because that is the part which grows with depth.
    const double panel_memory_MB =
        12.6 * static_cast<double>(panel.template_count()) * mean_bp / 1048576.0;
    // The phase audit is its own pass. It enlarges the CANDIDATE set without enlarging the haplotype
    // pool -- unique haplotypes are 139-140 at every audit size measured -- so neither term above sees
    // it, and pricing a run from K alone understated every audited run by more than half.
    // Two sweeps, each holding the other variable fixed: at 31,358 fragments the cost is affine in the
    // audit size (29.065 s at 0, then 35.539, 40.253 and 46.475 at 64, 128, 256), and at audit 256 the
    // cost above that fixed part is proportional to fragments (4.310 s at 2,073, 8.366 at 10,540,
    // 17.410 at 31,358). Both are reproduced within 0.21 s by a fixed setup plus a term in A x
    // fragments, which is the form below.
    const double audit_seconds =
        audit == 0 ? 0.0
                   : 3.36 + 0.0558 * static_cast<double>(audit) * (fragments / c4_fragments);

    // Each arm is priced on its own candidate set, and reported separately: a beam refusal says nothing
    // about whether the panel arm can run, and the panel arm does not use the beam at all.
    const auto arm_fits = [&](std::size_t candidates) {
        const double cells = static_cast<double>(distinct.size() + 2 * candidates) * fragments;
        const double pairs = static_cast<double>(panel_pairs + candidates) * fragments;
        return cells <= static_cast<double>(max_mass_cells) &&
               pairs <= static_cast<double>(max_panel_pair_evaluations);
    };
    const bool beam_fits = beam_overflow == "no" && beam_extensions <= max_beam_extensions;
    const bool model_runnable = round_trip == "ok" && representation == "ok";
    const bool panel_arm_runnable = model_runnable && arm_fits(0);
    const bool mosaic_arm_runnable = model_runnable && beam_fits && arm_fits(mosaics);
    const bool phase_audit_runnable = model_runnable && beam_fits && arm_fits(audit);
    const bool within_caps = beam_fits && arm_fits(scored_mosaics);
    const bool runs_as_is = model_runnable && within_caps;

    std::ostringstream out;
    out << std::setprecision(10) << "key\tvalue\n"
        << "locus\t" << model.locus_id() << "\nreference\t" << reference
        << "\ntemplates\t" << panel.template_count() << "\ndistinct_template_sequences\t" << distinct.size()
        << "\npanel_diplotypes\t" << panel_pairs
        // Requested versus actually applied. A curated name that is not a graph path removes nothing,
        // and a gate comparing only the requested count would not notice.
        << "\npanel_curation\t" << (curation_path.empty() ? "none" : curation_path)
        << "\npanel_curation_requested\t" << curated.size()
        << "\ntemplates_curated_out\t"
        << [&] {
               std::size_t applied = 0;
               for (const std::string& name : curated) {
                   const bool in_graph =
                       std::any_of(loaded.sorted.graph.paths.begin(), loaded.sorted.graph.paths.end(),
                                   [&](const panvar::PathRecord& path) { return path.name == name; });
                   const bool in_panel =
                       std::find(panel.template_names.begin(), panel.template_names.end(), name) !=
                       panel.template_names.end();
                   if (in_graph && !in_panel) ++applied;
               }
               return applied;
           }()
        << "\ntemplate_bp_min_mean_max\t" << (panel.template_count() ? min_bp : 0) << "," << mean_bp << "," << max_bp
        << "\nmodel_fingerprint\t" << loaded.fingerprint
        << "\npanel_round_trip\t" << round_trip << "\nrepresentation_independent\t" << representation
        << "\nruns_as_is\t" << (runs_as_is ? "true" : "false")
        << "\npredicted_work_within_caps\t" << (within_caps ? "true" : "false")
        << "\nmodel_runnable\t" << (model_runnable ? "true" : "false")
        << "\npanel_arm_runnable\t" << (panel_arm_runnable ? "true" : "false")
        << "\nmosaic_arm_runnable\t" << (mosaic_arm_runnable ? "true" : "false")
        << "\nphase_audit_runnable\t" << (phase_audit_runnable ? "true" : "false")
        << "\nbeam_extension_overflow\t" << beam_overflow
        << "\nblocks\t" << model.blocks().size() << "\nvariant_blocks\t" << report.variant_blocks
        << "\ncontext_blocks\t" << report.context_blocks
        << "\nsource_alleles\t" << report.total_source_alleles
        << "\ncanonical_alleles\t" << report.total_canonical_alleles
        << "\npartial_frames_refused\t" << report.refused_templates.size()
        << "\nblocks_with_bypass_templates\t" << bypass_blocks
        << "\nrepresentation_dependent_blocks\t" << report.representation_dependent_blocks.size()
        << "\nmarker_universe\t" << markers.codes.size() << "\nmarker_anchors\t" << markers.anchor_slots.size()
        << "\nblocks_without_markers\t" << blocks_without_markers
        << "\nmarkers_dropped_multiblock\t" << markers.dropped_multiblock_markers
        << "\nalleles_scored_of_total\t" << markers.scored_alleles << "/" << markers.total_alleles
        << "\nassumed_depth\t" << depth << "\nassumed_read_length\t" << read_length
        << "\nmax_mosaic_candidates\t" << mosaics << "\nphase_audit_candidates\t" << audit
        << "\nscored_mosaic_candidates\t" << scored_mosaics << "\nstate_limit\t" << state_limit
        << "\npredicted_fragments\t" << fragments
        << "\npredicted_unique_haplotypes\t" << unique_haplotypes
        << "\npredicted_mass_cells\t" << mass_cells << "\nmass_cells_cap\t" << max_mass_cells
        << "\nmass_cells_fraction_of_cap\t" << mass_cells / static_cast<double>(max_mass_cells)
        << "\npredicted_panel_pair_evaluations\t" << pair_evaluations
        << "\npanel_pair_evaluations_cap\t" << max_panel_pair_evaluations
        << "\npanel_pair_evaluations_fraction_of_cap\t" << pair_evaluations / static_cast<double>(max_panel_pair_evaluations)
        << "\nbeam_extension_upper_bound\t" << beam_extensions << "\nbeam_extensions_cap\t" << max_beam_extensions
        << "\nbeam_extensions_fraction_of_cap\t" << beam_extensions / static_cast<double>(max_beam_extensions)
        << "\npredicted_mass_matrix_MB\t" << mass_cells * 8.0 / 1048576.0
        << "\npredicted_panel_memory_MB\t" << panel_memory_MB
        << "\npredicted_peak_rss_MB\t" << panel_memory_MB + mass_cells * 8.0 / 1048576.0
        << "\nestimated_audit_seconds\t" << audit_seconds
        << "\nestimated_prior_seconds\t" << prior_seconds
        << "\nestimated_seconds_low\t"
        << fixed_seconds + audit_seconds + prior_seconds + std::min(by_pairs, by_mass)
        << "\nestimated_seconds_high\t"
        << fixed_seconds + audit_seconds + prior_seconds + std::max(by_pairs, by_mass)
        << "\nestimate_basis\tFour terms, each measured with the others held fixed. (1) Base, phase audit "
           "off, on c4: 4.91 s + 0.734 ms per fragment at 8,523 candidates, ~130 haplotypes, 226 kb, scaled by "
           "candidate pairs at the low end and haplotypes x length at the high end; fitted to 7.571 / 13.790 / "
           "29.065 s at 2,073 / 10,540 / 31,358 fragments with c4's own prior cost removed so (3) does not "
           "double-count. (2) Phase audit: 3.36 s + 55.8 ms per candidate per c4-30x-fragment, over audit "
           "0/64/128/256 at one coverage and audit 256 at three, within 0.21 s. (3) Exact-path prior: 12.2 us "
           "per candidate-block, INDEPENDENT of coverage -- 1.1 s on c4 and a third to a half of the run on a "
           "108,000-diplotype panel; measured on acot/gstm1/ankrd36c at 2x as 11.5/12.7/12.5 us and confirmed "
           "by sampling the running process. (4) Memory: 12.6 bytes per panel base for the spelled panel "
           "(13.8/12.6/11.5 measured), plus the mass matrix, which is the part that grows with depth. "
           "Constants for (1) and (2) still come from c4 alone; (3) and (4) come from three large loci at 2x, "
           "and no term has been checked above 2x on a large panel.\n";
    std::cout << out.str();
    if (!args.get("out", "").empty()) {
        std::ofstream file(args.get("out", ""));
        file << out.str();
    }
    return 0;
}

// ------------------------------------------------------------------------------------------------
// anchors: a donor's invariant-anchor counts, for the depth-fit oracle gate
// ------------------------------------------------------------------------------------------------

// Writes the histogram of (haploid multiplicity, observed count) over the invariant anchors as C++ initialiser
// rows, so a unit test can refit the real anchors without the reads.
int anchors(const Args& args) {
    const std::string out = args.required("out");
    // The staged graph is scratch: it goes to the system temporary directory, never beside the output, which
    // for a test fixture is the source tree.
    const std::string stage =
        (std::filesystem::temp_directory_path() / ("genotype_probe_anchors_" + std::to_string(::getpid()) + ".sorted.gfa")).string();
    const LoadedModel loaded = load_model(args.required("gfa"), args.required("reference"), stage, 50);
    std::filesystem::remove(stage);
    const std::vector<panvar::ReadPair> fragments = read_pairs(args.required("reads1"), args.required("reads2"));
    const panvar::PreparedMarkerPanel markers = panvar::build_marker_panel(*loaded.model, 31, 15);
    const panvar::MarkerObservations observations = panvar::count_marker_observations(markers, fragments);
    std::map<std::pair<std::uint32_t, std::uint64_t>, std::uint64_t> histogram;
    for (std::size_t k = 0; k < markers.anchor_slots.size(); ++k) {
        ++histogram[{markers.anchor_haploid_multiplicity[k], observations.counts[markers.anchor_slots[k]]}];
    }
    std::ofstream file(out);
    file << "// " << args.required("label") << ": invariant-anchor histogram from genotype_probe anchors\n"
         << "// {haploid multiplicity, observed count, anchors}; " << markers.anchor_slots.size() << " anchors\n";
    for (const auto& item : histogram) {
        file << "{" << item.first.first << ", " << item.first.second << ", " << item.second << "},\n";
    }
    std::cout << "anchors " << args.required("label") << ": " << markers.anchor_slots.size() << " anchors, "
              << histogram.size() << " histogram rows\n";
    return 0;
}

// ------------------------------------------------------------------------------------------------
// orientations: which templates the chain can splice, before any donor is drawn
// ------------------------------------------------------------------------------------------------

// mosaic-truth refuses to splice a reverse-oriented walk against a forward one, so a donor can source a
// synthetic mosaic only when BOTH of its panel rows are forward. That decides which donors a held-out
// draw may even consider, and it has to be known BEFORE the draw rather than discovered when staging
// throws -- a draw that silently skipped its refusals would no longer be the rule it claims to be.
//
// respell already computes this, but emits a row only for templates whose round trip FAILS, and on C4
// none do: it reports 131 templates, 0 failing, and not one orientation. So it cannot answer the
// question at all. This reports every template, and nothing else.
int orientations(const Args& args) {
    const std::string reference = args.required("reference");
    const std::size_t min_variant_bp = std::stoul(args.get("min-variant-bp", "50"));
    const std::string stage =
        (std::filesystem::temp_directory_path() /
         ("genotype_probe_orient_" + std::to_string(::getpid()) + ".sorted.gfa")).string();
    const LoadedModel loaded = load_model(args.required("gfa"), reference, stage, min_variant_bp);
    std::filesystem::remove(stage);
    const panvar::PanelAlleleMatrix& panel = loaded.model->panel();

    std::ostringstream out;
    out << "template\torientation\n";
    std::size_t forward = 0;
    for (std::size_t t = 0; t < panel.template_count(); ++t) {
        const bool is_forward =
            loaded.model->template_source_orientations()[t] == panvar::ChainOrientation::Forward;
        if (is_forward) ++forward;
        out << panel.template_names[t] << '\t' << (is_forward ? "forward" : "reverse") << '\n';
    }
    out << "#templates\t" << panel.template_count() << "\n#forward\t" << forward << "\n#reverse\t"
        << panel.template_count() - forward << '\n';
    if (args.get("out", "").empty()) {
        std::cout << out.str();
    } else {
        std::ofstream file(args.get("out", ""));
        file << out.str();
    }
    std::cout << "orientations: " << forward << " forward, " << panel.template_count() - forward
              << " reverse, of " << panel.template_count() << " templates\n";
    return 0;
}

// ------------------------------------------------------------------------------------------------
// panel-alleles: the panel allele matrix, so a draw can tell an off-panel plan from a collapsed one
// ------------------------------------------------------------------------------------------------

// A switch point next to a chain end takes only a block or two from donor A, and if donor B's row
// carries the same alleles there, the "mosaic" IS donor B's row: mosaic-truth then refuses it, rightly,
// as an answer reachable without recombining. Whether that happens is a property of the DATA, not of
// the donor names, so a draw over names alone cannot avoid it and will keep proposing plans that refuse.
//
// This dumps the one array the refusal is computed from -- panel().template_alleles -- so the draw can
// apply exactly the same test before it commits a plan, rather than discovering it at staging time.
// Frozen alongside the orientation census, for the same reason: a rule that needs the model to be
// rebuilt before it can be checked is not a rule anyone can check.
int panel_alleles(const Args& args) {
    const std::string reference = args.required("reference");
    const std::size_t min_variant_bp = std::stoul(args.get("min-variant-bp", "50"));
    const std::string stage =
        (std::filesystem::temp_directory_path() /
         ("genotype_probe_alleles_" + std::to_string(::getpid()) + ".sorted.gfa")).string();
    const LoadedModel loaded = load_model(args.required("gfa"), reference, stage, min_variant_bp);
    std::filesystem::remove(stage);
    const panvar::PanelAlleleMatrix& panel = loaded.model->panel();

    std::ostringstream out;
    out << "template\talleles\n";
    for (std::size_t t = 0; t < panel.template_count(); ++t) {
        out << panel.template_names[t] << '\t';
        for (std::size_t b = 0; b < panel.template_alleles[t].size(); ++b) {
            out << (b ? "," : "") << panel.template_alleles[t][b];
        }
        out << '\n';
    }
    out << "#templates\t" << panel.template_count() << "\n#blocks\t"
        << loaded.model->blocks().size() << '\n';
    if (args.get("out", "").empty()) {
        std::cout << out.str();
    } else {
        std::ofstream file(args.get("out", ""));
        file << out.str();
    }
    std::cout << "panel-alleles: " << panel.template_count() << " templates over "
              << loaded.model->blocks().size() << " blocks\n";
    return 0;
}

// ------------------------------------------------------------------------------------------------
// tied-consensus: what a refused run's tied winners agree on, and how far that agreement is from truth
// ------------------------------------------------------------------------------------------------

// A RESCORE_INCOMPLETE run emits status only -- no genotypes, no block genotypes, no FASTA -- and that
// is correct: a tie the scorer cannot separate is not a call. But the tied set itself is on disk
// (mosaic_scores.tsv flags it, mosaic_candidates.tsv carries the allele paths), so what those
// candidates agree on, and how far the agreement sits from the truth, is measurable after the fact.
//
// Everything this prints is a CONDITIONAL_CONSENSUS_DIAGNOSTIC. It is not a call, it changes no call,
// and it must never enter a strict call rate or a strict accuracy rate. The command is read-only: it
// reads a finished run and the model, and writes only its own tables.
int tied_consensus(const Args& args) {
    const std::string reference = args.required("reference");
    const std::string donor = args.required("donor");
    const std::string run = args.required("run");
    const std::string out = args.required("out");
    const std::string curation_path = args.get("panel-curation", "");
    const LoadedModel loaded = load_model(
        args.required("gfa"), reference, out + ".model.sorted.gfa", 50,
        curation_path.empty() ? std::vector<std::string>{} : read_panel_curation(curation_path));
    const panvar::PreparedGenotypeModel& model = *loaded.model;
    const std::size_t blocks = model.blocks().size();

    // Truth projected onto the model's block frame, exactly as evaluate does it -- which includes
    // placing the truth graph in the model's node frame before projecting anything.
    const SortedGraph truth_graph = load_truth_in_model_frame(
        args.required("truth-gfa"), reference, loaded.sorted, out + ".truth.sorted.gfa");
    std::vector<panvar::PathProjection> truth;
    for (const panvar::PathRecord* path : donor_paths(truth_graph.graph, donor)) {
        const auto projection = panvar::project_path_onto_locus(truth_graph.graph, loaded.bubbles,
                                                                loaded.built.report, path->name);
        if (!projection.has_value()) {
            throw std::runtime_error("truth path " + path->name + " does not span the block chain");
        }
        truth.push_back(*projection);
    }
    std::filesystem::remove(out + ".model.sorted.gfa");
    std::filesystem::remove(out + ".truth.sorted.gfa");

    const auto spelling = [&](std::size_t b, std::size_t allele) -> const std::string& {
        return model.blocks()[b].alleles[allele].sequence;
    };
    // The phase-invariant diploid distance at one block: the truth pair is unordered, so the predicted
    // pair is scored in whichever assignment fits better. -1 means the truth is undefined here, which
    // is reported as NA and never as zero.
    const auto block_distance = [&](std::size_t b, std::size_t a1, std::size_t a2) -> long {
        if (!truth[0].blocks[b].has_value() || !truth[1].blocks[b].has_value()) return -1;
        const std::string& t1 = *truth[0].blocks[b];
        const std::string& t2 = *truth[1].blocks[b];
        return std::min(edit_distance(t1, spelling(b, a1)) + edit_distance(t2, spelling(b, a2)),
                        edit_distance(t1, spelling(b, a2)) + edit_distance(t2, spelling(b, a1)));
    };

    std::set<std::string> winners;
    for (const auto& row : read_table(run + ".mosaic_scores.tsv")) {
        if (row.at("in_winner_set") == "true") winners.insert(row.at("representative_beam_candidate"));
    }
    struct Tied {
        std::string id;
        std::vector<std::size_t> hap1;
        std::vector<std::size_t> hap2;
    };
    std::vector<Tied> tied;
    for (const auto& row : read_table(run + ".mosaic_candidates.tsv")) {
        if (winners.count(row.at("candidate_id")) == 0) continue;
        tied.push_back({row.at("candidate_id"), parse_alleles(row.at("hap1_alleles")),
                        parse_alleles(row.at("hap2_alleles"))});
    }
    if (tied.empty()) throw std::runtime_error("no tied winner resolved from " + run);

    // Per candidate: distance to truth, and whether it is exactly right.
    std::ofstream candidates(out + ".candidates.tsv");
    candidates << "candidate_id\ttotal_distance\tmax_block_distance\tsv_subtotal\tcontext_subtotal\texact\n";
    std::vector<long> totals;
    std::size_t exact_candidates = 0;
    bool truth_pair_present = false;
    for (const Tied& candidate : tied) {
        long total = 0, worst = 0, sv = 0, context = 0;
        bool defined = true;
        for (std::size_t b = 0; b < blocks; ++b) {
            const long d = block_distance(b, candidate.hap1[b], candidate.hap2[b]);
            if (d < 0) { defined = false; continue; }
            total += d;
            worst = std::max(worst, d);
            (model.blocks()[b].kind == panvar::GenotypeBlockKind::Variant ? sv : context) += d;
        }
        if (defined && total == 0) { ++exact_candidates; truth_pair_present = true; }
        totals.push_back(total);
        candidates << candidate.id << '\t' << total << '\t' << worst << '\t' << sv << '\t' << context
                   << '\t' << (defined && total == 0 ? "true" : "false") << '\n';
    }
    std::sort(totals.begin(), totals.end());
    const long median = totals.size() % 2 == 1
                            ? totals[totals.size() / 2]
                            : (totals[totals.size() / 2 - 1] + totals[totals.size() / 2]) / 2;

    // Per block: do all tied winners predict the SAME unordered genotype, and how far is it from truth?
    std::ofstream per_block(out + ".blocks.tsv");
    per_block << "block_id\tkind\tunanimous\tconsensus_allele1\tconsensus_allele2\tdistance\texact\t"
                 "truth1_bp\ttruth2_bp\n";
    std::size_t unanimous_blocks = 0, exact_consensus_blocks = 0, undefined_blocks = 0;
    long consensus_total = 0, consensus_worst = 0, consensus_sv = 0, consensus_context = 0;
    for (std::size_t b = 0; b < blocks; ++b) {
        std::set<std::pair<std::size_t, std::size_t>> genotypes;
        for (const Tied& candidate : tied) {
            genotypes.insert(std::minmax(candidate.hap1[b], candidate.hap2[b]));
        }
        const bool unanimous = genotypes.size() == 1;
        const bool variant = model.blocks()[b].kind == panvar::GenotypeBlockKind::Variant;
        const std::string t1 = truth[0].blocks[b].has_value()
                                   ? std::to_string(truth[0].blocks[b]->size()) : std::string("NA");
        const std::string t2 = truth[1].blocks[b].has_value()
                                   ? std::to_string(truth[1].blocks[b]->size()) : std::string("NA");
        per_block << model.blocks()[b].block_id << '\t' << (variant ? "SV" : "context") << '\t'
                  << (unanimous ? "true" : "false");
        if (!unanimous) {
            per_block << "\tNA\tNA\tNA\tNA\t" << t1 << '\t' << t2 << '\n';
            continue;
        }
        ++unanimous_blocks;
        const std::size_t a1 = genotypes.begin()->first;
        const std::size_t a2 = genotypes.begin()->second;
        const long d = block_distance(b, a1, a2);
        per_block << '\t' << a1 << '\t' << a2 << '\t';
        if (d < 0) {
            ++undefined_blocks;
            per_block << "NA\tNA\t" << t1 << '\t' << t2 << '\n';
            continue;
        }
        consensus_total += d;
        consensus_worst = std::max(consensus_worst, d);
        (variant ? consensus_sv : consensus_context) += d;
        if (d == 0) ++exact_consensus_blocks;
        per_block << d << '\t' << (d == 0 ? "true" : "false") << '\t' << t1 << '\t' << t2 << '\n';
    }

    std::ofstream summary(out + ".summary.tsv");
    summary << "key\tvalue\nlabel\tCONDITIONAL_CONSENSUS_DIAGNOSTIC\nrun\t" << run << "\ndonor\t" << donor
            << "\nwinner_set_size\t" << tied.size()
            << "\ntrue_pair_in_tied_set\t" << (truth_pair_present ? "true" : "false")
            << "\nany_tied_candidate_exact\t" << (exact_candidates > 0 ? "true" : "false")
            << "\nexact_tied_candidates\t" << exact_candidates
            << "\nmin_whole_locus_distance\t" << totals.front()
            << "\nmedian_whole_locus_distance\t" << median
            << "\nmax_whole_locus_distance\t" << totals.back()
            << "\nblocks\t" << blocks
            << "\nunanimous_blocks\t" << unanimous_blocks
            << "\nexact_consensus_blocks\t" << exact_consensus_blocks
            << "\nconsensus_total_distance\t" << consensus_total
            << "\nconsensus_max_block_distance\t" << consensus_worst
            << "\nconsensus_sv_subtotal\t" << consensus_sv
            << "\nconsensus_context_subtotal\t" << consensus_context
            << "\nblocks_with_undefined_truth\t" << undefined_blocks << '\n';
    std::cout << "tied-consensus " << run << ": " << tied.size() << " tied winners, " << unanimous_blocks
              << "/" << blocks << " blocks unanimous, consensus distance " << consensus_total
              << " (DIAGNOSTIC, not a call)\n";
    return 0;
}

// ------------------------------------------------------------------------------------------------
// mosaic-truth: build a synthetic mosaic donor whose haplotypes are in no panel row
// ------------------------------------------------------------------------------------------------

// The existing off-panel fixture proves the machinery can recover a mosaic on a tiny locus. It does not
// show that proposal and scoring recover one on a real repetitive locus, because every large-locus
// experiment so far used donors whose haplotypes were intact panel rows: the correct answer needed zero
// switches, and mosaic candidates never changed a call.
//
// This builds the missing input. Two development donors, a per-block choice between them, and the
// result expressed TWICE by independent routes: as an allele path, which spell_mosaic turns into the
// truth sequence and the per-block truth alleles; and as a spliced node walk, which becomes a real GFA
// path reads can be simulated from. The two must agree byte for byte. If they do not, the splice is
// wrong, and the run refuses rather than scoring a truth nobody verified.
int mosaic_truth(const Args& args) {
    const std::string reference = args.required("reference");
    const std::string stage =
        (std::filesystem::temp_directory_path() /
         ("genotype_probe_mosaic_" + std::to_string(::getpid()) + ".sorted.gfa")).string();
    // The staged graph is kept, not deleted here: it is the base the mosaic paths are written onto.
    // Everything below -- the walks, the node sequences, the boundary nodes -- lives in SORTED node ids,
    // and reference sorting renumbers nodes, so writing those ids onto the original file would produce
    // paths whose ids name different nodes. Same id space in, same id space out.
    const LoadedModel loaded = load_model(args.required("gfa"), reference, stage,
                                          std::stoul(args.get("min-variant-bp", "50")));
    const panvar::PreparedGenotypeModel& model = *loaded.model;
    const std::size_t blocks = model.blocks().size();

    const std::string donor_a = args.required("donor-a");
    const std::string donor_b = args.required("donor-b");

    // Panel rows of each donor, in name order, so hap 1 and hap 2 are taken consistently.
    const auto rows_of = [&](const std::string& donor) {
        std::vector<std::size_t> rows;
        for (std::size_t r = 0; r < model.panel().template_count(); ++r) {
            if (model.panel().template_names[r].rfind(donor + "#", 0) == 0) rows.push_back(r);
        }
        if (rows.size() != 2) {
            throw std::runtime_error("donor " + donor + " must have exactly two panel rows, found " +
                                     std::to_string(rows.size()));
        }
        return rows;
    };
    const std::vector<std::size_t> rows_a = rows_of(donor_a);
    const std::vector<std::size_t> rows_b = rows_of(donor_b);

    // A reverse-oriented walk runs opposite to chain order, so splicing it against a forward one at a
    // shared boundary node would join two walks that traverse the locus in opposite directions. Refused
    // rather than handled: this experiment does not need it, and a silent wrong splice would look like
    // a caller failure.
    for (const std::size_t r : {rows_a[0], rows_a[1], rows_b[0], rows_b[1]}) {
        if (model.template_source_orientations()[r] != panvar::ChainOrientation::Forward) {
            throw std::runtime_error("mosaic donors must be forward-oriented; " +
                                     model.panel().template_names[r] + " is not");
        }
    }

    // Which of donor B's rows pairs with donor A's row h. Drawing this bit varies the donor/HAPLOTYPE
    // combination rather than only the switch points, so two plans over the same two donors need not
    // recombine the same two pairs of haplotypes.
    const bool swap_hap_b = args.get("swap-hap-b", "0") == "1";
    const auto row_of = [&](std::size_t h, int donor) {
        return donor == 0 ? rows_a[h] : rows_b[swap_hap_b ? 1 - h : h];
    };

    // Switch plan, PER HOMOLOGUE: which donor supplies each block of each haplotype. "every" alternates
    // at every block; a comma-separated list names the block indices at which the source flips; empty
    // leaves that homologue as donor A's row, unswitched. --switches sets both homologues together,
    // which is what the five development plans used and keeps them byte-identical; --switches-hap1 and
    // --switches-hap2 set them independently. The lockstep form cannot express a mosaic on one homologue
    // beside an intact one, which is a class the held-out protocol freezes, so it is not a convenience.
    const std::string spec = args.get("switches", "");
    const std::string spec_of_hap[2] = {args.get("switches-hap1", spec), args.get("switches-hap2", spec)};
    std::vector<std::vector<int>> donor_of_block(2, std::vector<int>(blocks, 0));
    for (std::size_t h = 0; h < 2; ++h) {
        const std::string& hap_spec = spec_of_hap[h];
        if (hap_spec == "every") {
            for (std::size_t b = 0; b < blocks; ++b) donor_of_block[h][b] = static_cast<int>(b % 2);
        } else if (!hap_spec.empty()) {
            std::set<std::size_t> flips;
            std::stringstream ss(hap_spec);
            std::string token;
            while (std::getline(ss, token, ',')) {
                if (!token.empty()) flips.insert(std::stoul(token));
            }
            if (flips.empty()) {
                throw std::runtime_error("a switch spec must name block indices, or 'every', or be empty");
            }
            if (*flips.rbegin() >= blocks) {
                throw std::runtime_error("switch spec names block " + std::to_string(*flips.rbegin()) +
                                         " but the chain has " + std::to_string(blocks) + " blocks");
            }
            int current = 0;
            for (std::size_t b = 0; b < blocks; ++b) {
                if (flips.count(b) != 0) current = 1 - current;
                donor_of_block[h][b] = current;
            }
        }
    }
    std::size_t switches_of_hap[2] = {0, 0};
    for (std::size_t h = 0; h < 2; ++h) {
        for (std::size_t b = 1; b < blocks; ++b) {
            switches_of_hap[h] += donor_of_block[h][b] != donor_of_block[h][b - 1];
        }
    }
    const std::size_t switches = switches_of_hap[0] + switches_of_hap[1];

    // How much of this diplotype the plan DECLARES to be off-panel, checked rather than assumed.
    // "both" is the five development plans' guarantee and stays the default: neither homologue is a
    // panel row. "one" is a switched homologue beside an intact one -- the diplotype is still off-panel,
    // because a panel diplotype is a pair of panel rows and one of these is not, so the correct answer
    // is still unreachable without recombining. "none" is the unchanged panel-haplotype control, which
    // is not a mosaic at all; it is declared, and then required to really be one, rather than being
    // quietly exempted from the check.
    const std::string require_off_panel = args.get("require-off-panel", "both");
    if (require_off_panel != "both" && require_off_panel != "one" && require_off_panel != "none") {
        throw std::runtime_error("--require-off-panel must be both, one or none");
    }
    if (switches == 0 && require_off_panel != "none") {
        throw std::runtime_error("the plan switches neither homologue: that is not a mosaic");
    }

    // Route one: the allele path, and the truth it defines.
    std::vector<panvar::MosaicPath> mosaic_paths(2);
    for (std::size_t h = 0; h < 2; ++h) {
        mosaic_paths[h].alleles.resize(blocks);
        for (std::size_t b = 0; b < blocks; ++b) {
            mosaic_paths[h].alleles[b] =
                model.panel().template_alleles[row_of(h, donor_of_block[h][b])][b];
        }
    }
    // Which homologues coincide with a panel row, against what the plan declared.
    std::size_t on_panel = 0;
    std::string on_panel_names;
    for (std::size_t h = 0; h < 2; ++h) {
        for (std::size_t r = 0; r < model.panel().template_count(); ++r) {
            if (model.panel().template_alleles[r] == mosaic_paths[h].alleles) {
                ++on_panel;
                if (!on_panel_names.empty()) on_panel_names += ",";
                on_panel_names += "hap" + std::to_string(h + 1) + "=" + model.panel().template_names[r];
                break;
            }
        }
    }
    if ((require_off_panel == "both" && on_panel != 0) || (require_off_panel == "one" && on_panel == 2)) {
        throw std::runtime_error("--require-off-panel " + require_off_panel + ", but " +
                                 std::to_string(on_panel) + " homologue(s) are exactly panel rows (" +
                                 on_panel_names + "): the correct answer is reachable without recombining");
    }
    if (require_off_panel == "none" && on_panel != 2) {
        throw std::runtime_error("--require-off-panel none declares an unchanged panel-haplotype control, "
                                 "but only " + std::to_string(on_panel) +
                                 " of 2 homologues are panel rows; that is not a control");
    }

    // Route two: splice the node walks, block by block, at the boundary nodes every template shares.
    const auto path_named = [&](const std::string& name) {
        for (const panvar::PathRecord& path : loaded.sorted.graph.paths) {
            if (path.name == name) return &path;
        }
        throw std::runtime_error("no graph path named " + name);
        return static_cast<const panvar::PathRecord*>(nullptr);
    };
    const auto index_of_node = [](const panvar::PathRecord& path, const std::string& node) {
        for (std::size_t i = 0; i < path.steps.size(); ++i) {
            if (path.steps[i].node_id == node) return i;
        }
        throw std::runtime_error("path " + path.name + " does not visit boundary node " + node);
        return std::size_t{0};
    };

    std::vector<std::vector<panvar::PathStep>> walks(2);
    for (std::size_t h = 0; h < 2; ++h) {
        for (std::size_t b = 0; b < blocks; ++b) {
            const std::size_t row = row_of(h, donor_of_block[h][b]);
            const panvar::PathRecord& path = *path_named(model.panel().template_names[row]);
            const std::string& source = loaded.built.report.blocks[b].source_node;
            const std::string& sink = loaded.built.report.blocks[b].sink_node;
            const std::size_t from = source.empty() ? 0 : index_of_node(path, source);
            const std::size_t to = sink.empty() ? path.steps.size() : index_of_node(path, sink);
            for (std::size_t i = from; i < to; ++i) walks[h].push_back(path.steps[i]);
        }
        // Every junction must be a real edge. A splice that is not is not a haplotype.
        for (std::size_t i = 1; i < walks[h].size(); ++i) {
            if (!loaded.sorted.graph.has_edge_between(walks[h][i - 1].node_id, walks[h][i].node_id)) {
                throw std::runtime_error("spliced walk leaves the graph between " +
                                         walks[h][i - 1].node_id + " and " + walks[h][i].node_id);
            }
        }
    }

    // The two routes must agree. This is the check that the splice reproduces the allele path, rather
    // than something plausible that merely looks like a haplotype.
    for (std::size_t h = 0; h < 2; ++h) {
        std::string walked;
        for (const panvar::PathStep& step : walks[h]) {
            const auto found = loaded.sorted.gfa.seq.find(step.node_id);
            if (found == loaded.sorted.gfa.seq.end()) {
                throw std::runtime_error("spliced walk names an unknown node " + step.node_id);
            }
            walked += found->second;
        }
        const std::string spelled = panvar::spell_mosaic(model, mosaic_paths[h]).bases;
        if (walked != spelled) {
            throw std::runtime_error(
                "mosaic haplotype " + std::to_string(h + 1) + ": the spliced walk spells " +
                std::to_string(walked.size()) + " bases but its allele path spells " +
                std::to_string(spelled.size()) + "; the splice does not reproduce the truth");
        }
    }

    // The staged graph: the original, plus the two mosaic paths, so simulate --donor MOSAIC can draw
    // reads from haplotypes that exist in no panel row.
    const std::string out_gfa = args.required("out-gfa");
    {
        // The SORTED graph, not the original: the walks above are in sorted node ids. Node 19 is a
        // different base in the two id spaces, so appending sorted-id paths to the original file
        // produces a path that leaves the graph at the first renumbered node -- which is exactly how
        // this failed the first time, with the model refusing "steps from 16+ to 19+ but the graph has
        // no such link". A later re-sort of this staged file is harmless, because a sort renames S, L
        // and P lines together.
        panvar::GzLineReader reader(stage);
        if (!reader.ok()) throw std::runtime_error("cannot reopen the sorted graph to stage the mosaic");
        std::ofstream out(out_gfa);
        std::string line;
        while (reader.getline(line)) out << line << '\n';
        for (std::size_t h = 0; h < 2; ++h) {
            out << "P\tMOSAIC#" << (h + 1) << "#synthetic:0-0\t";
            for (std::size_t i = 0; i < walks[h].size(); ++i) {
                out << (i ? "," : "") << walks[h][i].node_id << (walks[h][i].reverse ? '-' : '+');
            }
            out << "\t*\n";
        }
    }
    // Scratch, and only needed until the staged copy exists.
    std::filesystem::remove(stage);

    // The truth table. Source samples are reported as an EQUIVALENCE CLASS of every panel template
    // carrying that exact block sequence: several samples may carry an identical allele, which makes
    // the originating sample statistically unidentifiable and a "source sample" column a fiction.
    std::ostringstream truth;
    truth << "block_id\tkind\thap1_allele_bp\thap2_allele_bp\thap1_from\thap2_from\t"
             "hap1_carriers\thap2_carriers\n";
    for (std::size_t b = 0; b < blocks; ++b) {
        const panvar::PreparedGenotypeBlock& block = model.blocks()[b];
        truth << block.block_id << '\t'
              << (block.kind == panvar::GenotypeBlockKind::Variant ? "SV" : "context");
        for (std::size_t h = 0; h < 2; ++h) {
            truth << '\t' << block.alleles[mosaic_paths[h].alleles[b]].sequence.size();
        }
        for (std::size_t h = 0; h < 2; ++h) {
            truth << '\t' << (donor_of_block[h][b] == 0 ? donor_a : donor_b);
        }
        for (std::size_t h = 0; h < 2; ++h) {
            std::vector<std::string> carriers;
            for (std::size_t r = 0; r < model.panel().template_count(); ++r) {
                if (model.panel().template_alleles[r][b] == mosaic_paths[h].alleles[b]) {
                    carriers.push_back(model.panel().template_names[r]);
                }
            }
            truth << '\t' << carriers.size() << ":";
            for (std::size_t i = 0; i < carriers.size() && i < 4; ++i) {
                truth << (i ? "," : "") << carriers[i];
            }
            if (carriers.size() > 4) truth << ",...";
        }
        truth << '\n';
    }
    const std::string out_truth = args.required("out");
    { std::ofstream file(out_truth); file << truth.str(); }

    std::cout << "mosaic-truth " << donor_a << " x " << donor_b << ": " << switches_of_hap[0] << "+"
              << switches_of_hap[1] << " switch(es) over " << blocks << " blocks, " << on_panel
              << " homologue(s) on panel (declared " << require_off_panel
              << "); both routes agree; staged " << out_gfa << '\n';
    return 0;
}

// ------------------------------------------------------------------------------------------------
// repair-pilot: post-ranking block repair, scored in ONE pass over the union
//
// EXPERIMENTAL EVALUATION PATH, not the production caller. Production would need a persistent scoring
// context owning stable read forms, indices, haplotype buffers and cache across incremental scoring;
// this deliberately does not approximate that by keeping a pointer-based cache alive between calls.
// There is exactly one score_diplotypes invocation here, over the union of the original candidates and
// the repairs, which is why no cache has to survive anything.
//
// Seeds come from the SIGNATURE-BOUND campaign outputs rather than from a preliminary in-run ranking:
// the best panel candidate, and every member of the numerical winner-equivalence set. Reading them off
// disk is what dissolves the circularity -- seeds normally require a ranking, and a ranking requires
// the scores this pass produces.
//
// TRUTH IS NOT IN SCOPE during repair construction. The repairs are built from the model, the block
// proposal emissions and the campaign tables alone, in a scope entered before the truth graph is
// loaded, so "no truth field is reachable by the construction" is a property of the code's shape
// rather than a promise. Truth is loaded afterwards, and only to score the outcome.
int repair_pilot(const Args& args) {
    const std::string out = args.required("out");
    struct StagedGraphs {
        std::string prefix;
        ~StagedGraphs() {
            std::error_code ignored;
            std::filesystem::remove(prefix + ".model.sorted.gfa", ignored);
            std::filesystem::remove(prefix + ".truth.sorted.gfa", ignored);
        }
    } staged_graphs{out};
    const std::string reference = args.required("reference");
    const std::string donor = args.required("donor");
    const std::string run = args.required("run");
    const std::size_t local_pairs = std::stoul(args.get("local-pairs", "32"));
    const double switch_probability = std::stod(args.get("switch-probability", "0.01"));
    const std::uint64_t declared_cap = std::stoull(args.get("max-repair-candidates", "0"));

    const LoadedModel loaded =
        load_model(args.required("model-gfa"), reference, out + ".model.sorted.gfa", 50);
    const panvar::PreparedGenotypeModel& model = *loaded.model;
    const std::size_t blocks = model.blocks().size();
    const auto status = read_key_values(run + ".mosaic_status.tsv");
    if (status.at("model_fingerprint") != loaded.fingerprint) {
        throw std::runtime_error("the run was not scored on the rebuilt model");
    }

    const panvar::LibraryModel library = read_library_file(args.required("library"));
    const std::vector<panvar::ReadPair> fragments =
        read_pairs(args.required("reads1"), args.required("reads2"));

    // Block proposal emissions, rebuilt exactly as the run built them. The fitted background is not
    // declared in these libraries, so the fit runs; its rate and background are checked against the
    // run's recorded values, because emissions that differ would rank a different set of local pairs.
    const panvar::PreparedMarkerPanel markers = panvar::build_marker_panel(model, 31, 15);
    const panvar::MarkerObservations observations = panvar::count_marker_observations(markers, fragments);
    std::optional<double> declared_background;
    {
        std::ifstream file(args.required("library"));
        std::string key;
        double value = 0.0;
        for (std::string line; std::getline(file, line);) {
            std::istringstream fields(line);
            if (fields >> key >> value && key == "marker_background") declared_background = value;
        }
    }
    const panvar::MarkerDepthEstimate depth =
        panvar::estimate_marker_depth(markers, observations, declared_background);
    if (!depth.converged) throw std::runtime_error("the marker depth fit did not converge");
    const panvar::PreparedProposalEmissions emissions = panvar::build_proposal_emissions(
        markers, observations, panvar::proposal_marker_parameters(depth, 31, 15));

    std::ofstream summary(out + ".summary.tsv");
    summary << std::setprecision(12) << "key\tvalue\n";
    summary << "donor\t" << donor << "\nrun\t" << run << "\nlocal_pairs\t" << local_pairs << '\n';
    summary << "marker_depth_rate_per_copy\t" << depth.rate_per_copy << "\nmarker_background\t"
            << depth.background << "\nmarker_depth_termination\t" << depth.termination << '\n';
    summary << "recorded_marker_depth_rate_per_copy\t" << status.at("marker_depth_rate_per_copy")
            << "\nrecorded_marker_background\t" << status.at("marker_background") << '\n';

    // ---- the campaign's candidate set, read back ---------------------------------------------------
    struct Original {
        std::vector<std::size_t> first;
        std::vector<std::size_t> second;
        std::string source;
        // The row's own sequence-candidate id. The reproduction guard below compares each original
        // against ITS OWN recorded score: matching against the nearest score anywhere in the table
        // would pass whenever a candidate happened to coincide with some other candidate's value.
        std::string sequence_candidate_id;
    };
    std::vector<Original> originals;
    std::map<std::string, std::size_t> original_of_id;
    for (const auto& row : read_table(run + ".mosaic_candidates.tsv")) {
        Original entry;
        entry.first = parse_alleles(row.at("hap1_alleles"));
        entry.second = parse_alleles(row.at("hap2_alleles"));
        entry.source = row.count("source") != 0 ? row.at("source") : std::string("unknown");
        entry.sequence_candidate_id =
            row.count("sequence_candidate_id") != 0 ? row.at("sequence_candidate_id") : std::string();
        original_of_id.emplace(row.at("candidate_id"), originals.size());
        originals.push_back(std::move(entry));
    }
    if (originals.empty()) throw std::runtime_error("the run has no candidates");

    // Seeds: the best panel candidate, and every member of the winner-equivalence set.
    std::vector<std::size_t> seeds;
    std::set<std::size_t> seen_seed;
    std::string best_panel_id;
    double best_panel_score = 0.0;
    bool have_panel = false;
    std::map<std::string, double> recorded_score;
    for (const auto& row : read_table(run + ".mosaic_scores.tsv")) {
        const std::string rep = row.at("representative_beam_candidate");
        const double score = std::stod(row.at("score"));
        recorded_score.emplace(row.at("sequence_candidate_id"), score);
        const auto found = original_of_id.find(rep);
        if (found == original_of_id.end()) continue;
        if (originals[found->second].source == "panel" && (!have_panel || score > best_panel_score)) {
            have_panel = true;
            best_panel_score = score;
            best_panel_id = rep;
        }
        if (row.at("in_winner_set") == "true" && seen_seed.insert(found->second).second) {
            seeds.push_back(found->second);
        }
    }
    if (have_panel) {
        const std::size_t slot = original_of_id.at(best_panel_id);
        if (seen_seed.insert(slot).second) seeds.push_back(slot);
    }
    summary << "seeds\t" << seeds.size() << "\nbest_panel_candidate\t"
            << (have_panel ? best_panel_id : std::string("none")) << '\n';

    // ---- the declared bound, computed from the block shape rather than assumed ----------------------
    // Sum over blocks of min(local_pairs, unordered pairs available), times two homologue assignments,
    // times the seeds. A block with few alleles cannot supply local_pairs, so a flat seeds x blocks x
    // local_pairs x 2 would name a number no configuration can reach.
    std::uint64_t cap = 0;
    for (std::size_t b = 0; b < blocks; ++b) {
        const std::uint64_t a = model.blocks()[b].alleles.size();
        const std::uint64_t available = a * (a + 1) / 2;
        cap += std::min<std::uint64_t>(local_pairs, available) * 2;
    }
    cap *= seeds.size();
    summary << "repair_cap\t" << cap << "\ndeclared_cap\t" << declared_cap << '\n';
    if (declared_cap != 0 && cap > declared_cap) {
        summary << "refused\ttrue\nrefusal\tthe repair construction exceeds its declared bound\n";
        std::cout << "repair-pilot " << donor << ": REFUSED before allocation, " << cap
                  << " repair candidates exceeds the declared " << declared_cap << '\n';
        return 3;
    }
    // ---- repair construction. TRUTH IS NOT IN SCOPE: nothing below reads it, and the truth graph is
    // ---- not loaded until every repair has been built and scored.
    const std::vector<double> switches(blocks - 1, switch_probability);
    std::vector<panvar::MosaicPath> paths;
    std::vector<std::string> sequences;
    std::map<std::vector<std::size_t>, std::size_t> slot_of_path;
    const auto slot_for = [&](const std::vector<std::size_t>& alleles) {
        const auto found = slot_of_path.find(alleles);
        if (found != slot_of_path.end()) return found->second;
        const std::size_t slot = paths.size();
        slot_of_path.emplace(alleles, slot);
        paths.push_back(panvar::MosaicPath{alleles});
        sequences.push_back(panvar::spell_mosaic(model, panvar::MosaicPath{alleles}).bases);
        return slot;
    };
    struct Entry {
        std::size_t first_slot = 0;
        std::size_t second_slot = 0;
        std::size_t rank = 0;      // 0 for originals; the local-pair rank that admitted a repair
        bool repair = false;
        // Self-description, so a candidate table needs nothing reconstructed later. The previous
        // revision recorded only per-prefix outcomes, which made the R=16/R=32 comparison impossible
        // to answer from its own outputs.
        std::string seed_id;       // sequence-candidate id of the seed this repair came from
        std::size_t block = 0;     // the repaired block
        int assignment = 0;        // 0: (x,y) onto (hap1,hap2); 1: the swap
        std::size_t allele_x = 0;
        std::size_t allele_y = 0;
    };
    std::vector<Entry> entries;
    // Parallel to `entries` for the originals only: the sequence-candidate id each one was read from,
    // captured here because originals are deduplicated by allele-path pair on the way in, so indexing
    // `originals` afterwards would not line up.
    std::vector<std::string> original_entry_id;
    std::set<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> seen_path_pair;
    const auto key_of = [](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
        return a <= b ? std::make_pair(a, b) : std::make_pair(b, a);
    };
    for (const Original& entry : originals) {
        if (!seen_path_pair.insert(key_of(entry.first, entry.second)).second) continue;
        entries.push_back({slot_for(entry.first), slot_for(entry.second), 0, false, std::string(), 0, 0,
                           0, 0});
        original_entry_id.push_back(entry.sequence_candidate_id);
    }
    const std::size_t original_entries = entries.size();

    std::uint64_t generated = 0;
    std::uint64_t rejected_prior = 0;
    std::uint64_t duplicate = 0;
    for (const std::size_t seed : seeds) {
        const Original& base = originals[seed];
        for (std::size_t b = 0; b < blocks; ++b) {
            const panvar::BlockProposalEmissions& emission = emissions.blocks[b];
            const std::size_t alleles = emission.allele_count;
            std::vector<std::pair<double, std::pair<std::size_t, std::size_t>>> ranked;
            ranked.reserve(alleles * (alleles + 1) / 2);
            for (std::size_t x = 0; x < alleles; ++x) {
                for (std::size_t y = x; y < alleles; ++y) {
                    ranked.push_back({emission.at(x, y), {x, y}});
                }
            }
            // Deterministic: emission descending, then allele indices, so no tie depends on sort order.
            std::sort(ranked.begin(), ranked.end(),
                      [](const auto& lhs, const auto& rhs) {
                          if (lhs.first != rhs.first) return lhs.first > rhs.first;
                          return lhs.second < rhs.second;
                      });
            const std::size_t take = std::min(local_pairs, ranked.size());
            for (std::size_t r = 0; r < take; ++r) {
                const std::size_t x = ranked[r].second.first;
                const std::size_t y = ranked[r].second.second;
                for (int assignment = 0; assignment < 2; ++assignment) {
                    if (assignment == 1 && x == y) continue;   // the swap is the same diplotype
                    std::vector<std::size_t> first = base.first;
                    std::vector<std::size_t> second = base.second;
                    first[b] = assignment == 0 ? x : y;
                    second[b] = assignment == 0 ? y : x;
                    ++generated;
                    if (!seen_path_pair.insert(key_of(first, second)).second) {
                        ++duplicate;
                        continue;
                    }
                    const std::size_t s1 = slot_for(first);
                    const std::size_t s2 = slot_for(second);
                    // A finite switching prior is required: -inf means no allele path spells it.
                    const double prior = panvar::log_unordered_sequence_pair_prior(
                        model, sequences[s1], sequences[s2], switches);
                    if (!std::isfinite(prior)) {
                        ++rejected_prior;
                        continue;
                    }
                    entries.push_back({s1, s2, r + 1, true, base.sequence_candidate_id, b,
                                       assignment, x, y});
                }
            }
        }
    }
    const std::size_t repair_entries = entries.size() - original_entries;
    summary << "repairs_generated\t" << generated << "\nrepairs_duplicate\t" << duplicate
            << "\nrepairs_rejected_infinite_prior\t" << rejected_prior << "\nrepairs_kept\t"
            << repair_entries << "\noriginal_entries\t" << original_entries << '\n';

    // ---- ONE scoring pass over the union -----------------------------------------------------------
    const panvar::InternedHaplotypes interned = panvar::intern_haplotypes(paths, sequences);
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> scored_index;
    std::vector<panvar::HaplotypeCandidate> scored_pairs;
    std::vector<std::size_t> pair_of_entry(entries.size(), 0);
    for (std::size_t i = 0; i < entries.size(); ++i) {
        std::size_t x = interned.index_of[entries[i].first_slot];
        std::size_t y = interned.index_of[entries[i].second_slot];
        if (y < x) std::swap(x, y);
        const auto inserted = scored_index.emplace(std::make_pair(x, y), scored_pairs.size());
        if (inserted.second) scored_pairs.push_back({x, y});
        pair_of_entry[i] = inserted.first->second;
    }
    summary << "unique_haplotypes\t" << interned.sequences.size() << "\nscored_pairs\t"
            << scored_pairs.size() << '\n';

    panvar::FragmentScoringOptions scoring;
    scoring.threads = std::stoul(args.get("threads", "8"));
    scoring.mismatch_band = std::stoul(status.at("mismatch_band"));
    scoring.minimum_piece_bp = 12;
    scoring.max_piece_hits = std::stoull(status.at("max_piece_hits"));
    scoring.max_proposed_starts_per_haplotype = std::stoull(status.at("max_proposed_starts_per_haplotype"));
    scoring.max_physical_origins = std::stoull(status.at("max_physical_origins"));
    scoring.max_state_classes = std::stoull(status.at("max_state_classes"));
    scoring.max_mass_cells = std::stoull(status.at("max_mass_cells"));
    scoring.emission_cache_bytes = static_cast<std::size_t>(std::stoull(status.at("emission_cache_bytes")));
    const auto scoring_mark = std::chrono::steady_clock::now();
    const panvar::FragmentScoringResult scored =
        panvar::score_diplotypes(interned, scored_pairs, fragments, library, scoring);
    const double scoring_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - scoring_mark).count();
    summary << "scoring_seconds\t" << scoring_seconds << "\nemission_cache_dp_hits\t"
            << scored.emission_cache.dp_cache_hits << "\nemission_cache_dp_misses\t"
            << scored.emission_cache.dp_cache_misses << "\nemission_cache_executed_dps\t"
            << scored.emission_cache.executed_dps << "\nemission_cache_entries\t"
            << scored.emission_cache.entries << '\n';

    // Per-entry complete score. certify_scored_ranking adds the prior itself, so the reproduction
    // check below compares log_likelihood + prior -- adding the certified score would count it twice.
    // Once per SCORED PAIR, not once per entry: thousands of entries collapse onto far fewer pairs,
    // and this prior is a forward recursion over the panel, not a lookup.
    std::vector<double> prior_of_pair(scored_pairs.size(), 0.0);
    for (std::size_t slot = 0; slot < scored_pairs.size(); ++slot) {
        prior_of_pair[slot] = panvar::log_unordered_sequence_pair_prior(
            model, interned.sequences[scored_pairs[slot].first],
            interned.sequences[scored_pairs[slot].second], switches);
    }

    // ---- the guard: every ORIGINAL score must reproduce before any repair is interpreted -----------
    std::size_t compared = 0;
    std::size_t mismatched = 0;
    double worst = 0.0;
    std::size_t unkeyed = 0;
    for (std::size_t i = 0; i < original_entries; ++i) {
        const std::size_t slot = pair_of_entry[i];
        const double mine = scored.scores[slot].log_likelihood + prior_of_pair[slot];
        // Each original against ITS OWN recorded row, keyed by sequence-candidate id. Matching against
        // the nearest value anywhere in the table would pass whenever a candidate coincided with some
        // other candidate's score, and would cost 8,384^2 comparisons to do it.
        const auto found = recorded_score.find(original_entry_id[i]);
        ++compared;
        if (found == recorded_score.end()) {
            ++unkeyed;
            ++mismatched;
            continue;
        }
        const double delta = std::fabs(found->second - mine);
        if (delta > 1e-6) {
            ++mismatched;
            worst = std::max(worst, delta);
        }
    }
    summary << "originals_unkeyed\t" << unkeyed << '\n';
    summary << "originals_compared\t" << compared << "\noriginals_mismatched\t" << mismatched
            << "\nworst_original_delta\t" << worst << "\noriginals_reproduce\t"
            << (mismatched == 0 ? "true" : "false") << '\n';
    if (mismatched != 0) {
        summary << "refused\ttrue\nrefusal\toriginal scores did not reproduce; repairs not interpreted\n";
        std::cout << "repair-pilot " << donor << ": REFUSED, " << mismatched
                  << " original scores did not reproduce (worst delta " << worst << ")\n";
        return 4;
    }

    // ---- nested prefixes, all certified from the one scoring pass. NO TRUTH IN SCOPE YET. ----------
    std::vector<std::size_t> prefixes;
    {
        std::stringstream spec(args.get("prefixes", "8,16,32"));
        for (std::string cell; std::getline(spec, cell, ',');) prefixes.push_back(std::stoul(cell));
    }
    struct PrefixResult {
        std::size_t prefix = 0;
        std::size_t candidates = 0;
        std::size_t repairs_included = 0;
        std::size_t winner_entry = 0;
        bool complete = false;
        double winner_score = 0.0;
        double separation = 0.0;
        std::size_t fragment_contenders = 0;
        std::string phase_support;
        std::vector<std::pair<double, double>> block_confidence;   // posterior, GQ
        std::vector<bool> block_read_supported;
        std::vector<bool> block_separated;
    };
    // Certify one prefix from the shared scores. `order` lets the candidate-order invariance gate
    // re-certify a permuted presentation without a single new likelihood call.
    const auto certify_prefix = [&](std::size_t prefix, bool reversed) {
        std::vector<std::size_t> chosen;
        for (std::size_t i = 0; i < entries.size(); ++i) {
            if (!entries[i].repair || entries[i].rank <= prefix) chosen.push_back(i);
        }
        if (reversed) std::reverse(chosen.begin(), chosen.end());
        std::map<std::size_t, std::size_t> slot_of_local;
        std::vector<panvar::DiplotypeScore> scores;
        std::vector<double> priors;
        std::vector<std::size_t> entry_of_local;
        for (const std::size_t i : chosen) {
            const std::size_t slot = pair_of_entry[i];
            if (slot_of_local.count(slot) != 0) continue;
            slot_of_local.emplace(slot, scores.size());
            scores.push_back(scored.scores[slot]);
            priors.push_back(prior_of_pair[slot]);
            entry_of_local.push_back(i);
        }
        const panvar::CertifiedRanking ranking = panvar::certify_scored_ranking(
            scores, fragments.size(), scoring.mismatch_band, scored.counters,
            panvar::CertificationOptions{}, priors);
        PrefixResult result;
        result.prefix = prefix;
        result.candidates = scores.size();
        for (const std::size_t i : chosen) result.repairs_included += entries[i].repair ? 1 : 0;
        result.winner_entry = entry_of_local[ranking.winner_set.front()];
        result.complete = ranking.status == panvar::RescoreStatus::Complete;
        result.winner_score = ranking.candidates.front().score;
        result.separation = ranking.certified_gap;

        // Ranked order, winner first: block_genotype_confidence compares every row against row 0.
        std::vector<double> ranked_scores;
        std::vector<double> ranked_bounds;
        std::vector<std::vector<std::pair<std::size_t, std::size_t>>> ranked_genotypes;
        const auto genotype_row = [&](std::size_t entry) {
            const std::vector<std::size_t>& a = paths[entries[entry].first_slot].alleles;
            const std::vector<std::size_t>& b2 = paths[entries[entry].second_slot].alleles;
            std::vector<std::pair<std::size_t, std::size_t>> row(blocks);
            for (std::size_t b = 0; b < blocks; ++b) row[b] = std::minmax(a[b], b2[b]);
            return row;
        };
        for (const panvar::CertifiedCandidate& candidate : ranking.candidates) {
            ranked_scores.push_back(candidate.score);
            ranked_bounds.push_back(candidate.numeric_error_bound);
            ranked_genotypes.push_back(genotype_row(entry_of_local[candidate.candidate_id]));
        }
        const std::vector<panvar::BlockGenotypeConfidence> confidence =
            panvar::block_genotype_confidence(ranked_scores, ranked_bounds, ranked_genotypes);
        result.block_read_supported.assign(blocks, true);
        for (const auto& item : confidence) {
            result.block_confidence.push_back({item.posterior, item.quality});
            result.block_separated.push_back(item.numerically_separated);
        }
        // A fragment contender is a rival the READS do not put below the winner, which is production's
        // condition verbatim; a block on which any contender disagrees is not read-supported.
        const panvar::CertifiedCandidate& best = ranking.candidates.front();
        for (std::size_t c = 1; c < ranking.candidates.size(); ++c) {
            const panvar::CertifiedCandidate& rival = ranking.candidates[c];
            if (rival.fragment_term + rival.numeric_error_bound <
                best.fragment_term - best.numeric_error_bound) {
                continue;
            }
            ++result.fragment_contenders;
            for (std::size_t b = 0; b < blocks; ++b) {
                if (ranked_genotypes[c][b] != ranked_genotypes[0][b]) result.block_read_supported[b] = false;
            }
        }
        result.phase_support = !result.complete ? "UNRESOLVED"
                               : result.fragment_contenders == 0 ? "READS_WITHIN_CANDIDATES"
                                                                 : "PRIOR_ONLY";
        return result;
    };

    std::vector<PrefixResult> results;
    for (const std::size_t prefix : prefixes) results.push_back(certify_prefix(prefix, false));

    // ---- invariance gates, both from the already-computed scores: no new likelihood call ------------
    // Candidate order: the same prefix presented back to front must certify the same winner and score.
    bool order_invariant = true;
    for (const PrefixResult& forward : results) {
        const PrefixResult back = certify_prefix(forward.prefix, true);
        const std::vector<std::size_t>& f1 = paths[entries[forward.winner_entry].first_slot].alleles;
        const std::vector<std::size_t>& f2 = paths[entries[forward.winner_entry].second_slot].alleles;
        const std::vector<std::size_t>& b1 = paths[entries[back.winner_entry].first_slot].alleles;
        const std::vector<std::size_t>& b2 = paths[entries[back.winner_entry].second_slot].alleles;
        if (back.candidates != forward.candidates || back.winner_score != forward.winner_score ||
            std::minmax(f1, f2) != std::minmax(b1, b2)) {
            order_invariant = false;
        }
    }
    // Homologue swap: rebuilding the construction with each seed's two haplotypes exchanged must yield
    // the same candidate key set, because both assignments are generated for every local pair.
    bool swap_invariant = true;
    {
        std::set<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> swapped;
        for (const std::size_t seed : seeds) {
            const Original& base = originals[seed];
            for (std::size_t b = 0; b < blocks; ++b) {
                const panvar::BlockProposalEmissions& emission = emissions.blocks[b];
                const std::size_t alleles = emission.allele_count;
                std::vector<std::pair<double, std::pair<std::size_t, std::size_t>>> ranked;
                for (std::size_t x = 0; x < alleles; ++x) {
                    for (std::size_t y = x; y < alleles; ++y) ranked.push_back({emission.at(x, y), {x, y}});
                }
                std::sort(ranked.begin(), ranked.end(), [](const auto& lhs, const auto& rhs) {
                    if (lhs.first != rhs.first) return lhs.first > rhs.first;
                    return lhs.second < rhs.second;
                });
                const std::size_t take = std::min(local_pairs, ranked.size());
                for (std::size_t r = 0; r < take; ++r) {
                    const std::size_t x = ranked[r].second.first;
                    const std::size_t y = ranked[r].second.second;
                    for (int assignment = 0; assignment < 2; ++assignment) {
                        if (assignment == 1 && x == y) continue;
                        std::vector<std::size_t> first = base.second;   // homologues exchanged
                        std::vector<std::size_t> second = base.first;
                        first[b] = assignment == 0 ? x : y;
                        second[b] = assignment == 0 ? y : x;
                        swapped.insert(key_of(first, second));
                    }
                }
            }
        }
        std::set<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> forward;
        for (std::size_t i = original_entries; i < entries.size(); ++i) {
            forward.insert(key_of(paths[entries[i].first_slot].alleles,
                                  paths[entries[i].second_slot].alleles));
        }
        for (const auto& key : forward) {
            if (swapped.count(key) == 0) swap_invariant = false;
        }
    }
    summary << "order_invariant\t" << (order_invariant ? "true" : "false") << "\nswap_invariant\t"
            << (swap_invariant ? "true" : "false") << '\n';

    // ---- the candidate table: TRUTH-FREE, and written before truth is loaded -----------------------
    {
        std::ofstream candidates(out + ".candidates.tsv");
        candidates << std::setprecision(12)
                   << "candidate_index\tsource\tseed_id\trepaired_block\tlocal_pair_rank\tassignment\t"
                      "allele_x\tallele_y\thap1_alleles\thap2_alleles\tfragment_term\texposure_term\t"
                      "prior\ttotal";
        for (const PrefixResult& r : results) candidates << "\tin_R" << r.prefix;
        for (const PrefixResult& r : results) candidates << "\twinner_R" << r.prefix;
        candidates << '\n';
        const auto join = [](const std::vector<std::size_t>& v) {
            std::string text;
            for (std::size_t i = 0; i < v.size(); ++i) {
                if (i != 0) text += ',';
                text += std::to_string(v[i]);
            }
            return text;
        };
        for (std::size_t i = 0; i < entries.size(); ++i) {
            const Entry& entry = entries[i];
            const std::size_t slot = pair_of_entry[i];
            const double total = scored.scores[slot].log_likelihood + prior_of_pair[slot];
            candidates << i << '\t' << (entry.repair ? "repair" : "original") << '\t'
                       // original_entry_id holds only the originals prefix, so bound the index
                       // explicitly rather than leaning on the ternary's short circuit.
                       << (entry.repair || i >= original_entry_id.size() ? entry.seed_id
                                                                        : original_entry_id[i])
                       << '\t'
                       << (entry.repair ? model.blocks()[entry.block].block_id : std::string("NA")) << '\t'
                       << (entry.repair ? std::to_string(entry.rank) : std::string("NA")) << '\t'
                       << (entry.repair ? std::to_string(entry.assignment) : std::string("NA")) << '\t'
                       << (entry.repair ? std::to_string(entry.allele_x) : std::string("NA")) << '\t'
                       << (entry.repair ? std::to_string(entry.allele_y) : std::string("NA")) << '\t'
                       << join(paths[entry.first_slot].alleles) << '\t'
                       << join(paths[entry.second_slot].alleles) << '\t'
                       << scored.scores[slot].fragment_term << '\t' << scored.scores[slot].exposure_term
                       << '\t' << prior_of_pair[slot] << '\t' << total;
            for (const PrefixResult& r : results) {
                candidates << '\t' << ((!entry.repair || entry.rank <= r.prefix) ? "true" : "false");
            }
            for (const PrefixResult& r : results) {
                candidates << '\t' << (r.winner_entry == i ? "true" : "false");
            }
            candidates << '\n';
        }
    }

    // ---- truth enters ONLY here, after every candidate is built, scored and written ----------------
    const SortedGraph truth_graph = load_truth_in_model_frame(
        args.required("truth-gfa"), reference, loaded.sorted, out + ".truth.sorted.gfa");
    std::vector<panvar::PathProjection> truth;
    for (const panvar::PathRecord* path : donor_paths(truth_graph.graph, donor)) {
        const auto projection = panvar::project_path_onto_locus(truth_graph.graph, loaded.bubbles,
                                                                loaded.built.report, path->name);
        if (!projection.has_value()) throw std::runtime_error("truth does not span the block chain");
        truth.push_back(*projection);
    }
    const auto spelling = [&](std::size_t b, std::size_t allele) -> const std::string& {
        return model.blocks()[b].alleles[allele].sequence;
    };
    const auto block_distance = [&](std::size_t b, std::size_t a1, std::size_t a2) -> long {
        if (!truth[0].blocks[b].has_value() || !truth[1].blocks[b].has_value()) return -1;
        const std::string& t1 = *truth[0].blocks[b];
        const std::string& t2 = *truth[1].blocks[b];
        return std::min(edit_distance(t1, spelling(b, a1)) + edit_distance(t2, spelling(b, a2)),
                        edit_distance(t1, spelling(b, a2)) + edit_distance(t2, spelling(b, a1)));
    };
    std::vector<long> block_floor(blocks, 0);
    long floor_total = 0;
    for (std::size_t b = 0; b < blocks; ++b) {
        if (!truth[0].blocks[b].has_value() || !truth[1].blocks[b].has_value()) continue;
        bool available = false;
        block_floor[b] = closest_allele_distance(*truth[0].blocks[b], model.blocks()[b], available) +
                         closest_allele_distance(*truth[1].blocks[b], model.blocks()[b], available);
        floor_total += block_floor[b];
    }
    summary << "availability_floor_total\t" << floor_total << '\n';

    std::ofstream table(out + ".prefixes.tsv");
    table << std::setprecision(12)
          << "local_pairs\tcandidates\trepairs_included\tstatus\twinner_is_repair\twinner_entry\t"
             "winner_score\tfragment_term\texposure_term\tprior\tseparation\tfragment_contenders\t"
             "phase_support\tdistance\tfloor\texcess\n";
    std::ofstream per_block(out + ".truth.tsv");
    per_block << std::setprecision(12)
              << "local_pairs\tblock\tblock_id\tkind\twinner_allele1\twinner_allele2\tblock_distance\t"
                 "block_floor\tblock_excess\tgenotype_support\tposterior\tGQ\tnumerically_separated\n";
    for (const PrefixResult& r : results) {
        const std::size_t slot = pair_of_entry[r.winner_entry];
        const std::vector<std::size_t>& w1 = paths[entries[r.winner_entry].first_slot].alleles;
        const std::vector<std::size_t>& w2 = paths[entries[r.winner_entry].second_slot].alleles;
        long distance = 0;
        bool defined = true;
        for (std::size_t b = 0; b < blocks; ++b) {
            const long d = block_distance(b, w1[b], w2[b]);
            if (d < 0) {
                defined = false;
                continue;
            }
            distance += d;
            per_block << r.prefix << '\t' << b << '\t' << model.blocks()[b].block_id << '\t'
                      << (model.blocks()[b].kind == panvar::GenotypeBlockKind::Variant ? "SV" : "context")
                      << '\t' << w1[b] << '\t' << w2[b] << '\t' << d << '\t' << block_floor[b] << '\t'
                      << (d - block_floor[b]) << '\t'
                      << (r.block_read_supported[b] ? "READS" : "PRIOR_ONLY") << '\t'
                      << r.block_confidence[b].first << '\t' << r.block_confidence[b].second << '\t'
                      << (r.block_separated[b] ? "true" : "false") << '\n';
        }
        table << r.prefix << '\t' << r.candidates << '\t' << r.repairs_included << '\t'
              << (r.complete ? "COMPLETE" : "RESCORE_INCOMPLETE") << '\t'
              << (entries[r.winner_entry].repair ? "true" : "false") << '\t' << r.winner_entry << '\t'
              << r.winner_score << '\t' << scored.scores[slot].fragment_term << '\t'
              << scored.scores[slot].exposure_term << '\t' << prior_of_pair[slot] << '\t' << r.separation
              << '\t' << r.fragment_contenders << '\t' << r.phase_support << '\t'
              << (defined ? std::to_string(distance) : std::string("NA")) << '\t' << floor_total << '\t'
              << (defined ? std::to_string(distance - floor_total) : std::string("NA")) << '\n';
    }
    std::cout << "repair-pilot " << donor << ": " << repair_entries << " repairs kept of " << generated
              << " generated, " << scored_pairs.size() << " scored pairs, "
              << interned.sequences.size() << " unique haplotypes, " << scoring_seconds
              << "s scoring, order_invariant=" << (order_invariant ? "true" : "false")
              << " swap_invariant=" << (swap_invariant ? "true" : "false") << '\n';
    return 0;
}

// ------------------------------------------------------------------------------------------------
// repair-search: greedy coordinate repair, one exact scoring pass per round
//
// EXPERIMENTAL EVALUATION PATH. Not wired into the production caller, and deliberately not shaped
// like it: production would hold a persistent scoring context across rounds, whereas every round here
// is a cold, self-contained score_diplotypes call over the current candidate plus its neighbours.
// That is slower, and it is the only way to be sure a round's numbers owe nothing to a cache carried
// in from the round before.
//
// The move rule is coordinate ascent on the complete score. One block changes at a time, using the
// top local allele pairs the block's own proposal emissions rank, in both homologue assignments. The
// search moves only to a UNIQUE best neighbour that strictly improves, so a tie is a stop rather than
// an arbitrary choice.
//
// A block is never barred from being repaired again. Changing one homologue at block b can change
// which assignment block b' prefers, and that can in turn change what b wants -- barring revisits
// would make the reachable set depend on the order blocks happen to be numbered. Cycles are prevented
// instead by a visited set of allele-path pairs, which suppresses a neighbour that has already been
// occupied rather than forbidding a coordinate.
//
// TRUTH IS NOT IN SCOPE during the search. The rounds table is complete and written before the truth
// graph is opened; the distance tables are produced afterwards by the evaluator alone. --truth-gfa is
// optional precisely so this is testable: running with and without it must produce the identical
// rounds table.
int repair_search(const Args& args) {
    const std::string out = args.required("out");
    struct StagedGraphs {
        std::string prefix;
        ~StagedGraphs() {
            std::error_code ignored;
            std::filesystem::remove(prefix + ".model.sorted.gfa", ignored);
            std::filesystem::remove(prefix + ".truth.sorted.gfa", ignored);
        }
    } staged_graphs{out};
    const std::string reference = args.required("reference");
    const std::string donor = args.required("donor");
    const std::string run = args.required("run");
    const std::string truth_gfa = args.get("truth-gfa", "");
    const std::size_t local_pairs = std::stoul(args.get("local-pairs", "32"));
    const std::size_t max_rounds = std::stoul(args.get("rounds", "4"));
    // --context 1 drives the search through one PreparedScoringContext instead of a cold
    // score_diplotypes per round. The candidate search, the move rule and the model are untouched:
    // the only difference is who owns the read forms, the piece indexes and the haplotype columns.
    const bool use_context = args.get("context", "0") != "0";
    // --delta 1 scores each neighbour as a single-block replacement of the incumbent, reusing every
    // origin the replaced block cannot touch. It implies --context, because a delta needs somewhere
    // for the incumbent's census and columns to live between rounds.
    const bool use_delta = args.get("delta", "0") != "0";
    // Stage 6, diagnostic only: the window-origin cache and the preaggregated merge. Requires --delta.
    const bool use_origin_cache = args.get("origin-cache", "0") != "0";
    const std::size_t origin_cache_bytes =
        static_cast<std::size_t>(std::stoull(args.get("origin-cache-bytes", "0")));
    const std::uint64_t origin_cache_records = std::stoull(args.get("origin-cache-records", "0"));
    // Retention policy, experimental: first_fit (stage 6, the default) keeps any miss that fits;
    // planned selects each round's cache set from the round's planned windows by reuse per exact byte.
    const std::string origin_cache_policy = args.get("origin-cache-policy", "first_fit");
    if (origin_cache_policy != "first_fit" && origin_cache_policy != "planned") {
        throw std::runtime_error("--origin-cache-policy is first_fit or planned");
    }
    const bool planned_cache = origin_cache_policy == "planned";
    // ADVERSARIAL GATE HOOK: permute the neighbour list and the repair-request list every round. Under
    // planned retention nothing observable may change -- cache identities, accounting, scores, ranking.
    const std::uint64_t request_order_seed = std::stoull(args.get("request-order-seed", "0"));
    // DIAGNOSTIC: bounded, source-diverse multi-start. Seeds come from the exactly rescored panel+mosaic
    // union of --run: the top equivalence set, the best panel candidate with its ties, the best mosaic
    // candidate with its ties. Each is repaired independently; every scored candidate is pooled and the
    // pool is exactly reranked. Nothing here reads truth.
    const bool multi_start = args.get("multi-start", "0") != "0";
    const std::size_t max_repair_starts = std::stoul(args.get("max-repair-starts", "16"));
    const bool reverse_seeds = args.get("reverse-seeds", "0") != "0";
    const bool isolate_seeds = args.get("isolate-seeds", "0") != "0";
    const bool seeds_only = args.get("seeds-only", "0") != "0";
    const std::uint64_t shuffle_union = std::stoull(args.get("shuffle-union-input", "0"));
    const std::string explicit_seeds = args.get("seed-candidates", "");
    // Per-seed work bounds (0: none), checked at every committed round; the process-wide bound is the
    // context's own cumulative caps, enforced in flight. Crossing either refuses the whole search.
    const std::uint64_t seed_max_origins = std::stoull(args.get("seed-max-physical-origins", "0"));
    const std::uint64_t seed_max_classes = std::stoull(args.get("seed-max-state-classes", "0"));
    const std::uint64_t seed_max_piece_hits = std::stoull(args.get("seed-max-piece-hits", "0"));
    const std::uint64_t seed_max_latent = std::stoull(args.get("seed-max-latent-pairs", "0"));

    // MEASUREMENT ONLY: evaluate the necessary-condition piece screen in shadow on every verification
    // (it rejects nothing), and export named recorded windows for the enumeration benchmark.
    panvar::reset_peak_active_workers();
    const bool piece_screen_shadow = args.get("piece-screen-shadow", "0") != "0";
    if (piece_screen_shadow) panvar::set_piece_screen_shadow(true);
    std::vector<std::pair<std::string, std::uint64_t>> export_windows;
    {
        std::stringstream spec(args.get("export-windows", ""));
        for (std::string cell; std::getline(spec, cell, ',');) {
            const std::size_t eq = cell.find('=');
            if (eq == std::string::npos) throw std::runtime_error("--export-windows wants label=keyhash,...");
            export_windows.emplace_back(cell.substr(0, eq), std::stoull(cell.substr(eq + 1)));
        }
    }
    if (args.get("origin-self-check", "0") != "0") panvar::set_origin_cache_self_check(true);
    // Memory diagnostic: after the search stops without moving, re-score that same neighbourhood N more
    // times, scoring and releasing without promoting. Hold cycles never enter the rounds table.
    const std::size_t hold_cycles = std::stoul(args.get("hold-cycles", "0"));
    const bool pressure_relief = args.get("pressure-relief", "0") != "0";
    const double switch_probability = std::stod(args.get("switch-probability", "0.01"));
    const std::uint64_t declared_cap = std::stoull(args.get("max-neighbors", "0"));
    // Resumption. --start declares where to continue from and --start-score what it must score; the
    // incumbent guard in round one then IS the reproduction requirement, refusing unless the declared
    // value comes back exactly. Resuming without also restoring --visited would let the search walk
    // back onto candidates the earlier leg already occupied, so the caller restores both.
    const std::string declared_start = args.get("start", "");
    const bool resuming = !declared_start.empty();
    const std::size_t first_round = std::stoul(args.get("first-round", "1"));

    // ---- run state. A run is RUNNING until it is COMPLETE, TIME_LIMIT, RESOURCE_REFUSAL or FAILED.
    // Every completed round is committed atomically to <out>.*.partial.tsv as it finishes, and the
    // status file says what those rows are. Only a normal completion writes the final names and the
    // summary -- so no partial trajectory can be read as a completed call, whoever stops the run.
    const auto run_started = std::chrono::steady_clock::now();
    const auto elapsed = [&run_started]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - run_started).count();
    };
    const double time_limit = std::stod(args.get("time-limit-seconds", "0"));
    const std::uint64_t footprint_limit = std::stoull(args.get("footprint-limit-bytes", "0"));
    // TEST HOOK: stop through the time-limit path after this many committed rounds, so the gate can
    // exercise an interruption deterministically. It is reported as such in the status reason.
    const std::size_t interrupt_after = std::stoul(args.get("interrupt-after-rounds", "0"));
    std::mutex commit_mutex;
    std::atomic<bool> finished{false};
    std::size_t rounds_committed = 0;
    std::string status_extra;
    const auto write_status = [&](const std::string& state_name, const std::string& reason) {
        const FootprintSample footprint = physical_footprint();
        const bool complete = state_name == "COMPLETE";
        std::ostringstream text;
        text << std::setprecision(12) << "key\tvalue\nstatus\t" << state_name << "\nreason\t" << reason
             << "\ncomplete\t" << (complete ? "true" : "false") << "\nrounds_committed\t"
             << rounds_committed << "\nelapsed_seconds\t" << elapsed() << "\ntime_limit_seconds\t"
             << time_limit << "\nfootprint_limit_bytes\t" << footprint_limit
             << "\nphys_footprint_bytes\t" << footprint.current << "\nphys_footprint_peak_bytes\t"
             << footprint.peak << "\nrss_current_bytes\t" << current_rss_bytes()
             << "\nrss_peak_bytes\t" << static_cast<std::uint64_t>(peak_rss_gib() * 1073741824.0)
             << "\noutputs\t"
             << (complete ? "final: <out>.rounds.tsv, .summary.tsv and companions"
                          : "PARTIAL: <out>.*.partial.tsv hold committed rounds only; no summary, not a call")
             << '\n' << status_extra;
        atomic_write(out + ".run_status.tsv", text.str());
    };
    // Caller must NOT hold the lock. After this, the watchdog can no longer overwrite the status.
    const auto finish = [&](const std::string& state_name, const std::string& reason) {
        std::lock_guard<std::mutex> lock(commit_mutex);
        finished = true;
        write_status(state_name, reason);
    };
    // Stops the process where it stands. Rounds already committed stay on disk, intact.
    const auto stop_now = [&](const std::string& state_name, const std::string& reason, int code) {
        std::lock_guard<std::mutex> lock(commit_mutex);
        if (finished) return;
        finished = true;
        write_status(state_name, reason);
        std::cout << "repair-search: " << state_name << ": " << reason << " (" << rounds_committed
                  << " rounds committed)" << std::endl;
        // _Exit runs no destructor, so the staged graphs are removed here. The partial tables and the
        // status are the diagnostic record and stay.
        std::error_code ignored;
        std::filesystem::remove(out + ".model.sorted.gfa", ignored);
        std::filesystem::remove(out + ".truth.sorted.gfa", ignored);
        std::_Exit(code);
    };
    {
        std::lock_guard<std::mutex> lock(commit_mutex);
        write_status("RUNNING", "started");
    }
    std::thread watchdog;
    if (time_limit > 0.0 || footprint_limit != 0) {
        watchdog = std::thread([&]() {
            while (!finished.load()) {
                if (time_limit > 0.0 && elapsed() > time_limit) {
                    std::ostringstream why;
                    why << "wall " << elapsed() << " s exceeded the declared " << time_limit << " s";
                    stop_now("TIME_LIMIT", why.str(), 124);
                }
                if (footprint_limit != 0) {
                    const std::uint64_t now = physical_footprint().current;
                    if (now > footprint_limit) {
                        std::ostringstream why;
                        why << "physical footprint " << now << " bytes exceeded the declared "
                            << footprint_limit;
                        stop_now("RESOURCE_REFUSAL", why.str(), 125);
                    }
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
        });
    }
    struct WatchdogJoin {
        std::atomic<bool>& done;
        std::thread& thread;
        ~WatchdogJoin() {
            done = true;
            if (thread.joinable()) thread.join();
        }
    } watchdog_join{finished, watchdog};
    // Any exit that did not set a final status -- an exception thrown outside the search loop, or an
    // early return -- leaves FAILED rather than a stale RUNNING. Unwinding then runs StagedGraphs.
    struct UnfinishedGuard {
        const std::function<void()> on_exit;
        ~UnfinishedGuard() { on_exit(); }
    } unfinished_guard{[&]() {
        std::lock_guard<std::mutex> lock(commit_mutex);
        if (finished) return;
        finished = true;
        write_status("FAILED", std::uncaught_exceptions() > 0
                                   ? "an exception ended the run before a final status"
                                   : "the run returned without a final status");
    }};

    const std::string curation_path = args.get("panel-curation", "");
    const LoadedModel loaded =
        load_model(args.required("model-gfa"), reference, out + ".model.sorted.gfa", 50,
                   curation_path.empty() ? std::vector<std::string>{} : read_panel_curation(curation_path));
    const panvar::PreparedGenotypeModel& model = *loaded.model;
    const std::size_t blocks = model.blocks().size();
    const auto status = read_key_values(run + ".mosaic_status.tsv");
    if (status.at("model_fingerprint") != loaded.fingerprint) {
        throw std::runtime_error("the run was not scored on the rebuilt model");
    }

    const panvar::LibraryModel library = read_library_file(args.required("library"));
    const std::vector<panvar::ReadPair> fragments =
        read_pairs(args.required("reads1"), args.required("reads2"));

    // Block proposal emissions, rebuilt as the run built them; the fit is checked against the run's
    // recorded rate and background, because different emissions would rank different local pairs.
    const panvar::PreparedMarkerPanel markers = panvar::build_marker_panel(model, 31, 15);
    const panvar::MarkerObservations observations = panvar::count_marker_observations(markers, fragments);
    std::optional<double> declared_background;
    {
        std::ifstream file(args.required("library"));
        std::string key;
        double value = 0.0;
        for (std::string line; std::getline(file, line);) {
            std::istringstream fields(line);
            if (fields >> key >> value && key == "marker_background") declared_background = value;
        }
    }
    const panvar::MarkerDepthEstimate depth =
        panvar::estimate_marker_depth(markers, observations, declared_background);
    if (!depth.converged) throw std::runtime_error("the marker depth fit did not converge");
    const panvar::PreparedProposalEmissions emissions = panvar::build_proposal_emissions(
        markers, observations, panvar::proposal_marker_parameters(depth, 31, 15));

    // ---- the start: the best EXACTLY SCORED panel candidate, read off the signature-bound campaign
    struct Original {
        std::vector<std::size_t> first;
        std::vector<std::size_t> second;
        std::string source;
    };
    std::map<std::string, Original> original_of_id;
    std::vector<std::string> union_ids;   // the union's candidates in the call's own table order
    for (const auto& row : read_table(run + ".mosaic_candidates.tsv")) {
        Original entry;
        entry.first = parse_alleles(row.at("hap1_alleles"));
        entry.second = parse_alleles(row.at("hap2_alleles"));
        entry.source = row.count("source") != 0 ? row.at("source") : std::string("unknown");
        union_ids.push_back(row.at("candidate_id"));
        original_of_id.emplace(row.at("candidate_id"), std::move(entry));
    }
    std::string start_id;
    double start_recorded = 0.0;
    bool have_start = false;
    for (const auto& row : read_table(run + ".mosaic_scores.tsv")) {
        const auto found = original_of_id.find(row.at("representative_beam_candidate"));
        if (found == original_of_id.end() || found->second.source != "panel") continue;
        const double score = std::stod(row.at("score"));
        if (!have_start || score > start_recorded) {
            have_start = true;
            start_recorded = score;
            start_id = row.at("representative_beam_candidate");
        }
    }
    if (!have_start && !resuming) {
        throw std::runtime_error("the run has no exactly scored panel candidate");
    }

    std::vector<std::size_t> start_first;
    std::vector<std::size_t> start_second;
    double start_total = 0.0;
    if (resuming) {
        const std::size_t bar = declared_start.find('|');
        if (bar == std::string::npos) throw std::runtime_error("--start wants <hap1>|<hap2>");
        start_first = parse_alleles(declared_start.substr(0, bar));
        start_second = parse_alleles(declared_start.substr(bar + 1));
        if (start_first.size() != blocks || start_second.size() != blocks) {
            throw std::runtime_error("--start needs one allele per block");
        }
        start_total = std::stod(args.required("start-score"));
    } else {
        start_first = original_of_id.at(start_id).first;
        start_second = original_of_id.at(start_id).second;
        start_total = start_recorded;
    }

    // ---- the declared bound on ONE round's neighbourhood, from the block shape rather than assumed
    const std::uint64_t cap = panvar::repair_neighbour_bound(model, local_pairs);
    // Assembled in memory and written only by a normal completion: a partial summary is exactly the
    // kind of file that could be mistaken for a finished run.
    std::ostringstream summary;
    summary << std::setprecision(12) << "key\tvalue\n";
    summary << "donor\t" << donor << "\nrun\t" << run << "\nlocal_pairs\t" << local_pairs
            << "\nrounds_requested\t" << max_rounds << '\n';
    summary << "marker_depth_rate_per_copy\t" << depth.rate_per_copy << "\nmarker_background\t"
            << depth.background << "\nrecorded_marker_depth_rate_per_copy\t"
            << status.at("marker_depth_rate_per_copy") << "\nrecorded_marker_background\t"
            << status.at("marker_background") << '\n';
    summary << "resumed\t" << (resuming ? "true" : "false") << "\nfirst_round\t" << first_round
            << "\nstart_candidate\t" << (resuming ? std::string("declared") : start_id)
            << "\nstart_required_score\t" << start_total << '\n';
    {
        std::string a, b;
        for (std::size_t i = 0; i < blocks; ++i) {
            if (i != 0) { a += ','; b += ','; }
            a += std::to_string(start_first[i]);
            b += std::to_string(start_second[i]);
        }
        summary << "start_hap1\t" << a << "\nstart_hap2\t" << b << '\n';
    }
    summary << "neighbor_cap\t" << cap << "\ndeclared_cap\t" << declared_cap << '\n';
    summary << "origin_cache_policy\t" << origin_cache_policy << "\norigin_cache_budget_bytes\t"
            << origin_cache_bytes << '\n';
    if (declared_cap != 0 && cap > declared_cap) {
        summary << "refused\ttrue\nrefusal\tone round's neighbourhood exceeds its declared bound\n";
        std::cout << "repair-search " << donor << ": REFUSED before allocation, " << cap
                  << " neighbours per round exceeds the declared " << declared_cap << '\n';
        finish("RESOURCE_REFUSAL", "one round's neighbourhood (" + std::to_string(cap) +
                                       ") exceeds its declared bound " + std::to_string(declared_cap));
        return 3;
    }

    const auto key_of = [](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
        return a <= b ? std::make_pair(a, b) : std::make_pair(b, a);
    };
    std::vector<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> preoccupied_keys;
    // --visited pre-occupies candidates, so a gate can prove the set really suppresses a neighbour
    // that would otherwise have been selected. Nothing in the search itself writes to it this way.
    std::size_t preoccupied = 0;
    {
        // Semicolon-separated, because an allele path is itself comma-separated.
        std::stringstream spec(args.get("visited", ""));
        for (std::string cell; std::getline(spec, cell, ';');) {
            if (cell.empty()) continue;
            const std::size_t bar = cell.find('|');
            if (bar == std::string::npos) {
                throw std::runtime_error("--visited wants <hap1>|<hap2>, entries separated by ';'");
            }
            preoccupied_keys.push_back(key_of(parse_alleles(cell.substr(0, bar)), parse_alleles(cell.substr(bar + 1))));
            ++preoccupied;
        }
    }
    summary << "preoccupied_visited\t" << preoccupied << '\n';

    const std::vector<double> switches(blocks - 1, switch_probability);
    panvar::FragmentScoringOptions scoring;
    scoring.threads = std::stoul(args.get("threads", "8"));
    scoring.mismatch_band = std::stoul(status.at("mismatch_band"));
    scoring.minimum_piece_bp = 12;
    scoring.max_piece_hits = std::stoull(status.at("max_piece_hits"));
    scoring.max_proposed_starts_per_haplotype = std::stoull(status.at("max_proposed_starts_per_haplotype"));
    scoring.max_physical_origins = std::stoull(status.at("max_physical_origins"));
    scoring.max_state_classes = std::stoull(status.at("max_state_classes"));
    scoring.max_mass_cells = std::stoull(status.at("max_mass_cells"));
    // A run scored before the emission cache existed records no budget. The cache never changes a score,
    // but its budget is still declared, never guessed: such a run needs --emission-cache-bytes.
    if (args.get("emission-cache-bytes", "").empty() && status.count("emission_cache_bytes") == 0) {
        throw std::runtime_error("the run records no emission_cache_bytes; declare --emission-cache-bytes");
    }
    scoring.emission_cache_bytes = static_cast<std::size_t>(std::stoull(
        args.get("emission-cache-bytes", status.count("emission_cache_bytes") != 0 ? status.at("emission_cache_bytes") : "0")));
    // A CUMULATIVE budget is a different declared quantity from a per-call one. The campaign's caps
    // bound one scoring pass; a seven-round search through a persistent context sums seven of them
    // against a single budget, and inheriting the per-call figure silently would refuse a search that
    // no round of it exceeded. These overrides are how that difference gets declared out loud.
    // Named to avoid the neighbourhood bound already called `cap` in this function: that one bounds
    // how many candidates a round may construct, this one bounds the work the whole search may spend.
    const auto declared_budget = [&](const char* name, std::uint64_t current) {
        return std::stoull(args.get(name, std::to_string(current)));
    };
    scoring.max_physical_origins =
        declared_budget("max-physical-origins", scoring.max_physical_origins);
    scoring.max_state_classes = declared_budget("max-state-classes", scoring.max_state_classes);
    scoring.max_piece_hits = declared_budget("max-piece-hits", scoring.max_piece_hits);
    scoring.max_latent_alignment_pairs =
        declared_budget("max-latent-pairs", scoring.max_latent_alignment_pairs);
    scoring.max_mass_cells = declared_budget("max-mass-cells", scoring.max_mass_cells);

    // The search itself is the library's, shared with panvar genotype. Only the driver lives here.
    using Round = panvar::RepairRound;
    panvar::RepairSearchOptions search_options;
    search_options.local_pairs = local_pairs;
    search_options.max_rounds = max_rounds;
    search_options.first_round = first_round;
    search_options.mode = use_delta     ? panvar::RepairScoringMode::Delta
                          : use_context ? panvar::RepairScoringMode::Context
                                        : panvar::RepairScoringMode::OneShot;
    search_options.origin_cache = use_origin_cache;
    search_options.planned_cache = planned_cache;
    search_options.origin_cache_bytes = origin_cache_bytes;
    search_options.origin_cache_records = origin_cache_records;
    search_options.switch_probability = switch_probability;
    search_options.declared_neighbour_cap = declared_cap;
    search_options.request_order_seed = request_order_seed;
    search_options.hold_cycles = hold_cycles;
    search_options.pressure_relief = pressure_relief;
    panvar::CoordinateRepairSearch search(model, emissions, fragments, library, scoring, search_options);
    const std::size_t last_round = search.last_round();

    const auto join = [](const std::vector<std::size_t>& v) {
        std::string text;
        for (std::size_t i = 0; i < v.size(); ++i) {
            if (i != 0) text += ',';
            text += std::to_string(v[i]);
        }
        return text;
    };
    std::string rounds_prefix = out;          // where the current search's rounds are committed
    std::vector<Round> earlier_seed_rounds;   // multi-start: rounds of seeds already finished
    const auto margin = [](std::uint64_t cap_value, std::uint64_t used) -> std::string {
        if (cap_value == 0) return "none";
        return std::to_string(static_cast<std::int64_t>(cap_value) - static_cast<std::int64_t>(used));
    };

    // ---- the rounds table, rendered from the rounds that exist. TRUTH-FREE.
    const auto render_rounds = [&]() { return search.render_rounds(search.rounds()); };

    // ---- exact origin-cache demand, and its reconciliation against the work the rounds report
    const auto render_demand = [&](const std::vector<panvar::OriginDemandRow>& rows) {
        std::ostringstream table;
        table << std::setprecision(12)
              << "order\tkey_hash\twindow_bytes\tat_start\tat_end\trel_block_begin\trel_new_end\t"
                 "requests\trepeat_requests\thits\tmisses\tadmitted\tdeclined\tentry_bytes\trecords\t"
                 "origins_constructed\tseconds_enumeration\tseconds_reenumeration\t"
                 "aggregate_states\taggregate_fragments\taggregate_multiplicity\tretained_not_affected\t"
                 "aggregate_entry_bytes\taggregate_entry_bytes_packed\tgroup_reuses\tevictions\n";
        for (const panvar::OriginDemandRow& row : rows) {
            table << row.order << '\t' << row.key_hash << '\t' << row.window_bytes << '\t'
                  << (row.at_start ? 1 : 0) << '\t' << (row.at_end ? 1 : 0) << '\t' << row.rel_block_begin
                  << '\t' << row.rel_new_end << '\t' << row.requests << '\t' << (row.requests - 1) << '\t'
                  << row.hits << '\t' << row.misses << '\t' << (row.admitted ? "true" : "false") << '\t'
                  << row.declined << '\t' << row.entry_bytes << '\t' << row.records << '\t'
                  << row.origins_constructed << '\t' << row.seconds_enumeration << '\t'
                  << row.seconds_reenumeration << '\t' << row.aggregate_states << '\t'
                  << row.aggregate_fragments << '\t' << row.aggregate_multiplicity << '\t'
                  << row.retained_not_affected << '\t' << row.aggregate_entry_bytes << '\t'
                  << row.aggregate_entry_bytes_packed << '\t' << row.group_reuses << '\t' << row.evictions
                  << '\n';
        }
        return table.str();
    };
    const auto demand_report = [&](const std::vector<panvar::OriginDemandRow>& rows) {
        std::uint64_t requested_bytes = 0, requests = 0, repeats = 0, hits = 0, misses = 0;
        std::uint64_t admitted = 0, admitted_bytes = 0, declined = 0, declined_bytes = 0;
        std::uint64_t decline_events = 0, origins = 0, declined_repeats = 0, group_reuses = 0;
        double seconds = 0.0, reenumeration = 0.0, reenumeration_declined = 0.0;
        for (const panvar::OriginDemandRow& row : rows) {
            requested_bytes += row.entry_bytes;
            requests += row.requests;
            repeats += row.requests - 1;
            hits += row.hits;
            misses += row.misses;
            group_reuses += row.group_reuses;
            decline_events += row.declined;
            origins += row.origins_constructed;
            seconds += row.seconds_enumeration;
            reenumeration += row.seconds_reenumeration;
            if (row.admitted) {
                ++admitted;
                admitted_bytes += row.entry_bytes;
            } else if (row.declined != 0) {
                ++declined;
                declined_bytes += row.entry_bytes;
                declined_repeats += row.requests - 1;
                reenumeration_declined += row.seconds_reenumeration;
            }
        }
        std::ostringstream text;
        text << std::setprecision(12);
        if (rows.empty()) return std::string();
        text << "demand_distinct_windows\t" << rows.size() << "\ndemand_distinct_window_bytes\t"
             << requested_bytes << "\ndemand_requests\t" << requests << "\ndemand_repeat_requests\t"
             << repeats << "\ndemand_admitted_entries\t" << admitted << "\ndemand_admitted_bytes\t"
             << admitted_bytes << "\ndemand_declined_entries\t" << declined
             << "\ndemand_declined_bytes\t" << declined_bytes << "\ndemand_declined_repeat_requests\t"
             << declined_repeats << "\ndemand_decline_events\t" << decline_events
             << "\ndemand_seconds_enumeration\t" << seconds << "\ndemand_seconds_reenumeration\t"
             << reenumeration << "\ndemand_seconds_reenumeration_declined\t" << reenumeration_declined
             << '\n';
        // Reconciliation: the ledger is kept by the cache path, the round columns by the caller, and
        // the cumulative counters by the budget. Hold cycles add ledger traffic no round reports.
        if (search.hold_index() != 0) {
            text << "reconciled\tnot_applicable_hold_cycles\n";
            return text.str();
        }
        if (isolate_seeds) {
            text << "reconciled\tnot_applicable_isolated_contexts\n";
            return text.str();
        }
        std::uint64_t round_hits = 0, round_misses = 0, round_constructed = 0, round_declined = 0;
        std::uint64_t round_admitted_bytes = 0, round_evicted_bytes = 0, round_group_reuses = 0;
        std::vector<const Round*> every_round;
        for (const Round& r : earlier_seed_rounds) every_round.push_back(&r);
        for (const Round& r : search.rounds()) every_round.push_back(&r);
        for (const Round* each : every_round) {
            const Round& r = *each;
            round_hits += r.oc_hits_all;
            round_misses += r.oc_misses_all;
            round_constructed += r.origins_constructed_all;
            round_declined += r.oc_declined;
            round_admitted_bytes += r.oc_admitted_bytes;
            round_evicted_bytes += r.plan_evicted_bytes + r.oc_evicted_bytes;
            round_group_reuses += r.oc_group_reuses;
        }
        // Planned retention evicts, so what is cached is what was admitted less what was evicted; and
        // a ledger row marked admitted may since have been evicted, so the ledger's admitted total
        // bounds the cache only when nothing was.
        const std::uint64_t cached = search.context()->origin_cache_bytes();
        const std::uint64_t executed = search.context()->cumulative_counters().physical_origins - search.seeding_origins();
        const bool ok = round_hits == hits && round_misses == misses && round_constructed == origins &&
                        executed == origins && round_declined == decline_events &&
                        round_admitted_bytes - round_evicted_bytes == cached &&
                        (round_evicted_bytes != 0 || admitted_bytes == cached) &&
                        round_group_reuses == group_reuses && hits + misses + group_reuses == requests;
        text << "reconcile_group_reuses\t" << round_group_reuses << " rounds / " << group_reuses << " ledger\n";
        text << "reconcile_hits\t" << round_hits << " rounds / " << hits << " ledger\nreconcile_misses\t"
             << round_misses << " rounds / " << misses << " ledger\nreconcile_origins_constructed\t"
             << round_constructed << " rounds / " << origins << " ledger / " << executed
             << " cumulative-minus-seeding\nreconcile_declines\t" << round_declined << " rounds / "
             << decline_events << " ledger\nreconcile_admitted_bytes\t" << round_admitted_bytes
             << " admitted - " << round_evicted_bytes << " evicted (rounds) / " << admitted_bytes
             << " ledger-ever-admitted / " << cached << " cache\nreconciled\t" << (ok ? "true" : "false")
             << '\n';
        return text.str();
    };

    // ---- truth: opened only once a round is on disk, and only ever read to measure distance. Nothing
    // below is visible to the construction or the move rule above.
    std::optional<SortedGraph> truth_graph;
    std::vector<panvar::PathProjection> truth;
    std::vector<long> block_floor(blocks, 0);
    long floor_total = 0;
    const auto ensure_truth = [&]() {
        if (truth_graph.has_value()) return;
        truth_graph.emplace(
            load_truth_in_model_frame(truth_gfa, reference, loaded.sorted, out + ".truth.sorted.gfa"));
        for (const panvar::PathRecord* path : donor_paths(truth_graph->graph, donor)) {
            const auto projection = panvar::project_path_onto_locus(truth_graph->graph, loaded.bubbles,
                                                                    loaded.built.report, path->name);
            if (!projection.has_value()) throw std::runtime_error("truth does not span the block chain");
            truth.push_back(*projection);
        }
        for (std::size_t b = 0; b < blocks; ++b) {
            if (!truth[0].blocks[b].has_value() || !truth[1].blocks[b].has_value()) continue;
            bool available = false;
            block_floor[b] = closest_allele_distance(*truth[0].blocks[b], model.blocks()[b], available) +
                             closest_allele_distance(*truth[1].blocks[b], model.blocks()[b], available);
            floor_total += block_floor[b];
        }
    };
    struct TruthTables {
        std::string distance;
        std::string per_block;
        long start_distance = 0;
        long final_distance = 0;
    };
    const auto render_truth = [&]() {
        const auto spelling = [&](std::size_t b, std::size_t allele) -> const std::string& {
            return model.blocks()[b].alleles[allele].sequence;
        };
        const auto block_distance = [&](std::size_t b, std::size_t a1, std::size_t a2) -> long {
            if (!truth[0].blocks[b].has_value() || !truth[1].blocks[b].has_value()) return -1;
            const std::string& t1 = *truth[0].blocks[b];
            const std::string& t2 = *truth[1].blocks[b];
            return std::min(edit_distance(t1, spelling(b, a1)) + edit_distance(t2, spelling(b, a2)),
                            edit_distance(t1, spelling(b, a2)) + edit_distance(t2, spelling(b, a1)));
        };
        std::ostringstream distance;
        std::ostringstream per_block;
        distance << "round\thap1_alleles\thap2_alleles\tdistance\tfloor\texcess\n";
        per_block << "round\tblock\tblock_id\tkind\tallele1\tallele2\tblock_distance\tblock_floor\tblock_excess\n";
        const auto evaluate = [&](std::size_t index, const std::vector<std::size_t>& a,
                                  const std::vector<std::size_t>& b2) {
            long total = 0;
            for (std::size_t b = 0; b < blocks; ++b) {
                const long d = block_distance(b, a[b], b2[b]);
                if (d < 0) continue;
                total += d;
                per_block << index << '\t' << b << '\t' << model.blocks()[b].block_id << '\t'
                          << (model.blocks()[b].kind == panvar::GenotypeBlockKind::Variant ? "SV" : "context")
                          << '\t' << a[b] << '\t' << b2[b] << '\t' << d << '\t' << block_floor[b] << '\t'
                          << (d - block_floor[b]) << '\n';
            }
            distance << index << '\t' << join(a) << '\t' << join(b2) << '\t' << total << '\t' << floor_total
                     << '\t' << (total - floor_total) << '\n';
            return total;
        };
        TruthTables tables;
        tables.start_distance = evaluate(first_round - 1, start_first, start_second);
        tables.final_distance = tables.start_distance;
        for (const Round& r : search.rounds()) tables.final_distance = evaluate(r.index, r.first, r.second);
        tables.distance = distance.str();
        tables.per_block = per_block.str();
        return tables;
    };

    // ---- commit: every completed round, atomically, before the next one starts
    const auto commit = [&]() {
        std::vector<panvar::OriginDemandRow> demand_rows;
        if (search.context() && use_origin_cache) demand_rows = search.context()->origin_demand();
        const std::string report = demand_report(demand_rows);
        {
            std::lock_guard<std::mutex> lock(commit_mutex);
            if (finished) return;
            rounds_committed = search.rounds().size();
            atomic_write(rounds_prefix + ".rounds.partial.tsv", render_rounds());
            if (search.memory_on()) atomic_write(out + ".memory.partial.tsv", search.memory_text());
            if (!demand_rows.empty()) {
                atomic_write(out + ".origin_demand.partial.tsv", render_demand(demand_rows));
            }
            status_extra = report;
            write_status("RUNNING", "round " + std::to_string(search.rounds().back().index) + " committed");
        }
        // Multi-start evaluates truth once, after every seed has finished searching.
        if (truth_gfa.empty() || multi_start) return;
        ensure_truth();   // outside the lock: the round is already safe on disk
        const TruthTables tables = render_truth();
        std::lock_guard<std::mutex> lock(commit_mutex);
        if (finished) return;
        atomic_write(out + ".distance.partial.tsv", tables.distance);
        atomic_write(out + ".round_truth.partial.tsv", tables.per_block);
        status_extra = report + "last_committed_distance\t" + std::to_string(tables.final_distance) +
                       "\navailability_floor_total\t" + std::to_string(floor_total) +
                       "\nlast_committed_excess\t" + std::to_string(tables.final_distance - floor_total) +
                       '\n';
        write_status("RUNNING", "round " + std::to_string(search.rounds().back().index) + " committed, distance measured");
    };

    search.start(start_first, start_second, start_total);
    for (const auto& key : preoccupied_keys) search.add_visited(key.first, key.second);

    // ---- DIAGNOSTIC multi-start: the pool of every exactly scored candidate, keyed by complete identity
    using Identity = std::pair<std::vector<std::size_t>, std::vector<std::size_t>>;
    struct SeqId {
        std::uint64_t a = 0, b = 0;
        std::size_t length = 0;
        bool operator<(const SeqId& o) const { return std::tie(a, b, length) < std::tie(o.a, o.b, o.length); }
        bool operator==(const SeqId& o) const { return a == o.a && b == o.b && length == o.length; }
    };
    // Spelled sequences are grouped by a 128-bit content digest plus length; the allele paths are the
    // identity itself, compared exactly.
    const auto seq_id = [](const std::string& text) {
        SeqId id;
        id.length = text.size();
        std::uint64_t h1 = 1469598103934665603ULL;
        std::uint64_t h2 = 0x9e3779b97f4a7c15ULL;
        for (const char c : text) {
            const std::uint64_t v = static_cast<unsigned char>(c);
            h1 = (h1 ^ v) * 1099511628211ULL;
            h2 = (h2 + v + 1) * 0xbf58476d1ce4e5b9ULL;
            h2 ^= h2 >> 31;
        }
        id.a = h1;
        id.b = h2;
        return id;
    };
    struct PoolEntry {
        Identity key;
        SeqId first_seq, second_seq;
        panvar::DiplotypeScore score;
        double prior = 0.0;
        std::set<std::string> tags;
        bool union_member = false;
        std::string union_id;
        std::string union_source;
        std::size_t union_order = 0;
    };
    std::map<Identity, PoolEntry> pool;
    std::uint64_t pool_mismatches = 0;
    std::size_t current_seed = 0;
    panvar::FragmentScoringCounters seed_start;
    const auto same_score = [](const panvar::DiplotypeScore& x, const panvar::DiplotypeScore& y) {
        return std::memcmp(&x.log_likelihood, &y.log_likelihood, sizeof(double)) == 0 &&
               std::memcmp(&x.fragment_term, &y.fragment_term, sizeof(double)) == 0 &&
               std::memcmp(&x.exposure_term, &y.exposure_term, sizeof(double)) == 0 &&
               std::memcmp(&x.absolute_term_sum, &y.absolute_term_sum, sizeof(double)) == 0 &&
               x.fragments_with_mass == y.fragments_with_mass;
    };
    const auto pool_add = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b,
                              const std::string& sa, const std::string& sb, const panvar::DiplotypeScore& score,
                              double prior, const std::string& tag) -> PoolEntry& {
        const Identity key = key_of(a, b);
        auto found = pool.find(key);
        if (found == pool.end()) {
            PoolEntry entry;
            entry.key = key;
            const bool swapped = key.first != a;
            entry.first_seq = seq_id(swapped ? sb : sa);
            entry.second_seq = seq_id(swapped ? sa : sb);
            entry.score = score;
            entry.prior = prior;
            found = pool.emplace(key, std::move(entry)).first;
        } else if (!same_score(found->second.score, score) ||
                   std::memcmp(&found->second.prior, &prior, sizeof(double)) != 0) {
            ++pool_mismatches;   // one candidate scored twice must score identically
        }
        found->second.tags.insert(tag);
        return found->second;
    };
    const auto check_seed_caps = [&]() {
        if (!multi_start) return;
        const panvar::FragmentScoringCounters& c = search.context()->cumulative_counters();
        const auto over = [&](std::uint64_t now, std::uint64_t start, std::uint64_t cap, const char* what) {
            if (cap != 0 && now - start > cap) {
                throw panvar::ResourceExhausted(std::string("the per-seed ") + what + " bound was crossed by seed " +
                                                std::to_string(current_seed) +
                                                "; the whole adaptive search is refused");
            }
        };
        over(c.physical_origins, seed_start.physical_origins, seed_max_origins, "physical-origin");
        over(c.statistical_state_classes, seed_start.statistical_state_classes, seed_max_classes, "state-class");
        over(c.piece_hits, seed_start.piece_hits, seed_max_piece_hits, "piece-hit");
        over(c.latent_alignment_pairs, seed_start.latent_alignment_pairs, seed_max_latent, "latent-pair");
    };
    struct Seed {
        std::size_t index = 0;             // canonical position
        Identity key;
        std::string candidate_id;
        std::size_t union_order = 0;
        std::string union_source;
        std::size_t group_rank = 0;        // position of its sequence group in the exact union ranking
        std::string group_sources;         // panel, mosaic or both, over the group's members
        std::set<std::string> selected_as; // top_set, best_panel, best_mosaic
        panvar::CertifiedCandidate certified;
    };
    struct SeedOutcome {
        std::size_t seed = 0;
        std::vector<Round> rounds;
        std::string stop;
        Identity optimum;
        double optimum_total = 0.0;
        std::vector<Identity> incumbents;
        panvar::FragmentScoringCounters used;
        double seconds = 0.0;
    };
    struct PooledRanking {
        std::vector<std::vector<const PoolEntry*>> groups;
        std::vector<std::size_t> representative;
        panvar::CertifiedRanking ranking;
    };
    std::map<std::vector<std::size_t>, double> path_prior_cache;
    const auto path_prior = [&](const std::vector<std::size_t>& path) {
        const auto found = path_prior_cache.find(path);
        if (found != path_prior_cache.end()) return found->second;
        const double value = panvar::log_haploid_mosaic_prior(model.panel(), panvar::MosaicPath{path}, switches);
        path_prior_cache.emplace(path, value);
        return value;
    };
    // Sequence-level ranking, as the caller ranks: members spelling the same unordered sequence pair are
    // one state; the representative is the member with the highest exact path prior, ties to identity.
    const auto rank_pool = [&](bool union_only) {
        std::map<std::pair<SeqId, SeqId>, std::vector<const PoolEntry*>> by_sequence;
        for (const auto& item : pool) {
            const PoolEntry& e = item.second;
            if (union_only && !e.union_member) continue;
            SeqId x = e.first_seq;
            SeqId y = e.second_seq;
            if (y < x) std::swap(x, y);
            by_sequence[{x, y}].push_back(&e);
        }
        PooledRanking out;
        std::vector<panvar::DiplotypeScore> scores;
        std::vector<double> priors;
        for (auto& item : by_sequence) {
            const std::vector<const PoolEntry*>& members = item.second;
            std::size_t best = 0;
            double best_prior = -std::numeric_limits<double>::infinity();
            for (std::size_t m = 0; m < members.size(); ++m) {
                const Identity& k = members[m]->key;
                const double p = path_prior(k.first) + path_prior(k.second) + (k.first != k.second ? std::log(2.0) : 0.0);
                if (m == 0 || p > best_prior) {
                    best = m;
                    best_prior = p;
                }
                if (!same_score(members[m]->score, members.front()->score)) ++pool_mismatches;
            }
            out.groups.push_back(members);
            out.representative.push_back(best);
            scores.push_back(members.front()->score);
            priors.push_back(members.front()->prior);
        }
        out.ranking = panvar::certify_scored_ranking(scores, fragments.size(), scoring.mismatch_band,
                                                     panvar::FragmentScoringCounters{},
                                                     panvar::CertificationOptions{}, priors);
        return out;
    };
    const auto genotype_of = [&](const Identity& k) {
        std::string text;
        for (std::size_t b = 0; b < blocks; ++b) {
            if (b != 0) text += ',';
            const auto pair = std::minmax(k.first[b], k.second[b]);
            text += std::to_string(pair.first) + "/" + std::to_string(pair.second);
        }
        return text;
    };
    const auto classify = [&](const PooledRanking& r) -> std::string {
        if (r.ranking.status == panvar::RescoreStatus::Complete) return "SEPARATED_WINNER";
        std::set<std::string> genotypes;
        for (const std::size_t g : r.ranking.winner_set) {
            genotypes.insert(genotype_of(r.groups[g][r.representative[g]]->key));
        }
        return genotypes.size() == 1 ? "TIE_PHASE_EQUIVALENT" : "TIE_DISTINCT";
    };

    if (multi_start && (hold_cycles != 0 || resuming || !use_delta)) {
        throw std::runtime_error("--multi-start needs --delta 1 and excludes --hold-cycles and --start");
    }
    const auto run_rounds = [&]() -> int {
        panvar::RepairHooks hooks;
        if (multi_start) {
            hooks.on_scored = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b,
                                  const std::string& sa, const std::string& sb, const panvar::DiplotypeScore& score,
                                  double prior) {
                pool_add(a, b, sa, sb, score, prior, "seed" + std::to_string(current_seed) + "_visit");
            };
        }
        hooks.on_round = [&](const Round& rec) {
            commit();
            check_seed_caps();
            if (interrupt_after != 0 && search.rounds().size() >= interrupt_after && rec.moved && rec.index < last_round) {
                stop_now("TIME_LIMIT",
                         "interrupt-after-rounds " + std::to_string(interrupt_after) +
                             " (test hook: the time-limit path, triggered deterministically)",
                         124);
            }
        };
        if (!search.run(hooks)) {
            const std::size_t round = search.failed_round();
            summary << "refused\ttrue\nrefusal\tthe incumbent did not reproduce in round " << round
                    << "\nincumbent_expected\t" << search.incumbent_expected() << "\nincumbent_scored\t"
                    << search.incumbent_scored() << '\n';
            std::cout << "repair-search " << donor << ": REFUSED, the incumbent did not reproduce in round "
                      << round << " (expected " << search.incumbent_expected() << ", scored "
                      << search.incumbent_scored() << ")\n";
            finish("FAILED", "the incumbent did not reproduce in round " + std::to_string(round));
            return 4;
        }
        return 0;
    };
    std::vector<SeedOutcome> outcomes;
    std::vector<Seed> seeds;
    std::string seed_set_text;
    panvar::FragmentScoringCounters union_context_counters;
    std::size_t union_haplotypes = 0;
    double union_max_recorded_difference = 0.0;
    std::size_t union_recorded_compared = 0;
    try {
        if (!multi_start) {
            const int code = run_rounds();
            if (code != 0) return code;
        } else {
// ---- 1. the union, rescored exactly in this context. Truth is not in scope.
            const auto new_context = [&]() { search.new_context(); };
            std::vector<std::size_t> input_order(union_ids.size());
            std::iota(input_order.begin(), input_order.end(), std::size_t{0});
            if (shuffle_union != 0) {
                std::mt19937_64 rng(shuffle_union);
                std::shuffle(input_order.begin(), input_order.end(), rng);
            }
            std::map<std::vector<std::size_t>, std::size_t> spelled_index;
            std::vector<std::string> spelled;
            const auto spell_of = [&](const std::vector<std::size_t>& path) {
                const auto found = spelled_index.find(path);
                if (found != spelled_index.end()) return found->second;
                spelled.push_back(panvar::spell_mosaic(model, panvar::MosaicPath{path}).bases);
                spelled_index.emplace(path, spelled.size() - 1);
                return spelled.size() - 1;
            };
            for (const std::size_t i : input_order) {
                const Original& o = original_of_id.at(union_ids[i]);
                spell_of(o.first);
                spell_of(o.second);
            }
            union_haplotypes = spelled.size();
            search.open_transient_round();
            const std::uint64_t union_before = search.context()->cumulative_counters().physical_origins;
            const std::vector<std::size_t> slot_of_spelled = search.context()->add_haplotypes(spelled);
            std::map<std::size_t, std::size_t> spelled_of_slot;
            for (std::size_t i = 0; i < slot_of_spelled.size(); ++i) spelled_of_slot.emplace(slot_of_spelled[i], i);
            std::map<std::pair<std::size_t, std::size_t>, std::size_t> union_pair_index;
            std::vector<panvar::HaplotypeCandidate> union_pairs;
            for (const std::size_t i : input_order) {
                const Original& o = original_of_id.at(union_ids[i]);
                std::size_t x = slot_of_spelled[spell_of(o.first)];
                std::size_t y = slot_of_spelled[spell_of(o.second)];
                if (y < x) std::swap(x, y);
                if (union_pair_index.emplace(std::make_pair(x, y), union_pairs.size()).second) union_pairs.push_back({x, y});
            }
            const panvar::FragmentScoringResult union_scored = search.context()->score_candidates(union_pairs);
            search.add_seeding_origins(search.context()->cumulative_counters().physical_origins - union_before);
            for (const std::size_t i : input_order) {
                const Original& o = original_of_id.at(union_ids[i]);
                std::size_t x = slot_of_spelled[spell_of(o.first)];
                std::size_t y = slot_of_spelled[spell_of(o.second)];
                if (y < x) std::swap(x, y);
                const std::size_t at = union_pair_index.at({x, y});
                const std::string& sx = spelled[spelled_of_slot.at(x)];
                const std::string& sy = spelled[spelled_of_slot.at(y)];
                const double prior = panvar::log_unordered_sequence_pair_prior(model, sx, sy, switches);
                PoolEntry& e = pool_add(o.first, o.second, spelled[spell_of(o.first)], spelled[spell_of(o.second)],
                                        union_scored.scores[at], prior, "union_" + o.source);
                e.union_member = true;
                e.union_id = union_ids[i];
                e.union_source = o.source;
                e.union_order = i;
            }
            // The recorded call's own scores, at their printed precision, against the exact rescoring.
            for (const auto& row : read_table(run + ".mosaic_scores.tsv")) {
                const auto found = original_of_id.find(row.at("representative_beam_candidate"));
                if (found == original_of_id.end()) continue;
                const auto entry = pool.find(key_of(found->second.first, found->second.second));
                if (entry == pool.end()) continue;
                const double mine = entry->second.score.log_likelihood + entry->second.prior;
                union_max_recorded_difference =
                    std::max(union_max_recorded_difference, std::fabs(mine - std::stod(row.at("score"))));
                ++union_recorded_compared;
            }

            // ---- 2. the seed set, from the exact union ranking alone
            union_context_counters = search.context()->cumulative_counters();
            const PooledRanking initial = rank_pool(true);
            std::vector<std::size_t> position(initial.groups.size(), 0);
            for (std::size_t p = 0; p < initial.ranking.candidates.size(); ++p) {
                position[initial.ranking.candidates[p].candidate_id] = p;
            }
            const auto group_has = [&](std::size_t g, const std::string& source) {
                for (const PoolEntry* e : initial.groups[g]) {
                    if (e->union_source == source) return true;
                }
                return false;
            };
            const auto tied = [&](const panvar::CertifiedCandidate& a, const panvar::CertifiedCandidate& b) {
                return std::fabs(a.score - b.score) <= a.numeric_error_bound + b.numeric_error_bound;
            };
            std::map<std::size_t, std::set<std::string>> chosen;   // group -> why
            for (const std::size_t g : initial.ranking.winner_set) chosen[g].insert("top_set");
            for (const std::string source : {"panel", "mosaic"}) {
                const panvar::CertifiedCandidate* best = nullptr;
                for (const panvar::CertifiedCandidate& c : initial.ranking.candidates) {
                    if (group_has(c.candidate_id, source)) {
                        best = &c;
                        break;
                    }
                }
                if (best == nullptr) continue;
                for (const panvar::CertifiedCandidate& c : initial.ranking.candidates) {
                    if (tied(c, *best)) chosen[c.candidate_id].insert("best_" + source);
                }
            }
            std::set<std::string> explicit_ids;
            {
                std::stringstream list(explicit_seeds);
                for (std::string id; std::getline(list, id, ',');) {
                    if (!id.empty()) explicit_ids.insert(id);
                }
            }
            std::vector<Seed> gathered;
            std::set<Identity> seen_seed;
            const auto take = [&](std::size_t g, const PoolEntry* e, const std::set<std::string>& why) {
                if (!seen_seed.insert(e->key).second) return;
                Seed seed;
                seed.key = e->key;
                seed.candidate_id = e->union_id;
                seed.union_order = e->union_order;
                seed.union_source = e->union_source;
                seed.group_rank = position[g] + 1;
                const bool panel = group_has(g, "panel");
                const bool mosaic = group_has(g, "mosaic");
                seed.group_sources = panel && mosaic ? "both" : panel ? "panel" : "mosaic";
                seed.selected_as = why;
                seed.certified = initial.ranking.candidates[position[g]];
                gathered.push_back(seed);
            };
            if (!explicit_ids.empty()) {
                for (std::size_t g = 0; g < initial.groups.size(); ++g) {
                    for (const PoolEntry* e : initial.groups[g]) {
                        if (explicit_ids.count(e->union_id) != 0) take(g, e, {"explicit"});
                    }
                }
            } else {
                for (const auto& item : chosen) {
                    for (const PoolEntry* e : initial.groups[item.first]) take(item.first, e, item.second);
                }
            }
            // Canonical order: the exact ranking's total order, then the call's own candidate order.
            std::sort(gathered.begin(), gathered.end(), [](const Seed& x, const Seed& y) {
                return std::tie(x.group_rank, x.union_order) < std::tie(y.group_rank, y.union_order);
            });
            for (std::size_t i = 0; i < gathered.size(); ++i) gathered[i].index = i;
            seeds = gathered;
            {
                std::ostringstream text;
                text << std::setprecision(17)
                     << "seed\tcandidate_id\tunion_order\tunion_source\tselected_as\tgroup_rank\tgroup_sources\t"
                        "score\tfragment_term\texposure_term\tprior\tnumeric_error_bound\tunordered_genotype\t"
                        "hap1_alleles\thap2_alleles\n";
                for (const Seed& seed : seeds) {
                    std::string why;
                    for (const std::string& w : seed.selected_as) why += (why.empty() ? "" : ",") + w;
                    text << seed.index << '\t' << seed.candidate_id << '\t' << seed.union_order << '\t'
                         << seed.union_source << '\t' << why << '\t' << seed.group_rank << '\t' << seed.group_sources
                         << '\t' << seed.certified.score << '\t' << seed.certified.fragment_term << '\t'
                         << seed.certified.exposure_term << '\t' << seed.certified.prior << '\t'
                         << seed.certified.numeric_error_bound << '\t' << genotype_of(seed.key) << '\t'
                         << join(seed.key.first) << '\t' << join(seed.key.second) << '\n';
                }
                seed_set_text = text.str();
                // Partial until the search completes: a refused run leaves no final-named table.
                atomic_write(out + ".seed_set.partial.tsv", seed_set_text);
            }
            if (seeds.size() > max_repair_starts) {
                std::cout << "repair-search " << donor << ": REFUSED before repair, " << seeds.size()
                          << " seeds exceed max_repair_starts " << max_repair_starts << " (never truncated)\n";
                finish("RESOURCE_REFUSAL", std::to_string(seeds.size()) + " repair seeds exceed max_repair_starts " +
                                               std::to_string(max_repair_starts) + "; the set is not truncated");
                return 5;
            }
            if (!seeds_only) {
                // The union was scored in a context of its own. Its 250-odd full haplotypes leave slots,
                // a full emission cache and allocator residue that no seed needs, so the seeds get a
                // fresh shared context (or one each, isolated). Contexts never change a score.
                if (!isolate_seeds) new_context();
                // ---- 3. one independent search per seed, in canonical (or reversed) order
                std::vector<std::size_t> order(seeds.size());
                std::iota(order.begin(), order.end(), std::size_t{0});
                if (reverse_seeds) std::reverse(order.begin(), order.end());
                outcomes.resize(seeds.size());
                for (const std::size_t k : order) {
                    const Seed& seed = seeds[k];
                    if (isolate_seeds) new_context();
                    current_seed = k;
                    seed_start = search.context()->cumulative_counters();
                    rounds_prefix = out + ".seed" + std::to_string(k);
                    search.set_round_kind("seed" + std::to_string(k));
                    search.start(seed.key.first, seed.key.second, seed.certified.score);
                    const auto seed_mark = std::chrono::steady_clock::now();
                    const int code = run_rounds();
                    if (code != 0) return code;
                    check_seed_caps();
                    SeedOutcome& o = outcomes[k];
                    o.seed = k;
                    o.rounds = search.rounds();
                    o.stop = search.stop_reason();
                    o.optimum = key_of(search.incumbent_first(), search.incumbent_second());
                    o.optimum_total = search.incumbent_total();
                    o.incumbents.push_back(seed.key);
                    for (const Round& r : search.rounds()) {
                        if (r.moved) o.incumbents.push_back(key_of(r.first, r.second));
                    }
                    const panvar::FragmentScoringCounters& c = search.context()->cumulative_counters();
                    o.used.physical_origins = c.physical_origins - seed_start.physical_origins;
                    o.used.statistical_state_classes = c.statistical_state_classes - seed_start.statistical_state_classes;
                    o.used.piece_hits = c.piece_hits - seed_start.piece_hits;
                    o.used.latent_alignment_pairs = c.latent_alignment_pairs - seed_start.latent_alignment_pairs;
                    o.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - seed_mark).count();
                    pool.at(o.optimum).tags.insert("seed" + std::to_string(k) + "_optimum");
                    earlier_seed_rounds.insert(earlier_seed_rounds.end(), search.rounds().begin(), search.rounds().end());
                    if (!isolate_seeds) search.context()->release_transient({});
                }
            }
        }
    } catch (const panvar::ResourceExhausted& refusal) {
        finish("RESOURCE_REFUSAL", refusal.what());
        std::cout << "repair-search " << donor << ": RESOURCE_REFUSAL after " << search.rounds().size()
                  << " committed rounds: " << refusal.what() << '\n';
        return 5;
    } catch (const std::exception& failure) {
        finish("FAILED", failure.what());
        throw;
    }

    if (multi_start) {
        // ---- the multi-start report. Truth enters only here, after every search has finished.
        const PooledRanking initial = rank_pool(true);
        const PooledRanking final_ranking = rank_pool(false);
        const std::string initial_class = classify(initial);
        const std::string final_class = classify(final_ranking);
        const auto distance_of = [&](const Identity& k) -> long {
            if (truth_gfa.empty()) return -1;
            ensure_truth();
            long total = 0;
            for (std::size_t b = 0; b < blocks; ++b) {
                if (!truth[0].blocks[b].has_value() || !truth[1].blocks[b].has_value()) continue;
                const std::string& t1 = *truth[0].blocks[b];
                const std::string& t2 = *truth[1].blocks[b];
                const std::string& a1 = model.blocks()[b].alleles[k.first[b]].sequence;
                const std::string& a2 = model.blocks()[b].alleles[k.second[b]].sequence;
                total += std::min(edit_distance(t1, a1) + edit_distance(t2, a2),
                                  edit_distance(t1, a2) + edit_distance(t2, a1));
            }
            return total;
        };
        bool round_cap = false;
        for (const SeedOutcome& o : outcomes) round_cap = round_cap || o.stop == "round_cap";
        const std::string call_status =
            seeds_only ? "NOT_SEARCHED"
            : round_cap ? "NO_CALL_ROUND_CAP"
            : final_class == "SEPARATED_WINNER" ? "CALL_ELIGIBLE"
            : final_class == "TIE_PHASE_EQUIVALENT" ? "NO_CALL_TIE_PHASE_EQUIVALENT"
                                                    : "NO_CALL_TIE_DISTINCT";
        std::set<Identity> optima;
        for (const SeedOutcome& o : outcomes) optima.insert(o.optimum);
        std::ostringstream merges;
        for (std::size_t i = 0; i < outcomes.size(); ++i) {
            for (std::size_t j = i + 1; j < outcomes.size(); ++j) {
                std::size_t shared = 0;
                for (const Identity& a : outcomes[i].incumbents) {
                    shared += std::count(outcomes[j].incumbents.begin(), outcomes[j].incumbents.end(), a);
                }
                if (shared == 0) continue;
                merges << (merges.tellp() > 0 ? ";" : "") << i << "-" << j << ":" << shared
                       << (outcomes[i].optimum == outcomes[j].optimum ? "(same optimum)" : "");
            }
        }
        const auto render_ranking = [&](const PooledRanking& r, std::size_t limit) {
            std::ostringstream t;
            t << std::setprecision(17)
              << "rank\tgroup\tscore\tfragment_term\texposure_term\tprior\tnumeric_error_bound\tin_winner_set\t"
                 "members\tunion_members\ttags\trep_hap1\trep_hap2\tunordered_genotype\tdistance\n";
            const std::set<std::size_t> winners(r.ranking.winner_set.begin(), r.ranking.winner_set.end());
            for (std::size_t p = 0; p < r.ranking.candidates.size(); ++p) {
                const panvar::CertifiedCandidate& c = r.ranking.candidates[p];
                if (p >= limit && winners.count(c.candidate_id) == 0) continue;
                const std::vector<const PoolEntry*>& members = r.groups[c.candidate_id];
                const Identity& rep = members[r.representative[c.candidate_id]]->key;
                std::set<std::string> tags;
                std::size_t union_members = 0;
                for (const PoolEntry* e : members) {
                    tags.insert(e->tags.begin(), e->tags.end());
                    union_members += e->union_member ? 1 : 0;
                }
                std::string joined;
                for (const std::string& tag : tags) joined += (joined.empty() ? "" : ",") + tag;
                t << (p + 1) << '\t' << c.candidate_id << '\t' << c.score << '\t' << c.fragment_term << '\t'
                  << c.exposure_term << '\t' << c.prior << '\t' << c.numeric_error_bound << '\t'
                  << (winners.count(c.candidate_id) ? "true" : "false") << '\t' << members.size() << '\t'
                  << union_members << '\t' << joined << '\t' << join(rep.first) << '\t' << join(rep.second) << '\t'
                  << genotype_of(rep) << '\t' << distance_of(rep) << '\n';
            }
            return t.str();
        };
        const auto winners_text = [&](const PooledRanking& r) {
            std::string text;
            for (const std::size_t g : r.ranking.winner_set) {
                const Identity& rep = r.groups[g][r.representative[g]]->key;
                text += (text.empty() ? "" : ";") + join(rep.first) + "|" + join(rep.second);
            }
            return text;
        };
        std::ostringstream seed_table;
        seed_table << std::setprecision(17)
                   << "seed\tcandidate_id\tselected_as\tgroup_rank\tgroup_sources\tunion_source\tstart_score\t"
                      "unordered_genotype\thap1_alleles\thap2_alleles\trounds\tmoves\tstop_reason\toptimum_total\t"
                      "optimum_hap1\toptimum_hap2\toptimum_genotype\tstart_distance\toptimum_distance\tfloor\t"
                      "optimum_excess\tseconds\torigins_used\tclasses_used\tpiece_hits_used\tlatent_pairs_used\n";
        for (const SeedOutcome& o : outcomes) {
            const Seed& seed = seeds[o.seed];
            std::string why;
            for (const std::string& w : seed.selected_as) why += (why.empty() ? "" : ",") + w;
            const long d_start = distance_of(seed.key);
            const long d_opt = distance_of(o.optimum);
            seed_table << seed.index << '\t' << seed.candidate_id << '\t' << why << '\t' << seed.group_rank << '\t'
                       << seed.group_sources << '\t' << seed.union_source << '\t' << seed.certified.score << '\t'
                       << genotype_of(seed.key) << '\t' << join(seed.key.first) << '\t' << join(seed.key.second)
                       << '\t' << o.rounds.size() << '\t' << (o.incumbents.size() - 1) << '\t' << o.stop << '\t'
                       << o.optimum_total << '\t' << join(o.optimum.first) << '\t' << join(o.optimum.second) << '\t'
                       << genotype_of(o.optimum) << '\t' << d_start << '\t' << d_opt << '\t'
                       << (truth_gfa.empty() ? -1 : floor_total) << '\t'
                       << (truth_gfa.empty() ? -1 : d_opt - floor_total) << '\t' << o.seconds << '\t'
                       << o.used.physical_origins << '\t' << o.used.statistical_state_classes << '\t'
                       << o.used.piece_hits << '\t' << o.used.latent_alignment_pairs << '\n';
        }
        std::ostringstream report;
        report << std::setprecision(17) << "multi_start\ttrue\nseeds\t" << seeds.size() << "\nmax_repair_starts\t"
               << max_repair_starts << "\nseed_order\t" << (reverse_seeds ? "reversed" : "canonical")
               << "\nseed_contexts\t" << (isolate_seeds ? "isolated" : "shared") << "\nunion_candidates\t"
               << union_ids.size() << "\nunion_haplotypes\t" << union_haplotypes << "\nunion_recorded_compared\t"
               << union_recorded_compared << "\nunion_max_recorded_difference\t" << union_max_recorded_difference
               << "\npool_identities\t" << pool.size() << "\npool_score_mismatches\t" << pool_mismatches
               << "\ninitial_sequence_candidates\t" << initial.groups.size() << "\ninitial_class\t" << initial_class
               << "\ninitial_winner_set\t" << initial.ranking.winner_set.size() << "\ninitial_certified_gap\t"
               << initial.ranking.certified_gap << "\ninitial_winners\t" << winners_text(initial)
               << "\nfinal_sequence_candidates\t" << final_ranking.groups.size() << "\nfinal_class\t" << final_class
               << "\nfinal_winner_set\t" << final_ranking.ranking.winner_set.size() << "\nfinal_certified_gap\t"
               << final_ranking.ranking.certified_gap << "\nfinal_best_score\t"
               << (final_ranking.ranking.candidates.empty() ? 0.0 : final_ranking.ranking.candidates.front().score)
               << "\nfinal_winners\t" << winners_text(final_ranking) << "\ncall_status\t" << call_status
               << "\ndistinct_optima\t" << optima.size() << "\ntrajectory_merges\t"
               << (merges.tellp() > 0 ? merges.str() : std::string("none")) << '\n';
        if (!truth_gfa.empty() && !final_ranking.ranking.winner_set.empty()) {
            long best = std::numeric_limits<long>::max();
            long worst = 0;
            for (const std::size_t g : final_ranking.ranking.winner_set) {
                const long d = distance_of(final_ranking.groups[g][final_ranking.representative[g]]->key);
                best = std::min(best, d);
                worst = std::max(worst, d);
            }
            report << "final_winner_distance_min\t" << best << "\nfinal_winner_distance_max\t" << worst
                   << "\navailability_floor_total\t" << floor_total << "\nfinal_winner_excess_min\t"
                   << (best - floor_total) << '\n';
        }
        const panvar::FragmentScoringCounters& cumulative = search.context()->cumulative_counters();
        report << "union_context_physical_origins\t" << union_context_counters.physical_origins
               << "\nunion_context_state_classes\t" << union_context_counters.statistical_state_classes << '\n';
        report << "cumulative_physical_origins\t" << cumulative.physical_origins << "\ncumulative_state_classes\t"
               << cumulative.statistical_state_classes << "\ncumulative_piece_hits\t" << cumulative.piece_hits
               << "\ncumulative_latent_pairs\t" << cumulative.latent_alignment_pairs << "\nprocess_margin_origins\t"
               << margin(scoring.max_physical_origins, cumulative.physical_origins) << "\nprocess_margin_classes\t"
               << margin(scoring.max_state_classes, cumulative.statistical_state_classes)
               << "\nthreads_requested\t" << scoring.threads << "\npeak_active_workers\t"
               << panvar::peak_active_workers() << "\nphys_footprint_peak_bytes\t" << physical_footprint().peak << '\n';
        std::vector<panvar::OriginDemandRow> demand_rows;
        if (search.context() && use_origin_cache) demand_rows = search.context()->origin_demand();
        {
            std::lock_guard<std::mutex> lock(commit_mutex);
            if (finished) return 6;
            finished = true;
            for (const SeedOutcome& o : outcomes) {
                atomic_write(out + ".seed" + std::to_string(o.seed) + ".rounds.tsv", search.render_rounds(o.rounds));
                std::error_code ignored;
                std::filesystem::remove(out + ".seed" + std::to_string(o.seed) + ".rounds.partial.tsv", ignored);
            }
            atomic_write(out + ".seeds.tsv", seed_table.str());
            atomic_write(out + ".seed_set.tsv", seed_set_text);
            {
                std::error_code ignored;
                std::filesystem::remove(out + ".seed_set.partial.tsv", ignored);
            }
            atomic_write(out + ".initial_ranking.tsv", render_ranking(initial, 50));
            atomic_write(out + ".final_ranking.tsv", render_ranking(final_ranking, 50));
            if (search.memory_on()) atomic_write(out + ".memory.tsv", search.memory_text());
            if (!demand_rows.empty()) atomic_write(out + ".origin_demand.tsv", render_demand(demand_rows));
            atomic_write(out + ".summary.tsv", summary.str() + report.str());
            for (const char* partial : {"rounds", "memory", "origin_demand", "distance", "round_truth"}) {
                std::error_code ignored;
                std::filesystem::remove(out + "." + partial + ".partial.tsv", ignored);
            }
            status_extra = report.str();
            write_status("COMPLETE", seeds_only ? "seed set only" : call_status);
        }
        std::cout << "repair-search " << donor << ": multi-start, " << seeds.size() << " seeds, "
                  << optima.size() << " distinct optima, initial " << initial_class << " -> final " << final_class
                  << ", " << call_status << '\n';
        return 0;
    }

    bool order_all = true;
    bool swap_all = true;
    for (const Round& r : search.rounds()) {
        order_all = order_all && r.order_invariant;
        swap_all = swap_all && r.swap_invariant;
    }
    {
        double placement = 0.0;
        double combination = 0.0;
        std::size_t new_haps = 0;
        std::size_t reused_haps = 0;
        for (const Round& r : search.rounds()) {
            placement += r.placement_seconds;
            combination += r.combination_seconds;
            new_haps += r.haps_new;
            reused_haps += r.haps_reused;
        }
        std::uint64_t origins = 0;
        std::uint64_t classes = 0;
        for (const Round& r : search.rounds()) {
            origins += r.physical_origins;
            classes += r.state_classes;
        }
        if (search.context()) {
            origins = search.context()->cumulative_counters().physical_origins;
            classes = search.context()->cumulative_counters().statistical_state_classes;
        }
        summary << "cumulative_physical_origins\t" << origins << "\ncumulative_state_classes\t"
                << classes << "\ndeclared_max_physical_origins\t" << scoring.max_physical_origins
                << '\n';
        std::uint64_t reused = 0;
        std::uint64_t recomputed = 0;
        for (const Round& r : search.rounds()) {
            reused += r.origins_reused;
            recomputed += r.origins_recomputed;
        }
        summary << "origins_reused\t" << reused << "\norigins_recomputed\t" << recomputed << '\n';
        summary << "scoring_mode\t"
                << (use_delta ? "block_delta" : (search.context() ? "prepared_context" : "one_shot_per_round"))
                << "\nhaplotypes_placed\t" << new_haps << "\nhaplotypes_reused\t" << reused_haps
                << "\nseconds_placement\t" << placement << "\nseconds_combination\t" << combination
                << "\nfraction_scanning_new_haplotypes\t"
                << (search.cumulative_seconds() > 0.0 ? placement / search.cumulative_seconds() : 0.0) << '\n';
        if (search.context()) {
            summary << "context_unique_haplotypes\t" << search.context()->unique_haplotypes() << '\n';
        }
    }
    summary << "rounds_run\t" << search.rounds().size() << "\nstop_reason\t" << search.stop_reason()
            << "\nvisited_size\t" << search.visited_size() << "\norder_invariant\t"
            << (order_all ? "true" : "false") << "\nswap_invariant\t" << (swap_all ? "true" : "false")
            << "\ncumulative_scoring_seconds\t" << search.cumulative_seconds() << "\ncumulative_scored_pairs\t"
            << search.cumulative_pairs() << "\npeak_rss_gib\t" << peak_rss_gib() << "\nfinal_total\t" << search.incumbent_total()
            << "\nimprovement_total\t" << (search.incumbent_total() - start_total) << '\n';
    std::vector<panvar::OriginDemandRow> demand_rows;
    if (search.context() && use_origin_cache) demand_rows = search.context()->origin_demand();
    summary << demand_report(demand_rows);
    summary << "phys_footprint_peak_bytes\t" << physical_footprint().peak << '\n';
    // Threads: what was requested, what the window enumeration split into, and the most workers that
    // were ever active at once across every stage -- which may never exceed the request.
    summary << "threads_requested\t" << scoring.threads << "\npeak_active_workers\t"
            << panvar::peak_active_workers() << '\n';
    if (search.context()) {
        const panvar::ParallelEnumerationStats p = search.context()->parallel_enumeration_stats();
        summary << "enumeration_workers\t" << p.workers << "\nparallel_enumerations\t" << p.parallel_enumerations
                << "\nparallel_duplicate_dps\t" << p.duplicate_dps << "\nparallel_duplicate_scan_bases\t"
                << p.duplicate_scan_bases << '\n';
    }

    // ---- completion. Truth, if any, is loaded first (a committed round usually loaded it already);
    // then every final file is written under the commit lock, so a watchdog can neither interleave
    // with them nor stamp a finished run as stopped. The rounds table is TRUTH-FREE, as it always was.
    TruthTables tables;
    if (!truth_gfa.empty()) {
        ensure_truth();
        tables = render_truth();
    }
    {
        std::lock_guard<std::mutex> lock(commit_mutex);
        if (finished) return 6;   // a stop already won; its status stands
        finished = true;
        rounds_committed = search.rounds().size();
        atomic_write(out + ".rounds.tsv", render_rounds());
        if (search.memory_on()) atomic_write(out + ".memory.tsv", search.memory_text());
        if (!demand_rows.empty()) atomic_write(out + ".origin_demand.tsv", render_demand(demand_rows));
        if (truth_gfa.empty()) {
            summary << "truth_evaluated\tfalse\n";
        } else {
            atomic_write(out + ".distance.tsv", tables.distance);
            atomic_write(out + ".round_truth.tsv", tables.per_block);
            summary << "truth_evaluated\ttrue\navailability_floor_total\t" << floor_total
                    << "\nstart_distance\t" << tables.start_distance << "\nfinal_distance\t"
                    << tables.final_distance << "\nfinal_excess\t" << (tables.final_distance - floor_total)
                    << '\n';
        }
        atomic_write(out + ".summary.tsv", summary.str());
        for (const char* partial : {"rounds", "memory", "origin_demand", "distance", "round_truth"}) {
            std::error_code ignored;
            std::filesystem::remove(out + "." + partial + ".partial.tsv", ignored);
        }
        status_extra = demand_report(demand_rows);
        if (piece_screen_shadow) {
            panvar::set_piece_screen_shadow(false);
            const panvar::PieceScreenShadowReport report = panvar::piece_screen_shadow_report();
            std::ostringstream text;
            text << "pieces\tentering\tcount_pass\tadmitted\trejected\trejected_prefilter_passed\t"
                    "rejected_accepted\tscreen_ns\trejected_emission_ns\trejected_prefilter_ns\t"
                    "band\tstarts\tall_entering\tprefilter_passed\taccepted\temission_ns\tprefilter_ns\t"
                    "code_build_ns\tlayout_ns\n";
            for (const panvar::PieceScreenShadowRow& row : report.rows) {
                text << row.pieces << '\t' << row.entering << '\t' << row.count_pass << '\t' << row.admitted
                     << '\t' << row.rejected << '\t' << row.rejected_prefilter_passed << '\t'
                     << row.rejected_accepted << '\t' << row.screen_ns << '\t' << row.rejected_emission_ns
                     << '\t' << row.rejected_prefilter_ns << '\t' << report.band << '\t' << report.starts
                     << '\t' << report.entering << '\t' << report.prefilter_passed << '\t' << report.accepted
                     << '\t' << report.emission_ns << '\t' << report.prefilter_ns << '\t'
                     << report.code_build_ns << '\t' << report.layout_ns << '\n';
            }
            atomic_write(out + ".piece_screen.tsv", text.str());
        }
        for (const auto& [label, hash] : export_windows) {
            const panvar::DemandWindow w = search.context()->demand_window(hash);
            std::ostringstream text;
            text << "key\tvalue\nlabel\t" << label << "\nkey_hash\t" << hash << "\nfound\t"
                 << (w.found ? "true" : "false") << "\nat_start\t" << (w.at_start ? 1 : 0) << "\nat_end\t"
                 << (w.at_end ? 1 : 0) << "\nrel_block_begin\t" << w.rel_block_begin << "\nrel_new_end\t"
                 << w.rel_new_end << "\nwindow_bytes\t" << w.bytes.size() << "\nbytes\t" << w.bytes << '\n';
            atomic_write(out + ".window." + label + ".tsv", text.str());
        }
        write_status("COMPLETE", search.stop_reason());
    }
    if (truth_gfa.empty()) {
        std::cout << "repair-search " << donor << ": " << search.rounds().size() << " rounds, stop=" << search.stop_reason()
                  << ", improvement " << (search.incumbent_total() - start_total) << " nats, no truth evaluation\n";
        return 0;
    }
    std::cout << "repair-search " << donor << ": " << search.rounds().size() << " rounds, stop=" << search.stop_reason()
              << ", improvement " << (search.incumbent_total() - start_total) << " nats, distance " << tables.start_distance
              << " -> " << tables.final_distance << " (floor " << floor_total << "), " << search.cumulative_seconds()
              << "s scoring\n";
    return 0;
}

// ------------------------------------------------------------------------------------------------
// window-bench: MEASUREMENT ONLY. One recorded window's origin enumeration, with its fragments divided
// among --workers contiguous ranges, checked bit-identical against the production path.
int window_bench(const Args& args) {
    const std::string out = args.required("out");
    const auto status = read_key_values(args.required("run") + ".mosaic_status.tsv");
    panvar::FragmentScoringOptions scoring;
    scoring.threads = 1;
    scoring.mismatch_band = std::stoul(status.at("mismatch_band"));
    scoring.minimum_piece_bp = 12;
    scoring.max_piece_hits = std::stoull(status.at("max_piece_hits"));
    scoring.max_proposed_starts_per_haplotype = std::stoull(status.at("max_proposed_starts_per_haplotype"));
    scoring.max_physical_origins = std::stoull(status.at("max_physical_origins"));
    scoring.max_state_classes = std::stoull(status.at("max_state_classes"));
    scoring.max_mass_cells = std::stoull(status.at("max_mass_cells"));
    scoring.emission_cache_bytes = static_cast<std::size_t>(std::stoull(status.at("emission_cache_bytes")));
    const panvar::LibraryModel library = read_library_file(args.required("library"));
    const std::vector<panvar::ReadPair> fragments = read_pairs(args.required("reads1"), args.required("reads2"));
    const auto window = read_key_values(args.required("window"));
    const std::size_t workers = std::stoul(args.get("workers", "1"));
    const std::size_t repeats = std::stoul(args.get("repeats", "5"));
    const FootprintSample loaded = physical_footprint();
    panvar::PreparedScoringContext context(fragments, library, scoring);
    const FootprintSample prepared = physical_footprint();
    const panvar::WindowEnumerationBench b = context.benchmark_window_enumeration(
        window.at("bytes"), window.at("at_start") == "1", window.at("at_end") == "1", workers, repeats);
    const FootprintSample done = physical_footprint();
    std::vector<double> walls = b.wall_seconds;
    std::sort(walls.begin(), walls.end());
    const double median = walls[walls.size() / 2];
    double busy = 0.0;
    for (const double x : b.worker_busy_seconds) busy += x;
    std::ostringstream text;
    text << std::setprecision(9) << "key\tvalue\nlabel\t" << window.at("label") << "\nworkers\t" << b.workers
         << "\nfragments\t" << b.fragments << "\nwindow_bytes\t" << b.window_bytes << "\nat_start\t"
         << window.at("at_start") << "\nat_end\t" << window.at("at_end") << "\nrepeats\t" << walls.size()
         << "\nreference_cold_seconds\t" << b.reference_cold_seconds << "\nproduction_warm_seconds\t"
         << b.production_warm_seconds << "\nwall_seconds_median\t" << median << "\nwall_seconds_min\t" << walls.front()
         << "\nwall_seconds_max\t" << walls.back() << "\nworker_busy_seconds_sum\t" << busy << "\nutilization\t"
         << (busy / (static_cast<double>(b.workers) * b.wall_seconds.back())) << "\nprefilter_seconds\t"
         << b.prefilter_seconds << "\nemission_seconds\t" << b.emission_seconds << "\norigins_identical\t"
         << (b.origins_identical ? "true" : "false") << "\ncounters_identical\t"
         << (b.counters_identical ? "true" : "false") << "\nfirst_difference\t"
         << (b.first_difference.empty() ? "none" : b.first_difference) << "\norigins\t" << b.origins
         << "\npiece_hits\t" << b.piece_hits << "\nlatent_pairs\t" << b.latent_pairs << "\nexecuted_dps\t"
         << b.executed_dps << "\nexecuted_prefilters\t" << b.executed_prefilters << "\nverifications\t"
         << b.verifications << "\nduplicate_scan_bases\t" << b.duplicate_scan_bases
         << "\nphys_footprint_loaded\t" << loaded.current << "\nphys_footprint_prepared\t" << prepared.current
         << "\nphys_footprint_done\t" << done.current << "\nphys_footprint_peak\t" << done.peak << '\n';
    atomic_write(out, text.str());
    std::cout << "window-bench " << window.at("label") << " workers " << b.workers << ": median " << median
              << " s (production warm " << b.production_warm_seconds << " s), origins "
              << (b.origins_identical ? "identical" : "DIFFER") << ", counters "
              << (b.counters_identical ? "identical" : "DIFFER") << '\n';
    return b.origins_identical && b.counters_identical ? 0 : 7;
}

// ------------------------------------------------------------------------------------------------
// score-set: exactly score a handful of declared diplotypes, and nothing else
//
// A candidate's score depends only on its own two haplotypes, the fragments and the library -- never
// on which other candidates are scored beside it. That is what makes this legitimate: six candidates
// scored alone carry the same fragment and exposure terms they carried inside a 9,028-candidate
// union, so a compound candidate can be tested without rebuilding that union.
//
// The declared file is <label> <hap1 alleles> <hap2 alleles> [expected total]. Where an expected
// total is given it MUST reproduce, and the command refuses before reporting anything if it does
// not: a compound result is only interpretable once the candidates it is compared against are shown
// to be the same ones that were scored before.
//
// Truth is loaded after every candidate is scored, and only to measure distance.
int score_set(const Args& args) {
    const std::string out = args.required("out");
    struct StagedGraphs {
        std::string prefix;
        ~StagedGraphs() {
            std::error_code ignored;
            std::filesystem::remove(prefix + ".model.sorted.gfa", ignored);
            std::filesystem::remove(prefix + ".truth.sorted.gfa", ignored);
        }
    } staged_graphs{out};
    const std::string reference = args.required("reference");
    const std::string donor = args.required("donor");
    const std::string run = args.required("run");
    const double switch_probability = std::stod(args.get("switch-probability", "0.01"));

    const LoadedModel loaded =
        load_model(args.required("model-gfa"), reference, out + ".model.sorted.gfa", 50);
    const panvar::PreparedGenotypeModel& model = *loaded.model;
    const std::size_t blocks = model.blocks().size();
    const auto status = read_key_values(run + ".mosaic_status.tsv");
    if (status.at("model_fingerprint") != loaded.fingerprint) {
        throw std::runtime_error("the run was not scored on the rebuilt model");
    }

    struct Declared {
        std::string label;
        std::vector<std::size_t> first;
        std::vector<std::size_t> second;
        bool has_expected = false;
        double expected = 0.0;
    };
    std::vector<Declared> declared;
    {
        std::ifstream file(args.required("candidates"));
        if (!file) throw std::runtime_error("cannot read " + args.required("candidates"));
        for (std::string line; std::getline(file, line);) {
            if (line.empty() || line[0] == '#') continue;
            std::stringstream fields(line);
            std::string label, first, second, expected;
            if (!std::getline(fields, label, '\t') || !std::getline(fields, first, '\t') ||
                !std::getline(fields, second, '\t')) {
                throw std::runtime_error("each candidate needs <label> <hap1> <hap2> [expected total]");
            }
            Declared entry;
            entry.label = label;
            entry.first = parse_alleles(first);
            entry.second = parse_alleles(second);
            if (std::getline(fields, expected, '\t') && expected != "NA" && !expected.empty()) {
                entry.has_expected = true;
                entry.expected = std::stod(expected);
            }
            if (entry.first.size() != blocks || entry.second.size() != blocks) {
                throw std::runtime_error("candidate " + label + " needs one allele per block");
            }
            declared.push_back(std::move(entry));
        }
    }
    if (declared.empty()) throw std::runtime_error("no candidates declared");

    const panvar::LibraryModel library = read_library_file(args.required("library"));
    const std::vector<panvar::ReadPair> fragments =
        read_pairs(args.required("reads1"), args.required("reads2"));

    // --rank-blocks: every unordered allele pair at the named blocks, ranked by proposal emission exactly as
    // CoordinateRepairSearch ranks its local pairs (emission descending, then indices), so a pair's rank says
    // whether it was in the repair's local-pair list. Truth-free; the command's stage 2 with its defaults.
    const std::string rank_blocks = args.get("rank-blocks", "");
    if (!rank_blocks.empty()) {
        const panvar::PreparedMarkerPanel markers = panvar::build_marker_panel(model, 31, 15);
        const panvar::MarkerObservations observations = panvar::count_marker_observations(markers, fragments);
        std::optional<double> declared_background;
        {
            std::ifstream library_file(args.required("library"));
            std::string key;
            double value = 0.0;
            for (std::string line; std::getline(library_file, line);) {
                std::istringstream fields(line);
                if (fields >> key >> value && key == "marker_background") declared_background = value;
            }
        }
        const panvar::MarkerDepthEstimate depth =
            panvar::estimate_marker_depth(markers, observations, declared_background);
        if (!depth.converged) throw std::runtime_error("the depth fit did not converge: " + depth.termination);
        const panvar::PreparedProposalEmissions emissions = panvar::build_proposal_emissions(
            markers, observations, panvar::proposal_marker_parameters(depth, 31, 15));
        const std::size_t local_pairs = std::stoul(status.at("repair_local_pairs"));
        std::ofstream ranks(out + ".local_ranks.tsv");
        ranks << std::setprecision(12)
              << "block_id\trank\tallele_x\tallele_y\tx_bp\ty_bp\temission\tbelow_top\tin_local_list\n";
        std::stringstream wanted(rank_blocks);
        for (std::string id; std::getline(wanted, id, ',');) {
            std::size_t b = blocks;
            for (std::size_t i = 0; i < blocks; ++i) {
                if (model.blocks()[i].block_id == id) b = i;
            }
            if (b == blocks) throw std::runtime_error("no block " + id);
            const panvar::BlockProposalEmissions& emission = emissions.blocks[b];
            std::vector<std::pair<double, std::pair<std::size_t, std::size_t>>> ranked;
            for (std::size_t x = 0; x < emission.allele_count; ++x) {
                for (std::size_t y = x; y < emission.allele_count; ++y) ranked.push_back({emission.at(x, y), {x, y}});
            }
            std::sort(ranked.begin(), ranked.end(), [](const auto& lhs, const auto& rhs) {
                if (lhs.first != rhs.first) return lhs.first > rhs.first;
                return lhs.second < rhs.second;
            });
            for (std::size_t r = 0; r < ranked.size(); ++r) {
                const auto& alleles = model.blocks()[b].alleles;
                ranks << id << '\t' << (r + 1) << '\t' << ranked[r].second.first << '\t' << ranked[r].second.second
                      << '\t' << alleles[ranked[r].second.first].sequence.size() << '\t'
                      << alleles[ranked[r].second.second].sequence.size() << '\t' << ranked[r].first << '\t' << (ranked.front().first - ranked[r].first) << '\t'
                      << (r < local_pairs ? "true" : "false") << '\n';
            }
        }
        ranks.close();
        std::cout << "score-set " << donor << ": local ranks written for " << rank_blocks << '\n';
        if (args.get("emissions-only", "0") == "1") return 0;
    }

    std::vector<panvar::MosaicPath> paths;
    std::vector<std::string> sequences;
    std::map<std::vector<std::size_t>, std::size_t> slot_of_path;
    const auto slot_for = [&](const std::vector<std::size_t>& alleles) {
        const auto found = slot_of_path.find(alleles);
        if (found != slot_of_path.end()) return found->second;
        const std::size_t slot = paths.size();
        slot_of_path.emplace(alleles, slot);
        paths.push_back(panvar::MosaicPath{alleles});
        sequences.push_back(panvar::spell_mosaic(model, panvar::MosaicPath{alleles}).bases);
        return slot;
    };
    std::vector<std::pair<std::size_t, std::size_t>> slots;
    for (const Declared& entry : declared) slots.push_back({slot_for(entry.first), slot_for(entry.second)});

    const panvar::InternedHaplotypes interned = panvar::intern_haplotypes(paths, sequences);
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> scored_index;
    std::vector<panvar::HaplotypeCandidate> scored_pairs;
    std::vector<std::size_t> pair_of_declared(declared.size(), 0);
    for (std::size_t i = 0; i < declared.size(); ++i) {
        std::size_t x = interned.index_of[slots[i].first];
        std::size_t y = interned.index_of[slots[i].second];
        if (y < x) std::swap(x, y);
        const auto inserted = scored_index.emplace(std::make_pair(x, y), scored_pairs.size());
        if (inserted.second) scored_pairs.push_back({x, y});
        pair_of_declared[i] = inserted.first->second;
    }

    panvar::FragmentScoringOptions scoring;
    scoring.threads = std::stoul(args.get("threads", "8"));
    scoring.mismatch_band = std::stoul(status.at("mismatch_band"));
    scoring.minimum_piece_bp = 12;
    scoring.max_piece_hits = std::stoull(status.at("max_piece_hits"));
    scoring.max_proposed_starts_per_haplotype = std::stoull(status.at("max_proposed_starts_per_haplotype"));
    scoring.max_physical_origins = std::stoull(status.at("max_physical_origins"));
    scoring.max_state_classes = std::stoull(status.at("max_state_classes"));
    scoring.max_mass_cells = std::stoull(status.at("max_mass_cells"));
    scoring.emission_cache_bytes = static_cast<std::size_t>(std::stoull(status.at("emission_cache_bytes")));
    // --decompose A,B: the scorer's own accepted origins for every fragment, recorded (not re-derived).
    const std::string decompose = args.get("decompose", "");
    if (!decompose.empty()) {
        for (std::size_t f = 0; f < fragments.size(); ++f) scoring.record_origins_for_fragments.push_back(f);
    }
    const auto mark = std::chrono::steady_clock::now();
    const panvar::FragmentScoringResult scored =
        panvar::score_diplotypes(interned, scored_pairs, fragments, library, scoring);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - mark).count();

    const std::vector<double> switches(blocks - 1, switch_probability);
    std::vector<double> prior_of_pair(scored_pairs.size(), 0.0);
    for (std::size_t slot = 0; slot < scored_pairs.size(); ++slot) {
        prior_of_pair[slot] = panvar::log_unordered_sequence_pair_prior(
            model, interned.sequences[scored_pairs[slot].first],
            interned.sequences[scored_pairs[slot].second], switches);
    }

    // Reproduction first: a declared expectation that does not reproduce refuses the whole report.
    std::size_t checked = 0;
    std::size_t failed = 0;
    double worst = 0.0;
    for (std::size_t i = 0; i < declared.size(); ++i) {
        if (!declared[i].has_expected) continue;
        const std::size_t slot = pair_of_declared[i];
        const double total = scored.scores[slot].log_likelihood + prior_of_pair[slot];
        const double delta = std::fabs(total - declared[i].expected);
        ++checked;
        if (delta > 1e-6) {
            ++failed;
            worst = std::max(worst, delta);
        }
    }
    if (failed != 0) {
        std::cout << "score-set " << donor << ": REFUSED, " << failed << " of " << checked
                  << " declared expectations did not reproduce (worst delta " << worst << ")\n";
        return 4;
    }

    // The truth-free score table is written, and its digest printed, before truth is loaded.
    {
        std::ostringstream free_table;
        free_table << std::setprecision(17)
                   << "label\tsequence_pair\thap1_alleles\thap2_alleles\tfragment_term\texposure_term\tprior\ttotal\t"
                      "expected_total\n";
        for (std::size_t i = 0; i < declared.size(); ++i) {
            const std::size_t slot = pair_of_declared[i];
            std::string h1, h2;
            for (std::size_t b = 0; b < blocks; ++b) {
                h1 += (b ? "," : "") + std::to_string(declared[i].first[b]);
                h2 += (b ? "," : "") + std::to_string(declared[i].second[b]);
            }
            free_table << declared[i].label << '\t' << slot << '\t' << h1 << '\t' << h2 << '\t'
                       << scored.scores[slot].fragment_term << '\t' << scored.scores[slot].exposure_term << '\t'
                       << prior_of_pair[slot] << '\t' << (scored.scores[slot].log_likelihood + prior_of_pair[slot])
                       << '\t' << (declared[i].has_expected ? std::to_string(declared[i].expected) : std::string("NA"))
                       << '\n';
        }
        atomic_write(out + ".scores_truth_free.tsv", free_table.str());
        std::cout << "score-set " << donor << ": truth-free table sha256 "
                  << panvar::sha256_file(out + ".scores_truth_free.tsv") << '\n';
    }
    // ---- --decompose A,B (truth-free): the fragment-term difference A - B, fragment by fragment, with each
    // haplotype's origins, collapsed states, placement class against --decompose-block, orientation and insert
    // geometry, and each mass cell rebuilt from its recorded origins to check they are counted exactly once.
    if (!decompose.empty()) {
        const std::size_t comma = decompose.find(',');
        const std::string label_a = decompose.substr(0, comma);
        const std::string label_b = decompose.substr(comma + 1);
        std::size_t ia = declared.size(), ib = declared.size();
        for (std::size_t i = 0; i < declared.size(); ++i) {
            if (declared[i].label == label_a) ia = i;
            if (declared[i].label == label_b) ib = i;
        }
        if (ia == declared.size() || ib == declared.size()) throw std::runtime_error("--decompose names an undeclared label");
        const std::string block_id = args.required("decompose-block");
        std::size_t block = blocks;
        for (std::size_t b = 0; b < blocks; ++b) {
            if (model.blocks()[b].block_id == block_id) block = b;
        }
        if (block == blocks) throw std::runtime_error("no block " + block_id);
        struct Hap {
            std::string name;
            std::size_t unique = 0;
            std::size_t block_begin = 0, block_end = 0, length = 0;
        };
        std::vector<Hap> haps;
        for (const auto& [name, i, side] : std::vector<std::tuple<std::string, std::size_t, int>>{
                 {"A1", ia, 0}, {"A2", ia, 1}, {"B1", ib, 0}, {"B2", ib, 1}}) {
            const std::vector<std::size_t>& path = side == 0 ? declared[i].first : declared[i].second;
            const panvar::HaplotypeSequence spelled = panvar::spell_mosaic(model, panvar::MosaicPath{path});
            Hap hap;
            hap.name = name;
            hap.unique = interned.index_of[slot_of_path.at(path)];
            hap.block_begin = spelled.block_intervals[block].begin;
            hap.block_end = spelled.block_intervals[block].end;
            hap.length = spelled.bases.size();
            haps.push_back(hap);
        }
        const panvar::DiscreteInsertDistribution insert = panvar::build_library_insert_distribution(library);
        std::ostringstream rows;
        rows << std::setprecision(17) << "fragment\tlog_background";
        for (const Hap& h : haps) {
            rows << '\t' << h.name << "_log_mass\t" << h.name << "_rebuilt_equal\t" << h.name << "_origins\t" << h.name
                 << "_collapsed_states\t" << h.name << "_exact_duplicates\t" << h.name << "_in_block\t" << h.name
                 << "_junction\t" << h.name << "_outside\t" << h.name << "_first_forward\t" << h.name << "_mean_insert";
        }
        rows << "\tcontribution_A\tcontribution_B\tdelta\tclass\n";
        struct ClassSum {
            std::size_t fragments = 0;
            double delta = 0.0;
            std::vector<double> deltas;
        };
        std::map<std::string, ClassSum> classes;
        std::vector<double> all_deltas;
        std::size_t rebuilt_mismatches = 0, duplicate_origins = 0;
        double sum_a = 0.0, sum_b = 0.0;
        for (std::size_t f = 0; f < fragments.size(); ++f) {
            const double background = panvar::log_background_probability(fragments[f]);
            rows << f << '\t' << background;
            std::vector<std::size_t> in_block(4, 0), junction(4, 0);
            for (std::size_t k = 0; k < haps.size(); ++k) {
                const Hap& h = haps[k];
                const double mass = scored.mass_of(f, h.unique);
                const auto found = scored.recorded_origins.find({f, h.unique});
                const std::vector<panvar::FragmentPlacement> none;
                const std::vector<panvar::FragmentPlacement>& origins =
                    found == scored.recorded_origins.end() ? none : found->second;
                const double rebuilt = origins.empty() ? -std::numeric_limits<double>::infinity()
                                                       : panvar::log_placement_mass(origins, fragments[f], insert, library);
                const bool equal = std::memcmp(&rebuilt, &mass, sizeof(double)) == 0 ||
                                   (std::isinf(rebuilt) && std::isinf(mass));
                if (!equal) ++rebuilt_mismatches;
                std::set<std::tuple<std::size_t, std::size_t, int>> positions;
                std::set<std::tuple<std::size_t, int, std::uint64_t>> states;
                std::size_t duplicates = 0, outside = 0, forward = 0;
                double insert_sum = 0.0;
                for (const panvar::FragmentPlacement& o : origins) {
                    if (!positions.insert({o.forward_start, o.insert_bp, static_cast<int>(o.orientation)}).second) ++duplicates;
                    std::uint64_t bits = 0;
                    std::memcpy(&bits, &o.log_latent_mass, sizeof(bits));
                    states.insert({o.insert_bp, static_cast<int>(o.orientation), bits});
                    const std::size_t left = o.forward_start, right = o.forward_start + o.insert_bp;
                    if (left >= h.block_begin && right <= h.block_end) ++in_block[k];
                    else if (right <= h.block_begin || left >= h.block_end) ++outside;
                    else ++junction[k];
                    if (o.orientation == panvar::MateOrientation::FirstForward) ++forward;
                    insert_sum += static_cast<double>(o.insert_bp);
                }
                duplicate_origins += duplicates;
                rows << '\t' << mass << '\t' << (equal ? "true" : "false") << '\t' << origins.size() << '\t'
                     << states.size() << '\t' << duplicates << '\t' << in_block[k] << '\t' << junction[k] << '\t'
                     << outside << '\t' << forward << '\t'
                     << (origins.empty() ? 0.0 : insert_sum / static_cast<double>(origins.size()));
            }
            const double ma = panvar::log_add_exp(scored.mass_of(f, haps[0].unique), scored.mass_of(f, haps[1].unique));
            const double mb = panvar::log_add_exp(scored.mass_of(f, haps[2].unique), scored.mass_of(f, haps[3].unique));
            const double ca = panvar::log_fragment_contribution(ma, background, library);
            const double cb = panvar::log_fragment_contribution(mb, background, library);
            sum_a += ca;
            sum_b += cb;
            const std::size_t block_a = in_block[0] + in_block[1] + junction[0] + junction[1];
            const std::size_t block_b = in_block[2] + in_block[3] + junction[2] + junction[3];
            std::string kind;
            if (block_a == 0 && block_b == 0) kind = "outside_block_on_all";
            else if (block_b == 0) kind = "block_on_A_only";
            else if (block_a == 0) kind = "block_on_B_only";
            else if (junction[0] + junction[1] + junction[2] + junction[3] == 0) kind = "block_interior_on_both";
            else kind = "block_junction_on_both";
            ClassSum& sum = classes[kind];
            ++sum.fragments;
            sum.delta += ca - cb;
            sum.deltas.push_back(ca - cb);
            all_deltas.push_back(ca - cb);
            rows << '\t' << ca << '\t' << cb << '\t' << (ca - cb) << '\t' << kind << '\n';
        }
        atomic_write(out + ".fragment_decomposition.tsv", rows.str());
        std::ostringstream summary;
        summary << std::setprecision(12) << "key\tvalue\nlabel_A\t" << label_a << "\nlabel_B\t" << label_b << "\nblock\t"
                << block_id << "\nfragments\t" << fragments.size() << "\nfragment_term_A\t" << sum_a
                << "\nfragment_term_A_scorer\t" << scored.scores[pair_of_declared[ia]].fragment_term << "\nfragment_term_B\t"
                << sum_b << "\nfragment_term_B_scorer\t" << scored.scores[pair_of_declared[ib]].fragment_term
                << "\nfragment_delta\t" << (sum_a - sum_b) << "\nrebuilt_mass_mismatches\t" << rebuilt_mismatches
                << "\nexact_duplicate_origins\t" << duplicate_origins << "\nintensity_per_bp_copy\t"
                << library.fragment_start_intensity << "\noutlier_weight\t" << library.outlier_weight << '\n';
        for (const Hap& h : haps) {
            summary << h.name << "_length\t" << h.length << '\n' << h.name << "_block_bp\t" << (h.block_end - h.block_begin)
                    << '\n' << h.name << "_exposure\t" << scored.exposure[h.unique] << '\n' << h.name
                    << "_expected_signal_fragments\t" << panvar::expected_signal_fragments(scored.exposure[h.unique], library)
                    << '\n';
        }
        summary << "exposure_term_A\t" << scored.scores[pair_of_declared[ia]].exposure_term << "\nexposure_term_B\t"
                << scored.scores[pair_of_declared[ib]].exposure_term << '\n';
        std::vector<double> magnitude = all_deltas;
        std::sort(magnitude.begin(), magnitude.end(), [](double x, double y) { return std::fabs(x) > std::fabs(y); });
        for (const std::size_t top : {1, 10, 50, 200}) {
            double part = 0.0;
            for (std::size_t i = 0; i < std::min(top, magnitude.size()); ++i) part += magnitude[i];
            summary << "top" << top << "_fragments_delta\t" << part << '\n';
        }
        for (auto& [kind, sum] : classes) {
            std::sort(sum.deltas.begin(), sum.deltas.end());
            summary << "class_" << kind << "_fragments\t" << sum.fragments << "\nclass_" << kind << "_delta\t" << sum.delta
                    << "\nclass_" << kind << "_median_delta\t" << sum.deltas[sum.deltas.size() / 2] << '\n';
        }
        atomic_write(out + ".decomposition_summary.tsv", summary.str());
        std::cout << "score-set " << donor << ": decomposition " << label_a << " - " << label_b << " fragment delta "
                  << (sum_a - sum_b) << ", rebuilt-mass mismatches " << rebuilt_mismatches << ", duplicate origins "
                  << duplicate_origins << '\n';
    }
    if (args.get("truth-gfa", "").empty()) return 0;

    // Truth enters only now.
    const SortedGraph truth_graph = load_truth_in_model_frame(
        args.required("truth-gfa"), reference, loaded.sorted, out + ".truth.sorted.gfa");
    std::vector<panvar::PathProjection> truth;
    for (const panvar::PathRecord* path : donor_paths(truth_graph.graph, donor)) {
        const auto projection = panvar::project_path_onto_locus(truth_graph.graph, loaded.bubbles,
                                                                loaded.built.report, path->name);
        if (!projection.has_value()) throw std::runtime_error("truth does not span the block chain");
        truth.push_back(*projection);
    }
    const auto spelling = [&](std::size_t b, std::size_t allele) -> const std::string& {
        return model.blocks()[b].alleles[allele].sequence;
    };
    const auto block_distance = [&](std::size_t b, std::size_t a1, std::size_t a2) -> long {
        if (!truth[0].blocks[b].has_value() || !truth[1].blocks[b].has_value()) return -1;
        const std::string& t1 = *truth[0].blocks[b];
        const std::string& t2 = *truth[1].blocks[b];
        return std::min(edit_distance(t1, spelling(b, a1)) + edit_distance(t2, spelling(b, a2)),
                        edit_distance(t1, spelling(b, a2)) + edit_distance(t2, spelling(b, a1)));
    };
    std::vector<long> block_floor(blocks, 0);
    long floor_total = 0;
    for (std::size_t b = 0; b < blocks; ++b) {
        if (!truth[0].blocks[b].has_value() || !truth[1].blocks[b].has_value()) continue;
        bool available = false;
        block_floor[b] = closest_allele_distance(*truth[0].blocks[b], model.blocks()[b], available) +
                         closest_allele_distance(*truth[1].blocks[b], model.blocks()[b], available);
        floor_total += block_floor[b];
    }

    std::ofstream table(out + ".scores.tsv");
    table << std::setprecision(12)
          << "label\tsequence_pair\thap1_alleles\thap2_alleles\tfragment_term\texposure_term\tprior\t"
             "total\texpected_total\treproduced\tdistance\tfloor\texcess\n";
    std::ofstream per_block(out + ".blocks.tsv");
    per_block << "label\tblock\tblock_id\tkind\tallele1\tallele2\tblock_distance\tblock_floor\tblock_excess\n";
    const auto join = [](const std::vector<std::size_t>& v) {
        std::string text;
        for (std::size_t i = 0; i < v.size(); ++i) {
            if (i != 0) text += ',';
            text += std::to_string(v[i]);
        }
        return text;
    };
    for (std::size_t i = 0; i < declared.size(); ++i) {
        const std::size_t slot = pair_of_declared[i];
        const double total = scored.scores[slot].log_likelihood + prior_of_pair[slot];
        long distance = 0;
        for (std::size_t b = 0; b < blocks; ++b) {
            const long d = block_distance(b, declared[i].first[b], declared[i].second[b]);
            if (d < 0) continue;
            distance += d;
            per_block << declared[i].label << '\t' << b << '\t' << model.blocks()[b].block_id << '\t'
                      << (model.blocks()[b].kind == panvar::GenotypeBlockKind::Variant ? "SV" : "context")
                      << '\t' << declared[i].first[b] << '\t' << declared[i].second[b] << '\t' << d << '\t'
                      << block_floor[b] << '\t' << (d - block_floor[b]) << '\n';
        }
        table << declared[i].label << '\t' << slot << '\t' << join(declared[i].first) << '\t'
              << join(declared[i].second) << '\t' << scored.scores[slot].fragment_term << '\t'
              << scored.scores[slot].exposure_term << '\t' << prior_of_pair[slot] << '\t' << total << '\t'
              << (declared[i].has_expected ? std::to_string(declared[i].expected) : std::string("NA"))
              << '\t' << (declared[i].has_expected ? "true" : "NA") << '\t' << distance << '\t'
              << floor_total << '\t' << (distance - floor_total) << '\n';
    }
    std::cout << "score-set " << donor << ": " << declared.size() << " declared, " << scored_pairs.size()
              << " distinct diplotypes, " << interned.sequences.size() << " unique haplotypes, " << checked
              << " expectations reproduced, " << seconds << "s scoring\n";
    return 0;
}

// ------------------------------------------------------------------------------------------------
// form-census: how much placement work is spent re-scanning read sequence that already occurred
// ------------------------------------------------------------------------------------------------

// A measurement taken BEFORE deciding whether to deduplicate oriented read forms, and deliberately not
// the optimisation itself. The tempting number -- what fraction of forms are duplicates -- is the wrong
// one: what costs is piece hits, and a duplicated form is plausibly a repetitive one contributing far
// more hits than average. So this reports the form collapse, the seed-entry collapse, and the exact
// weighted prediction, and lets them disagree.
int form_census(const Args& args) {
    const std::string reference = args.required("reference");
    const std::string stage =
        (std::filesystem::temp_directory_path() /
         ("genotype_probe_census_" + std::to_string(::getpid()) + ".sorted.gfa")).string();
    const std::string curation_path = args.get("panel-curation", "");
    const LoadedModel loaded =
        load_model(args.required("gfa"), reference, stage,
                   std::stoul(args.get("min-variant-bp", "50")),
                   curation_path.empty() ? std::vector<std::string>{}
                                         : read_panel_curation(curation_path));
    std::filesystem::remove(stage);
    const panvar::PreparedGenotypeModel& model = *loaded.model;

    // The panel is what placement scans: one haplotype per panel row.
    std::vector<panvar::MosaicPath> paths;
    std::vector<std::string> sequences;
    for (std::size_t row = 0; row < model.panel().template_count(); ++row) {
        paths.push_back(panvar::MosaicPath{model.panel().template_alleles[row]});
        sequences.push_back(panvar::spell_mosaic(model, paths.back()).bases);
    }
    const panvar::InternedHaplotypes interned = panvar::intern_haplotypes(paths, sequences);

    const std::vector<panvar::ReadPair> fragments =
        read_pairs(args.required("reads1"), args.required("reads2"));
    panvar::FragmentScoringOptions options;
    options.mismatch_band = std::stoul(args.get("band", "6"));
    options.minimum_piece_bp = std::stoul(args.get("minimum-piece-bp", "12"));

    const panvar::PieceHitCensus census = panvar::census_piece_hits(interned, fragments, options);

    const auto saving = [](std::uint64_t kept, std::uint64_t whole) {
        return whole == 0 ? 0.0
                          : 100.0 * (1.0 - static_cast<double>(kept) / static_cast<double>(whole));
    };
    std::ostringstream out;
    out << std::fixed << std::setprecision(2);
    out << "key\tvalue\n"
        << "locus\t" << loaded.sorted.resolved << '\n'
        << "panel_curation\t" << (curation_path.empty() ? "none" : curation_path) << '\n'
        << "panel_rows\t" << model.panel().template_count() << '\n'
        << "unique_haplotypes\t" << interned.sequences.size() << '\n'
        << "fragments\t" << fragments.size() << '\n'
        << "total_forms\t" << census.total_forms << '\n'
        << "unique_forms\t" << census.unique_forms << '\n'
        << "form_collapse_percent\t" << saving(census.unique_forms, census.total_forms) << '\n'
        << "total_seed_entries\t" << census.total_entries << '\n'
        << "unique_seed_entries\t" << census.unique_entries << '\n'
        << "seed_entry_collapse_percent\t" << saving(census.unique_entries, census.total_entries) << '\n'
        << "piece_hits_all_forms\t" << census.predicted_piece_hits_all << '\n'
        << "piece_hits_unique_forms\t" << census.predicted_piece_hits_unique << '\n'
        // The number the deduplication decision rests on, and the only one of the three that does.
        << "predicted_piece_hit_reduction_percent\t"
        << saving(census.predicted_piece_hits_unique, census.predicted_piece_hits_all) << '\n'
        << "max_piece_hits_at_one_seed\t" << census.max_piece_hits_at_one_seed << '\n'
        << "max_piece_hits_one_haplotype\t" << census.max_piece_hits_one_haplotype << '\n';
    std::cout << out.str();
    if (!args.get("out", "").empty()) {
        std::ofstream file(args.get("out", ""));
        file << out.str();
    }
    return 0;
}

// ------------------------------------------------------------------------------------------------
// trace: follow the leave-zero-out truth pair through every proposal layer
// ------------------------------------------------------------------------------------------------

// Rebuilds the command's model, marker panel, depth fit and proposal emissions, then runs the proposal
// with the donor's own two panel rows declared as the trace target. The traced search is the production
// search (gated): its candidates are compared with the run's manifest.
int trace(const Args& args) {
    const std::string out = args.required("out");
    // Same staging leak compare had: the sorted graphs are an artefact of loading, 3.3 MB each, and
    // trace removed neither. A guard rather than a call at the end, because the target and fingerprint
    // checks below throw after the model graph is staged. (evaluate leaks these too, deliberately left
    // alone: its staged truth graph is the evidence a node-frame diagnosis reads.)
    struct StagedGraphs {
        std::string prefix;
        ~StagedGraphs() {
            std::error_code ignored;
            std::filesystem::remove(prefix + ".model.sorted.gfa", ignored);
            std::filesystem::remove(prefix + ".truth.sorted.gfa", ignored);
        }
    } staged_graphs{out};
    const std::string reference = args.required("reference");
    const std::string donor = args.required("donor");
    const std::string run = args.required("run");
    const std::size_t state_limit = std::stoul(args.get("state-limit", "1024"));
    const std::size_t target = std::stoul(args.get("candidates", "256"));
    const std::string curation_path = args.get("panel-curation", "");
    const LoadedModel loaded =
        load_model(args.required("gfa"), reference, out + ".model.sorted.gfa", 50,
                   curation_path.empty() ? std::vector<std::string>{}
                                         : read_panel_curation(curation_path));
    const panvar::PreparedGenotypeModel& model = *loaded.model;
    if (read_key_values(run + ".mosaic_status.tsv").at("model_fingerprint") != loaded.fingerprint) {
        throw std::runtime_error("the run was not scored on the rebuilt model");
    }
    const panvar::PanelAlleleMatrix& panel = model.panel();
    std::vector<std::size_t> rows;
    for (std::size_t row = 0; row < panel.template_count(); ++row) {
        if (panel.template_names[row].rfind(donor + "#", 0) == 0) rows.push_back(row);
    }
    std::sort(rows.begin(), rows.end(), [&](std::size_t lhs, std::size_t rhs) {
        return panel.template_names[lhs] < panel.template_names[rhs];
    });

    // Two ways to name the target. A leave-zero-out donor IS two panel rows, and its alleles are read
    // straight off them. A synthetic mosaic is deliberately in no panel row -- that is the whole point
    // of the experiment -- so its target comes from projecting the truth paths onto the block chain,
    // exactly as evaluate derives truth, and there are no declared templates to follow alongside it.
    panvar::MosaicDiplotype truth;
    // A declared target, given as "<hap1 allele list>|<hap2 allele list>". Needed to trace a candidate
    // that is not the donor's truth -- a repaired candidate, for instance, which is the winner with one
    // block replaced. It also sidesteps two dead ends for a leave-one-out donor: its panel rows are
    // gone, so rows.size() is 0 rather than 2, and the truth branch refuses outright when a projected
    // allele is not an allele of the model ("the target cannot be traced through a beam that could
    // never spell it"), which is exactly the case for a donor whose truth the panel cannot spell.
    const std::string target_spec = args.get("target-alleles", "");
    if (!target_spec.empty()) {
        const std::size_t bar = target_spec.find('|');
        if (bar == std::string::npos) {
            throw std::runtime_error("--target-alleles expects \"<hap1 alleles>|<hap2 alleles>\"");
        }
        panvar::MosaicPath first{parse_alleles(target_spec.substr(0, bar))};
        panvar::MosaicPath second{parse_alleles(target_spec.substr(bar + 1))};
        const std::size_t blocks = model.blocks().size();
        if (first.alleles.size() != blocks || second.alleles.size() != blocks) {
            throw std::runtime_error("--target-alleles needs one allele per block (" +
                                     std::to_string(blocks) + " per haplotype)");
        }
        for (std::size_t b = 0; b < blocks; ++b) {
            const std::size_t count = model.blocks()[b].alleles.size();
            if (first.alleles[b] >= count || second.alleles[b] >= count) {
                throw std::runtime_error("--target-alleles names an allele outside block " +
                                         model.blocks()[b].block_id);
            }
        }
        truth = panvar::MosaicDiplotype{std::move(first), std::move(second)};
    } else if (rows.size() == 2) {
        truth = panvar::MosaicDiplotype{
            panvar::MosaicPath{std::vector<std::size_t>(panel.template_alleles[rows[0]].begin(),
                                                        panel.template_alleles[rows[0]].end())},
            panvar::MosaicPath{std::vector<std::size_t>(panel.template_alleles[rows[1]].begin(),
                                                        panel.template_alleles[rows[1]].end())}};
    } else {
        const std::string truth_gfa = args.get("truth-gfa", "");
        if (truth_gfa.empty()) {
            throw std::runtime_error("donor " + donor + " has " + std::to_string(rows.size()) +
                                     " panel rows; pass --truth-gfa to trace a donor that is not in "
                                     "the panel, such as a curated-out synthetic mosaic");
        }
        const SortedGraph truth_graph =
            load_truth_in_model_frame(truth_gfa, reference, loaded.sorted, out + ".truth.sorted.gfa");
        const std::vector<const panvar::PathRecord*> paths = donor_paths(truth_graph.graph, donor);
        std::vector<panvar::MosaicPath> built;
        for (const panvar::PathRecord* path : paths) {
            const auto projection = panvar::project_path_onto_locus(
                truth_graph.graph, loaded.bubbles, loaded.built.report, path->name);
            if (!projection.has_value()) {
                throw std::runtime_error("truth path " + path->name + " does not span the block chain");
            }
            panvar::MosaicPath spelled;
            for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                if (!projection->blocks[b].has_value()) {
                    throw std::runtime_error("truth path " + path->name + " has no allele at block " +
                                             model.blocks()[b].block_id);
                }
                const auto& alleles = model.blocks()[b].alleles;
                std::size_t found = alleles.size();
                for (std::size_t a = 0; a < alleles.size(); ++a) {
                    if (alleles[a].sequence == *projection->blocks[b]) {
                        found = a;
                        break;
                    }
                }
                if (found == alleles.size()) {
                    throw std::runtime_error("truth allele at " + model.blocks()[b].block_id +
                                             " is not an allele of the model: the target cannot be "
                                             "traced through a beam that could never spell it");
                }
                spelled.alleles.push_back(found);
            }
            built.push_back(std::move(spelled));
        }
        truth = panvar::MosaicDiplotype{built[0], built[1]};
    }
    const auto winner = winner_paths(run);
    const panvar::MosaicDiplotype winner_diplotype{panvar::MosaicPath{winner.first}, panvar::MosaicPath{winner.second}};

    // The command's stages 2 and 3, with its defaults.
    const std::vector<panvar::ReadPair> fragments = read_pairs(args.required("reads1"), args.required("reads2"));
    const panvar::PreparedMarkerPanel markers = panvar::build_marker_panel(model, 31, 15);
    const panvar::MarkerObservations observations = panvar::count_marker_observations(markers, fragments);
    std::optional<double> declared_background;
    {
        std::ifstream library(args.required("library"));
        std::string key;
        double value = 0.0;
        for (std::string line; std::getline(library, line);) {
            std::istringstream fields(line);
            if (fields >> key >> value && key == "marker_background") declared_background = value;
        }
    }
    const panvar::MarkerDepthEstimate depth =
        panvar::estimate_marker_depth(markers, observations, declared_background);
    if (!depth.converged) throw std::runtime_error("the depth fit did not converge: " + depth.termination);
    const panvar::PreparedProposalEmissions emissions = panvar::build_proposal_emissions(
        markers, observations, panvar::proposal_marker_parameters(depth, 31, 15));
    const std::vector<double> switches(model.blocks().size() - 1, 0.01);
    panvar::MosaicBeamParameters parameters;
    // As the command does: return every diplotype found, so off-panel mosaics can be compared exactly.
    parameters.unique_candidate_target = std::max(state_limit, target + panel.template_count() * panel.template_count());
    parameters.state_limit = state_limit;
    parameters.max_extension_evaluations = 2000000000ULL;
    parameters.trace_target = truth;
    // Only a panel donor has declared templates to follow; a curated-out mosaic has none.
    if (rows.size() == 2) parameters.trace_templates = std::make_pair(rows[0], rows[1]);
    // --grouped selects the EXPERIMENTAL grouped-prefix proposer instead of the existing one. Both
    // return the same MosaicBeamResult, so every summary below is computed the same way; the grouped
    // one fills only the layer fields it can, so its per-layer trace is thinner.
    const bool grouped = args.get("grouped", "0") != "0";
    const panvar::MosaicBeamResult result =
        grouped ? panvar::propose_mosaic_diplotypes_grouped(panel, emissions, switches, parameters)
                : panvar::propose_mosaic_diplotypes(panel, emissions, switches, parameters);

    // The traced search must be the run's search.
    // A manifest without a source column predates the panel baseline: every row is a beam candidate, and the
    // beam's own order is compared. Otherwise only mosaic rows, against beam candidates that are not panel pairs.
    std::vector<std::map<std::string, std::string>> manifest;
    bool has_source = false;
    for (const auto& row : read_table(run + ".mosaic_candidates.tsv")) {
        has_source = has_source || row.count("source") != 0;
        if (row.count("source") == 0 || row.at("source") == "mosaic") manifest.push_back(row);
    }
    const panvar::MosaicBeamResult compared_beam =
        has_source ? panvar::off_panel_mosaics(result, panvar::panel_diplotypes(panel)) : result;
    std::vector<const panvar::MosaicCandidate*> beam_mosaics;
    for (const auto& candidate : compared_beam.candidates) beam_mosaics.push_back(&candidate);
    // Agreement between this standalone proposal and the production run's own manifest, on all three
    // of the things that identify a candidate: WHICH diplotype it is, WHERE it ranks (the comparison is
    // index by index, so agreeing in order is agreeing on rank), and WHAT it scored. Identity and rank
    // alone would pass an implementation that ordered the same candidates by a different quantity.
    //
    // The score is compared with a relative tolerance because the manifest is a written decimal, not
    // the double the search held; production and this path run the same code on the same inputs, so
    // anything beyond that tolerance is a real disagreement.
    bool matches_run = !manifest.empty();
    bool matches_run_score = !manifest.empty();
    const std::size_t compared = std::min(manifest.size(), beam_mosaics.size());
    for (std::size_t i = 0; i < compared; ++i) {
        const bool same_identity =
            parse_alleles(manifest[i].at("hap1_alleles")) == beam_mosaics[i]->diplotype.first.alleles &&
            parse_alleles(manifest[i].at("hap2_alleles")) == beam_mosaics[i]->diplotype.second.alleles;
        matches_run = matches_run && same_identity;
        bool same_score = same_identity;
        if (manifest[i].count("retained_logsum_proposal_score") != 0) {
            const double run_score = std::stod(manifest[i].at("retained_logsum_proposal_score"));
            const double here = beam_mosaics[i]->retained_logsum_proposal_score;
            same_score = same_score &&
                         std::abs(run_score - here) <= 1e-9 * std::max(1.0, std::abs(here));
        }
        matches_run_score = matches_run_score && same_score;
    }

    std::ofstream layers(out + ".trace.tsv");
    layers << std::setprecision(10)
           << "block\tblock_id\tkind\ttruth_alleles\ttruth_emission\twinner_alleles\twinner_emission\t"
              "best_emission_in_block\tstates_offered\tstates_retained\tstates_pruned\t"
              "prefix_pairs_offered\tprefix_pairs_retained\tmax_states_per_prefix\tretained_min_log_mass\t"
              "truth_prefixes_present\ttruth_states_offered\ttruth_states_retained\ttruth_best_log_mass\t"
              "truth_logsum_log_mass\ttruth_retained_logsum_log_mass\ttruth_lost_log_mass\t"
              "truth_group_rank_by_mass\ttruth_rank\tdonor_templates_offered\tdonor_templates_log_mass\t"
              "donor_templates_retained\ttruth_pair_local_rank\tlocal_pairs_ranked\t"
              "truth_allele1_local_rank\ttruth_allele2_local_rank\tlocal_alleles\t"
              "local_alleles_scored\n";
    double truth_emission_total = 0.0;
    double winner_emission_total = 0.0;
    std::string first_loss = "never";
    for (const panvar::MosaicBeamTraceLayer& layer : result.trace) {
        const std::size_t b = layer.block;
        const panvar::BlockProposalEmissions& emission = emissions.blocks[b];
        const double truth_e = emission.at(truth.first.alleles[b], truth.second.alleles[b]);
        const double winner_e = emission.at(winner.first[b], winner.second[b]);
        truth_emission_total += truth_e;
        winner_emission_total += winner_e;
        const double best_e = *std::max_element(emission.log_score.begin(), emission.log_score.end());
        // The TRUTH-FREE local standing of the declared target's alleles in this block. Nothing here
        // reads truth: it ranks the target's allele pair among every unordered pair the block offers,
        // by the same proposal emission the beam ranks on, and ranks each of its two alleles by the
        // best pairing that allele can reach. That is what a generic block-diversity rule would have to
        // see to reach these alleles without being told the answer.
        //
        // log_score is a symmetric allele_count x allele_count matrix, so best_e above is a PAIR score
        // and is commensurable with truth_e. Unscored alleles keep 0.0, which at these magnitudes ranks
        // them last rather than first; local_alleles_scored is reported so that stays auditable.
        const std::size_t local_alleles = emission.allele_count;
        std::size_t pair_rank = 1;
        std::size_t pairs_ranked = 0;
        std::vector<double> allele_best(local_alleles, -std::numeric_limits<double>::infinity());
        for (std::size_t x = 0; x < local_alleles; ++x) {
            for (std::size_t y = x; y < local_alleles; ++y) {
                const double e = emission.at(x, y);
                ++pairs_ranked;
                if (e > truth_e) ++pair_rank;
                if (e > allele_best[x]) allele_best[x] = e;
                if (e > allele_best[y]) allele_best[y] = e;
            }
        }
        const auto marginal_rank = [&](std::size_t allele) {
            std::size_t rank = 1;
            for (std::size_t x = 0; x < local_alleles; ++x) {
                if (allele_best[x] > allele_best[allele]) ++rank;
            }
            return rank;
        };
        if (first_loss == "never" && layer.target_states_retained == 0) first_loss = model.blocks()[b].block_id;
        layers << b << '\t' << model.blocks()[b].block_id << '\t'
               << (model.blocks()[b].kind == panvar::GenotypeBlockKind::Variant ? "SV" : "context") << '\t'
               << truth.first.alleles[b] << '/' << truth.second.alleles[b] << '\t' << truth_e << '\t'
               << winner.first[b] << '/' << winner.second[b] << '\t' << winner_e << '\t' << best_e << '\t'
               << layer.states_offered << '\t' << layer.states_retained << '\t' << layer.states_pruned << '\t'
               << layer.prefix_pairs_offered << '\t' << layer.prefix_pairs_retained << '\t'
               << layer.max_states_per_prefix << '\t'
               << layer.retained_min_log_mass << '\t'
               << (layer.target_prefixes_present ? "true" : "false") << '\t' << layer.target_states_offered << '\t'
               << layer.target_states_retained << '\t' << layer.target_best_log_mass << '\t'
               << layer.target_logsum_log_mass << '\t'
               << layer.target_retained_logsum_log_mass << '\t' << layer.target_lost_log_mass << '\t'
               << (layer.target_group_rank_by_offered_mass > 0
                       ? std::to_string(layer.target_group_rank_by_offered_mass)
                       : "NA")
               << '\t'
               << (layer.target_states_offered > 0 ? std::to_string(layer.offered_above_target_best + 1) : "NA")
               << '\t' << (layer.declared_templates_offered ? "true" : "false") << '\t'
               << layer.declared_templates_log_mass << '\t' << (layer.declared_templates_retained ? "true" : "false")
               << '\t' << pair_rank << '\t' << pairs_ranked << '\t'
               << marginal_rank(truth.first.alleles[b]) << '\t'
               << marginal_rank(truth.second.alleles[b]) << '\t' << local_alleles << '\t'
               << emission.scored_alleles
               << '\n';
    }
    bool truth_proposed = false;
    for (const auto& candidate : result.candidates) {
        const auto& d = candidate.diplotype;
        if ((d.first == truth.first && d.second == truth.second) || (d.first == truth.second && d.second == truth.first)) {
            truth_proposed = true;
        }
    }
    // The truth's rank in the FULL final beam, which is the question a ranked prefix cannot answer:
    // candidate_prefix() truncates to K before the run's manifest is written, so "absent from the
    // manifest" conflates ranked-beyond-K, pruned-at-a-layer, and never-generated. This scans
    // result.candidates, the complete retained set, in the beam's own order.
    std::string truth_full_beam_rank = "absent";
    for (const auto& candidate : result.candidates) {
        const auto& d = candidate.diplotype;
        if ((d.first == truth.first && d.second == truth.second) ||
            (d.first == truth.second && d.second == truth.first)) {
            truth_full_beam_rank = std::to_string(candidate.unique_rank);
            break;
        }
    }

    std::ofstream summary(out + ".trace.summary.tsv");
    summary << std::setprecision(12) << "key\tvalue\nproposer\t" << (grouped ? "grouped" : "state_only")
            << "\ndonor\t" << donor << "\nstate_limit\t" << state_limit
            << "\ncandidate_target\t" << target << "\nunique_candidates_found\t" << result.unique_candidates_found
            << "\nstates_pruned\t" << result.states_pruned
            << "\nextensions_evaluated\t" << result.extensions_evaluated
            << "\nproposal_payload_predicted_bytes\t" << result.allocation_predicted.total_bytes
            << "\nproposal_payload_observed_bytes\t" << result.allocation_observed.total_bytes
            << "\nproposal_payload_predicted_rolling_bytes\t" << result.allocation_predicted.rolling_vector_bytes
            << "\nproposal_payload_predicted_ancestry_bytes\t" << result.allocation_predicted.ancestry_bytes
            << "\nproposal_payload_predicted_pair_bytes\t" << result.allocation_predicted.pair_bytes
            << "\nproposal_payload_predicted_index_bytes\t" << result.allocation_predicted.index_bytes
            << "\nproposal_payload_predicted_child_records\t" << result.allocation_predicted.child_records
            << "\nproposal_payload_observed_child_records\t" << result.allocation_observed.child_records
            << "\nproposal_matches_run\t"
            << (matches_run ? "true" : "false") << "\nproposal_matches_run_score\t"
            << (matches_run_score ? "true" : "false") << "\ncandidates_compared_with_run\t" << compared
            // Empty for a synthetic mosaic: it is curated out of the panel by construction, so it has no
            // panel rows. Indexing rows[0] here regardless is what segfaulted the first mosaic trace --
            // an out-of-bounds read on an empty vector, silent, with no exception to catch.
            << "\ntruth_templates\t"
            << (rows.size() == 2 ? panel.template_names[rows[0]] + "," + panel.template_names[rows[1]]
                                 : std::string("none: donor is not a panel row"))
            << "\ntruth_in_full_beam\t" << (truth_proposed ? "true" : "false")
            << "\ntruth_full_beam_rank\t" << truth_full_beam_rank
            << "\nfirst_block_truth_not_retained\t"
            << first_loss << "\ntruth_emission_total\t" << truth_emission_total << "\nwinner_emission_total\t"
            << winner_emission_total << "\ntruth_path_prior\t"
            << panvar::log_unordered_diplotype_prior(panel, truth, switches) << "\nwinner_path_prior\t"
            << panvar::log_unordered_diplotype_prior(panel, winner_diplotype, switches) << "\nmarker_depth_rate\t"
            << depth.rate_per_copy << "\nmarker_depth_termination\t" << depth.termination << '\n';
    std::cout << "traced " << donor << ": truth proposed " << (truth_proposed ? "yes" : "no") << ", first lost at "
              << first_loss << ", proposal matches run " << (matches_run ? "yes" : "no") << '\n';
    return 0;
}

// ------------------------------------------------------------------------------------------------
// genotype: the command with the development candidate-prefix harness, which the CLI does not expose
// ------------------------------------------------------------------------------------------------

int genotype(const Args& args) {
    panvar::GenotypeCommandOptions options;
    options.gfa_path = args.required("gfa");
    options.reference_path = args.required("reference");
    options.reads1_path = args.required("reads1");
    options.reads2_path = args.required("reads2");
    options.library_path = args.required("library");
    options.out_prefix = args.required("out-prefix");
    options.proposal_state_limit = std::stoul(args.get("state-limit", "1024"));
    options.phase_audit_candidates = std::stoul(args.get("phase-audit", "0"));
    // Workers for the diploid combination phase. Absent until now, so the probe accepted --threads as an
    // unknown key and dropped it: a benchmark invoked with --threads 2 recorded threads_requested 1 and
    // measured the serial path four times over.
    options.threads = std::stoul(args.get("threads", "1"));
    // Optional curation manifest. The curated panel is a different model and says so: the status
    // records which templates left, why, and a digest of the decision, and the fingerprint changes.
    options.curation_manifest_path = args.get("panel-curation", "");
    // Explicit scorer-cap overrides. The defaults are unchanged and stay the defaults; these exist
    // so a declared performance experiment can raise them without moving them for everyone, and
    // every value lands in the status so the run signature records what it actually ran under.
    options.max_state_classes = std::stoull(args.get("max-state-classes", "200000000"));
    options.max_physical_origins = std::stoull(args.get("max-physical-origins", "500000000"));
    options.emission_cache_bytes = std::stoull(args.get("emission-cache-bytes", "0"));
    for (const std::size_t k : parse_alleles(args.required("prefixes"))) {
        options.development_candidate_prefixes.push_back(k);
    }
    if (options.development_candidate_prefixes.empty()) throw std::runtime_error("--prefixes is empty");
    options.max_mosaic_candidates = options.development_candidate_prefixes.back();
    options.quiet = true;
    // Opt-in census of the verification stream. It answers one question -- would exact reuse of
    // identical emission calculations pay -- and is off unless a destination is named.
    const std::string census_path = args.get("emission-census", "");
    if (!census_path.empty()) panvar::begin_emission_census();
    const panvar::GenotypeRunStatus status = panvar::run_genotype(options);
    if (!census_path.empty()) {
        const panvar::EmissionCensus census = panvar::end_emission_census();
        const auto ratio = [](std::uint64_t requests, std::uint64_t distinct) {
            return distinct == 0 ? 0.0 : static_cast<double>(requests) / static_cast<double>(distinct);
        };
        std::ofstream out(census_path);
        out << std::setprecision(6) << "key\tvalue\n"
            << "emission_requests\t" << census.requests << '\n'
            << "distinct_emission_keys\t" << census.distinct << '\n'
            << "reuse_factor_all\t" << ratio(census.requests, census.distinct) << '\n'
            << "prefilter_refusal_requests\t" << census.prefilter_refusal_requests << '\n'
            << "prefilter_refusal_distinct\t" << census.prefilter_refusal_distinct << '\n'
            << "reuse_factor_prefilter_refusals\t"
            << ratio(census.prefilter_refusal_requests, census.prefilter_refusal_distinct) << '\n'
            << "dp_requests\t" << census.dp_requests << '\n'
            << "dp_distinct\t" << census.dp_distinct << '\n'
            // The decision number: the factor by which exact reuse would reduce executed DP calls.
            << "reuse_factor_dp\t" << ratio(census.dp_requests, census.dp_distinct) << '\n'
            << "predicted_key_bytes\t" << census.key_bytes << '\n'
            << "predicted_value_bytes\t" << census.value_bytes << '\n'
            << "predicted_cache_payload_bytes\t" << census.key_bytes + census.value_bytes << '\n';
        std::cout << "emission census: " << census.requests << " requests, " << census.distinct
                  << " distinct (reuse " << ratio(census.requests, census.distinct)
                  << "x; DP reuse " << ratio(census.dp_requests, census.dp_distinct) << "x)\n";
    }
    std::cout << "genotype " << options.out_prefix << " state limit " << options.proposal_state_limit
              << ": largest prefix " << panvar::genotype_run_status_name(status) << '\n';
    return 0;
}

// ------------------------------------------------------------------------------------------------
// tie-branch: DIAGNOSTIC ONLY (v2/post-hoc). The production call under a profile, with the improving-move-tie
// branching diagnostic switched on. Production outputs are written exactly as the plain call writes them;
// the diagnostic adds <out-prefix>.diag.* files. Truth is never an input here: evaluate separately.
int tie_branch(const Args& args) {
    panvar::GenotypeCommandOptions options;
    panvar::apply_genotype_profile(args.required("profile"), options);
    options.threads = std::stoul(args.required("threads"));
    options.overridden_settings.push_back("threads=" + std::to_string(options.threads));
    options.gfa_path = args.required("gfa");
    options.reference_path = args.required("reference");
    options.reads1_path = args.required("reads1");
    options.reads2_path = args.required("reads2");
    options.library_path = args.required("library");
    options.curation_manifest_path = args.get("panel-curation", "");
    options.out_prefix = args.required("out-prefix");
    options.quiet = true;
    options.diagnostic_tie_branching = args.get("diagnostic", "1") != "0";
    options.diagnostic_max_branches = std::stoul(args.get("max-branches", "16"));
    options.diagnostic_time_limit_seconds = std::stod(args.get("diagnostic-time-limit-seconds", "3600"));
    options.diagnostic_footprint_limit_bytes = std::stoull(args.get("diagnostic-footprint-limit-bytes", "3758096384"));
    options.diagnostic_max_total_physical_origins =
        std::stoull(args.get("diagnostic-max-total-physical-origins", "19600000000"));
    options.diagnostic_reverse_branch_order = args.get("reverse-branches", "0") != "0";
    options.diagnostic_reuse_trace = args.get("reuse-trace", "0") != "0";
    options.diagnostic_warm_branches = args.get("warm-branches", "0") != "0";
    options.diagnostic_digest_memo = args.get("digest-memo", "0") != "0";
    options.diagnostic_digest_memo_check = args.get("digest-memo-check", "0") != "0";
    options.diagnostic_branch_request_order_seed = std::stoull(args.get("branch-request-order-seed", "0"));
    options.diagnostic_pool_dump = args.get("pool-dump", "0") != "0";
    // v2 production gates: overrides are recorded; the reversal is a test hook, never a setting.
    options.test_reverse_tie_branch_order = args.get("tie-reverse", "0") != "0";
    for (const char* key : {"repair-max-physical-origins", "repair-tie-max-branches"}) {
        const std::string value = args.get(key, "");
        if (value.empty()) continue;
        if (std::string(key) == "repair-max-physical-origins") options.repair_max_physical_origins = std::stoull(value);
        else options.repair_tie_max_branches = std::stoul(value);
        std::string name = key;
        std::replace(name.begin(), name.end(), '-', '_');
        options.overridden_settings.push_back(name + "=" + value);
    }
    // EXPERIMENT ONLY: a higher union physical-origin cap for one diagnosed refusal. Recorded as an override;
    // no profile or default moves.
    if (!args.get("union-max-physical-origins", "").empty()) {
        options.max_physical_origins = std::stoull(args.get("union-max-physical-origins", ""));
        options.overridden_settings.push_back("max_physical_origins=" + std::to_string(options.max_physical_origins));
    }
    const panvar::GenotypeRunStatus status = panvar::run_genotype(options);
    const auto diag = read_key_values(options.out_prefix + ".diag.status.tsv");
    std::cout << "tie-branch " << options.out_prefix << ": production " << panvar::genotype_run_status_name(status)
              << ", diagnostic " << (diag.count("state") ? diag.at("state") : std::string("not run")) << ' '
              << (diag.count("decision") ? diag.at("decision") : std::string("")) << '\n';
    return 0;
}

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            std::cerr << "usage: genotype_probe simulate|loo-gfa|evaluate|compare|genotype|trace|anchors|subsample|preflight|respell|coverage-audit|form-census|mosaic-truth|orientations|panel-alleles|tied-consensus|score-set|repair-search --key value ...\n";
            return 2;
        }
        const std::string command = argv[1];
        const Args args = parse_args(argc, argv, 2);
        if (command == "simulate") return simulate(args);
        if (command == "loo-gfa") return loo_gfa(args);
        if (command == "evaluate") return evaluate(args);
        if (command == "compare") return compare(args);
        if (command == "genotype") return genotype(args);
        if (command == "trace") return trace(args);
        if (command == "anchors") return anchors(args);
        if (command == "subsample") return subsample(args);
        if (command == "preflight") return preflight(args);
        if (command == "respell") return respell(args);
        if (command == "coverage-audit") return coverage_audit(args);
        if (command == "form-census") return form_census(args);
        if (command == "mosaic-truth") return mosaic_truth(args);
        if (command == "orientations") return orientations(args);
        if (command == "panel-alleles") return panel_alleles(args);
        if (command == "tied-consensus") return tied_consensus(args);
        if (command == "repair-pilot") return repair_pilot(args);
        if (command == "score-set") return score_set(args);
        if (command == "repair-search") return repair_search(args);
        if (command == "window-bench") return window_bench(args);
        if (command == "tie-branch") return tie_branch(args);
        std::cerr << "unknown probe command: " << command << '\n';
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "genotype_probe: " << error.what() << '\n';
        return 1;
    }
}
