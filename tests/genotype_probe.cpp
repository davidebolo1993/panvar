// Developer probe for the mosaic genotyper: a seeded paired-read simulator, a leave-one-out graph
// filter, and a truth evaluator. It is not part of the product surface and is never installed.
//
// Everything it measures is computed against a truth the product never sees. The evaluator rebuilds
// the model from the same GFA the run used and refuses unless that model reproduces the run's model
// fingerprint; truth paths are then projected from the full graph onto that model's block frame.

#include "edlib.h"
#include "panvar/bubbles.hpp"
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_command.hpp"
#include "panvar/genotype_locus.hpp"
#include "panvar/genotype_model.hpp"
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
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <map>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
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
};

SortedGraph load_sorted(const std::string& gfa_path, const std::string& reference,
                        const std::string& stage_path) {
    SortedGraph out;
    out.gfa = panvar::read_gfa_model(gfa_path);
    panvar::GraphSortOptions sort_options;
    sort_options.reference_path = reference;
    sort_options.flip = true;
    out.resolved = panvar::sort_graph_reference(out.gfa, sort_options).resolved_reference;
    panvar::write_gfa_model(stage_path, out.gfa);
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
                             bool& exactly_available) {
    std::vector<const std::string*> alleles;
    for (const auto& allele : block.alleles) alleles.push_back(&allele.sequence);
    std::sort(alleles.begin(), alleles.end(), [&](const std::string* lhs, const std::string* rhs) {
        const long dl = std::labs(static_cast<long>(lhs->size()) - static_cast<long>(truth.size()));
        const long dr = std::labs(static_cast<long>(rhs->size()) - static_cast<long>(truth.size()));
        return dl < dr;
    });
    long best = -1;
    exactly_available = false;
    for (const std::string* allele : alleles) {
        const long length_gap =
            std::labs(static_cast<long>(allele->size()) - static_cast<long>(truth.size()));
        if (best >= 0 && length_gap >= best) break;   // length difference bounds edit distance below
        const long distance = edit_distance(truth, *allele, best < 0 ? -1 : best - 1);
        if (distance >= 0 && (best < 0 || distance < best)) best = distance;
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

    // Truth comes from the full graph. Sorting depends only on the reference path, so a node id names
    // the same sequence in both graphs; that is checked here rather than assumed.
    const SortedGraph truth_graph =
        load_sorted(args.required("truth-gfa"), reference, out + ".truth.sorted.gfa");
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
             "truth2_available\tgenotype_support\n";

    const bool complete = status.count("status") && status.at("status") == "COMPLETE";
    std::vector<std::size_t> hap1;
    std::vector<std::size_t> hap2;
    std::map<std::string, std::string> genotype_support;
    const std::string phase_support = status.count("phase_support") ? status.at("phase_support") : "NA";
    // Phase-aware accuracy only where the reads chose the phase: within the candidates, or audited.
    const bool phase_claimed = phase_support == "READS_WITHIN_CANDIDATES" || phase_support == "READS_AUDITED";
    if (complete) {
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
                  << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
            continue;
        }
        bool available1 = false;
        bool available2 = false;
        const long floor1 = closest_allele_distance(*t1, block, available1);
        const long floor2 = closest_allele_distance(*t2, block, available2);
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
              << (genotype_support.count(block.block_id) ? genotype_support.at(block.block_id) : "NA") << '\n';
    }

    // The only no-call: no genotype result exists.
    summary << "run_no_call\t" << (complete ? "false" : "true") << '\n';
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
    } else {
        summary << "accuracy\texcluded: run status is not COMPLETE\n";
    }
    std::cout << "evaluated " << scheme << " " << donor << " K=" << status.at("max_mosaic_candidates") << '\n';
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
std::vector<panvar::FragmentPlacement> scan_origins(const std::string& haplotype, const panvar::ReadPair& fragment,
                                                    const panvar::InsertSupport& support, std::size_t band) {
    const auto starts = [&](const std::string& read) {
        std::vector<std::pair<std::size_t, std::size_t>> out;
        if (haplotype.size() < read.size()) return out;
        for (std::size_t p = 0; p + read.size() <= haplotype.size(); ++p) {
            const auto mismatches = panvar::count_mismatches(haplotype, read, p, band);
            if (mismatches.has_value()) out.emplace_back(p, *mismatches);
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
        const auto r_starts = starts(reverse);
        for (const auto& fs : f_starts) {
            for (const auto& rs : r_starts) {
                const auto insert = panvar::fr_insert_length(fs.first, forward.size(), rs.first, reverse.size(), support);
                if (!insert.has_value()) continue;
                panvar::FragmentPlacement placement;
                placement.forward_start = fs.first;
                placement.reverse_start = rs.first;
                placement.insert_bp = *insert;
                placement.forward_mismatches = fs.second;
                placement.reverse_mismatches = rs.second;
                placement.orientation = arm == 0 ? panvar::MateOrientation::FirstForward : panvar::MateOrientation::SecondForward;
                out.push_back(placement);
            }
        }
    }
    return out;
}

auto placement_key(const panvar::FragmentPlacement& p) {
    return std::make_tuple(static_cast<int>(p.orientation), p.forward_start, p.reverse_start, p.insert_bp,
                           p.forward_mismatches, p.reverse_mismatches);
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
    const LoadedModel loaded = load_model(args.required("model-gfa"), reference, out + ".model.sorted.gfa", 50);
    const panvar::PreparedGenotypeModel& model = *loaded.model;
    for (const char* run : {"run-a", "run-b"}) {
        if (read_key_values(args.required(run) + ".mosaic_status.tsv").at("model_fingerprint") != loaded.fingerprint) {
            throw std::runtime_error(std::string(run) + " was not scored on the rebuilt model");
        }
    }
    const auto a = winner_paths(args.required("run-a"));
    const auto b = winner_paths(args.required("run-b"));
    std::size_t block = model.blocks().size();
    for (std::size_t i = 0; i < model.blocks().size(); ++i) {
        if (model.blocks()[i].block_id == block_id) block = i;
    }
    if (block == model.blocks().size()) throw std::runtime_error("unknown block " + block_id);

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
    const SortedGraph truth_graph = load_sorted(args.required("truth-gfa"), reference, out + ".truth.sorted.gfa");
    std::vector<std::string> truth_blocks;
    for (const panvar::PathRecord* path : donor_paths(truth_graph.graph, args.required("donor"))) {
        const auto projection = panvar::project_path_onto_locus(truth_graph.graph, loaded.bubbles, loaded.built.report, path->name);
        if (!projection.has_value() || !projection->blocks[block].has_value()) throw std::runtime_error("truth has no frame at the block");
        truth_blocks.push_back(*projection->blocks[block]);
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
            auto scanned = scan_origins(all.sequences[h], subset[i], support, band);
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
            const double scanned_mass = panvar::log_placement_mass(scanned, subset[i], insert, library, band);
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
        const double calibration_gap = scored.scores[1].fragment_term - scored.scores[0].fragment_term;
        const bool expected = item.first == "from_a" ? calibration_gap < 0.0 : calibration_gap > 0.0;
        summary << "calibration_" << item.first << "_fragments\t" << reads.size() << "\ncalibration_" << item.first
                << "_gap_b_minus_a\t" << calibration_gap << "\ncalibration_" << item.first << "_prefers_generator\t"
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
    LoadedModel loaded;
    std::string model_refusal;
    try {
        loaded = load_model(args.required("gfa"), reference, stage, min_variant_bp);
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

    // Switch plan: which donor supplies each block. "every" alternates at every block; otherwise a
    // comma-separated list of block indices at which the source flips.
    const std::string spec = args.get("switches", "");
    std::vector<int> donor_of_block(blocks, 0);
    if (spec == "every") {
        for (std::size_t b = 0; b < blocks; ++b) donor_of_block[b] = static_cast<int>(b % 2);
    } else {
        std::set<std::size_t> flips;
        std::stringstream ss(spec);
        std::string token;
        while (std::getline(ss, token, ',')) {
            if (!token.empty()) flips.insert(std::stoul(token));
        }
        if (flips.empty()) throw std::runtime_error("--switches must name at least one block, or 'every'");
        int current = 0;
        for (std::size_t b = 0; b < blocks; ++b) {
            if (flips.count(b) != 0) current = 1 - current;
            donor_of_block[b] = current;
        }
    }
    const std::size_t switches = [&] {
        std::size_t n = 0;
        for (std::size_t b = 1; b < blocks; ++b) n += donor_of_block[b] != donor_of_block[b - 1];
        return n;
    }();
    if (switches == 0) throw std::runtime_error("the plan produces no switch: that is not a mosaic");

    // Route one: the allele path, and the truth it defines.
    std::vector<panvar::MosaicPath> mosaic_paths(2);
    for (std::size_t h = 0; h < 2; ++h) {
        mosaic_paths[h].alleles.resize(blocks);
        for (std::size_t b = 0; b < blocks; ++b) {
            const std::size_t row = donor_of_block[b] == 0 ? rows_a[h] : rows_b[h];
            mosaic_paths[h].alleles[b] = model.panel().template_alleles[row][b];
        }
    }
    // Neither mosaic may coincide with a panel row: if it does, the correct answer is reachable without
    // recombining and the experiment tests nothing it claims to.
    for (std::size_t h = 0; h < 2; ++h) {
        for (std::size_t r = 0; r < model.panel().template_count(); ++r) {
            if (model.panel().template_alleles[r] == mosaic_paths[h].alleles) {
                throw std::runtime_error("mosaic haplotype " + std::to_string(h + 1) +
                                         " is exactly panel row " + model.panel().template_names[r] +
                                         "; choose different donors or switch points");
            }
        }
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
            const std::size_t row = donor_of_block[b] == 0 ? rows_a[h] : rows_b[h];
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
        for (std::size_t h = 0; h < 2; ++h) truth << '\t' << (donor_of_block[b] == 0 ? donor_a : donor_b);
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

    std::cout << "mosaic-truth " << donor_a << " x " << donor_b << ": " << switches
              << " switch(es) over " << blocks << " blocks; both routes agree; staged " << out_gfa
              << '\n';
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
    if (rows.size() == 2) {
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
        const SortedGraph truth_graph = load_sorted(truth_gfa, reference, out + ".truth.sorted.gfa");
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
    const panvar::MosaicBeamResult result =
        panvar::propose_mosaic_diplotypes(panel, emissions, switches, parameters);

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
    bool matches_run = !manifest.empty();
    const std::size_t compared = std::min(manifest.size(), beam_mosaics.size());
    for (std::size_t i = 0; i < compared && matches_run; ++i) {
        matches_run = parse_alleles(manifest[i].at("hap1_alleles")) == beam_mosaics[i]->diplotype.first.alleles &&
                      parse_alleles(manifest[i].at("hap2_alleles")) == beam_mosaics[i]->diplotype.second.alleles;
    }

    std::ofstream layers(out + ".trace.tsv");
    layers << std::setprecision(10)
           << "block\tblock_id\tkind\ttruth_alleles\ttruth_emission\twinner_alleles\twinner_emission\t"
              "best_emission_in_block\tstates_offered\tstates_retained\tstates_pruned\t"
              "prefix_pairs_offered\tprefix_pairs_retained\tmax_states_per_prefix\tretained_min_log_mass\t"
              "truth_prefixes_present\ttruth_states_offered\ttruth_states_retained\ttruth_best_log_mass\t"
              "truth_logsum_log_mass\ttruth_rank\tdonor_templates_offered\tdonor_templates_log_mass\t"
              "donor_templates_retained\n";
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
               << (layer.target_states_offered > 0 ? std::to_string(layer.offered_above_target_best + 1) : "NA")
               << '\t' << (layer.declared_templates_offered ? "true" : "false") << '\t'
               << layer.declared_templates_log_mass << '\t' << (layer.declared_templates_retained ? "true" : "false")
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
    summary << std::setprecision(12) << "key\tvalue\ndonor\t" << donor << "\nstate_limit\t" << state_limit
            << "\ncandidate_target\t" << target << "\nunique_candidates_found\t" << result.unique_candidates_found
            << "\nstates_pruned\t" << result.states_pruned
            << "\nextensions_evaluated\t" << result.extensions_evaluated << "\nproposal_matches_run\t"
            << (matches_run ? "true" : "false") << "\ncandidates_compared_with_run\t" << compared
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
    for (const std::size_t k : parse_alleles(args.required("prefixes"))) {
        options.development_candidate_prefixes.push_back(k);
    }
    if (options.development_candidate_prefixes.empty()) throw std::runtime_error("--prefixes is empty");
    options.max_mosaic_candidates = options.development_candidate_prefixes.back();
    options.quiet = true;
    const panvar::GenotypeRunStatus status = panvar::run_genotype(options);
    std::cout << "genotype " << options.out_prefix << " state limit " << options.proposal_state_limit
              << ": largest prefix " << panvar::genotype_run_status_name(status) << '\n';
    return 0;
}

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            std::cerr << "usage: genotype_probe simulate|loo-gfa|evaluate|compare|genotype|trace|anchors|subsample|preflight --key value ...\n";
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
        std::cerr << "unknown probe command: " << command << '\n';
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "genotype_probe: " << error.what() << '\n';
        return 1;
    }
}
