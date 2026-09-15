// Developer probe for the mosaic genotyper: a seeded paired-read simulator, a leave-one-out graph
// filter, and a truth evaluator. It is not part of the product surface and is never installed.
//
// Everything it measures is computed against a truth the product never sees. The evaluator rebuilds
// the model from the same GFA the run used and refuses unless that model reproduces the run's model
// fingerprint; truth paths are then projected from the full graph onto that model's block frame.

#include "edlib.h"
#include "panvar/bubbles.hpp"
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_locus.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/gfa.hpp"
#include "panvar/gfa_io.hpp"
#include "panvar/graph_sort.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/gz_reader.hpp"
#include "panvar/integrated_snarls.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
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

LoadedModel load_model(const std::string& gfa_path, const std::string& reference,
                       const std::string& stage_path, std::size_t min_variant_bp) {
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
    const std::uint64_t seed = std::stoull(args.get("seed", "7"));

    panvar::LibraryModel library;
    library.insert_mean_bp = std::stod(args.get("insert-mean", "350"));
    library.insert_sd_bp = std::stod(args.get("insert-sd", "40"));
    library.support_sigmas = std::stod(args.get("sigmas", "4"));
    library.substitution_error = std::stod(args.get("error", "0.002"));
    library.outlier_weight = std::stod(args.get("outlier", "0.001"));
    // Templates per base per chromosome copy, summed over both orientations: per-copy coverage is
    // depth/2, and each template contributes two mates of read_length bases.
    library.fragment_start_intensity = (depth / 2.0) / (2.0 * static_cast<double>(read_length));
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
    manifest << "key\tvalue\n" << "donor\t" << donor << "\nseed\t" << seed << "\ndepth\t" << depth
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
        std::poisson_distribution<long> count_draw(library.fragment_start_intensity * exposure);
        std::discrete_distribution<std::size_t> length_draw(weights.begin(), weights.end());
        std::bernoulli_distribution orientation_draw(0.5);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        const long templates = count_draw(rng);
        manifest << "haplotype_" << p << "\t" << paths[p]->name << "\nhaplotype_" << p << "_bp\t"
                 << haplotype.size() << "\nhaplotype_" << p << "_expected_templates\t"
                 << library.fragment_start_intensity * exposure << "\nhaplotype_" << p
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
    lib << "# declared simulation parameters, genotype_probe simulate seed " << seed << '\n'
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
    for (const char* key : {"status", "max_candidates", "mismatch_band", "unique_candidates",
                            "sequence_candidates", "winner_decided_by", "certified_gap",
                            "fragment_term_gap_to_runner_up", "winner_candidate_set_probability",
                            "seconds_total", "seconds_rescoring", "model_fingerprint",
                            "marker_depth_rate_per_copy", "marker_background"}) {
        const auto found = status.find(key);
        summary << key << '\t' << (found == status.end() ? "NA" : found->second) << '\n';
    }

    const LoadedModel loaded = load_model(args.required("model-gfa"), reference,
                                          out + ".model.sorted.gfa", 50);
    if (status.count("model_fingerprint") == 0 || status.at("model_fingerprint") != loaded.fingerprint) {
        throw std::runtime_error("the rebuilt model does not reproduce the run's model fingerprint");
    }
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
    table << "scheme\tdonor\tmax_candidates\tblock_id\tkind\ttruth1_bp\ttruth2_bp\tcalled1_bp\t"
             "called2_bp\tunordered_distance\tavailability_floor\texcess\ttruth1_available\t"
             "truth2_available\n";

    const bool complete = status.count("status") && status.at("status") == "COMPLETE";
    std::vector<std::size_t> hap1;
    std::vector<std::size_t> hap2;
    if (complete) {
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
            table << scheme << '\t' << donor << '\t' << status.at("max_candidates") << '\t'
                  << block.block_id << '\t' << (variant ? "SV" : "context")
                  << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
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
        table << scheme << '\t' << donor << '\t' << status.at("max_candidates") << '\t' << block.block_id
              << '\t' << (variant ? "SV" : "context") << '\t' << t1->size() << '\t' << t2->size() << '\t'
              << called1_bp << '\t' << called2_bp << '\t' << distance_cell << '\t' << floor << '\t'
              << excess_cell << '\t' << (available1 ? "true" : "false") << '\t'
              << (available2 ? "true" : "false") << '\n';
    }

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
                << "\nphased_distance_total\t" << std::min(order_a, order_b)
                << "\nphase_switch_blocks\t" << switch_blocks
                << "\nblocks_exact\t" << all_exact << "\nsv_blocks_exact\t" << sv_exact << '\n';

        // Truth in candidate set: some retained beam candidate spells both truth haplotypes.
        std::string truth_proposal_rank = "absent";
        std::string truth_sequence_candidate = "absent";
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
        summary << "truth_in_candidates\t" << (truth_proposal_rank != "absent" ? "true" : "false")
                << "\ntruth_proposal_rank\t" << truth_proposal_rank
                << "\ntruth_fragment_score_rank\t" << truth_score_rank << '\n';
    } else {
        summary << "accuracy\texcluded: run status is not COMPLETE\n";
    }
    std::cout << "evaluated " << scheme << " " << donor << " K=" << status.at("max_candidates") << '\n';
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            std::cerr << "usage: genotype_probe simulate|loo-gfa|evaluate --key value ...\n";
            return 2;
        }
        const std::string command = argv[1];
        const Args args = parse_args(argc, argv, 2);
        if (command == "simulate") return simulate(args);
        if (command == "loo-gfa") return loo_gfa(args);
        if (command == "evaluate") return evaluate(args);
        std::cerr << "unknown probe command: " << command << '\n';
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "genotype_probe: " << error.what() << '\n';
        return 1;
    }
}
