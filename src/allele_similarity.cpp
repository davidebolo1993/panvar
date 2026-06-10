#include "allele_internal.hpp"

#include "call_internal.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panvar {
namespace {

void write_similarity_alleles_tsv(
    const std::string& output_path,
    const Bubble& bubble,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<std::size_t>& cluster_of_unique) {

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write similarity allele table: " + output_path);
    }

    out << "bubble_id\tallele_id\tcluster_id\tpath_support\tsequence_length\tmode\tsignature\n";
    for (std::size_t i = 0; i < unique_alleles.size(); ++i) {
        const auto& allele = unique_alleles[i];
        out << bubble.id << '\t'
            << allele.allele_id << '\t'
            << cluster_of_unique[i] << '\t'
            << allele.path_support << '\t'
            << allele.sequence_length << '\t'
            << (allele.uses_sequence_similarity ? "sequence" : "walk") << '\t'
            << allele.signature << '\n';
    }
}

void write_similarity_distance_matrix_tsv(
    const std::string& output_path,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<std::vector<double>>& dist) {

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write similarity distance matrix: " + output_path);
    }

    out << std::fixed << std::setprecision(6);
    out << "allele_id";
    for (const auto& allele : unique_alleles) {
        out << '\t' << allele.allele_id;
    }
    out << '\n';

    for (std::size_t i = 0; i < unique_alleles.size(); ++i) {
        out << unique_alleles[i].allele_id;
        for (std::size_t j = 0; j < unique_alleles.size(); ++j) {
            out << '\t' << dist[i][j];
        }
        out << '\n';
    }
}

void write_similarity_distance_matrix_abs_tsv(
    const std::string& output_path,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<std::vector<int>>& dist_abs) {

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write absolute distance matrix: " + output_path);
    }

    out << "allele_id";
    for (const auto& allele : unique_alleles) {
        out << '\t' << allele.allele_id;
    }
    out << '\n';

    for (std::size_t i = 0; i < unique_alleles.size(); ++i) {
        out << unique_alleles[i].allele_id;
        for (std::size_t j = 0; j < unique_alleles.size(); ++j) {
            out << '\t' << dist_abs[i][j];
        }
        out << '\n';
    }
}

std::vector<double> compute_silhouette_scores(
    const std::vector<std::size_t>& cluster_of_unique,
    const std::vector<std::vector<double>>& dist_norm) {

    const std::size_t n = cluster_of_unique.size();
    std::vector<double> scores(n, 0.0);
    if (n <= 1 || dist_norm.size() != n) {
        return scores;
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> members;
    members.reserve(n * 2);
    for (std::size_t i = 0; i < n; ++i) {
        members[cluster_of_unique[i]].push_back(i);
    }
    if (members.size() <= 1) {
        return scores;
    }

    for (std::size_t i = 0; i < n; ++i) {
        const auto& own = members.at(cluster_of_unique[i]);
        double a = 0.0;
        if (own.size() > 1) {
            for (const std::size_t j : own) {
                if (j != i) {
                    a += dist_norm[i][j];
                }
            }
            a /= static_cast<double>(own.size() - 1);
        }

        double b = std::numeric_limits<double>::infinity();
        for (const auto& [cluster_id, idxs] : members) {
            if (cluster_id == cluster_of_unique[i] || idxs.empty()) {
                continue;
            }
            double mean = 0.0;
            for (const std::size_t j : idxs) {
                mean += dist_norm[i][j];
            }
            mean /= static_cast<double>(idxs.size());
            b = std::min(b, mean);
        }
        if (!std::isfinite(b)) {
            scores[i] = 0.0;
            continue;
        }
        const double denom = std::max(a, b);
        scores[i] = (denom <= 0.0) ? 0.0 : ((b - a) / denom);
    }
    return scores;
}

void write_similarity_cluster_stats_tsv(
    const std::string& output_path,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<std::size_t>& cluster_of_unique,
    const std::vector<std::vector<double>>& dist_norm,
    const std::vector<double>& silhouette_scores,
    const std::vector<std::vector<int>>& dist_abs,
    double min_similarity) {

    std::unordered_map<std::size_t, std::vector<std::size_t>> members;
    for (std::size_t i = 0; i < cluster_of_unique.size(); ++i) {
        members[cluster_of_unique[i]].push_back(i);
    }

    std::vector<std::size_t> cluster_ids;
    cluster_ids.reserve(members.size());
    for (const auto& [cluster_id, _idxs] : members) {
        cluster_ids.push_back(cluster_id);
    }
    std::sort(cluster_ids.begin(), cluster_ids.end());

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write similarity cluster stats: " + output_path);
    }
    out << std::fixed << std::setprecision(6);
    out << "cluster_id\tmember_alleles\tpath_support_sum\t"
           "mean_silhouette\tpath_weighted_mean_silhouette\t"
           "mean_intra_distance_norm\tmean_intra_similarity\tmean_intra_edit_distance\t"
           "max_intra_distance_norm\tmax_intra_edit_distance\t"
           "mean_to_other_clusters_distance_norm\tmean_to_other_clusters_similarity\t"
           "mean_to_other_clusters_edit_distance\t"
           "min_to_other_clusters_distance_norm\tmax_to_other_clusters_distance_norm\t"
           "min_to_other_clusters_edit_distance\tmax_to_other_clusters_edit_distance\t"
           "nearest_other_cluster_id\tnearest_other_min_distance_norm\tnearest_other_min_similarity\t"
           "separation_margin_norm\tthreshold_margin_norm\tseparation_quality\n";

    for (const std::size_t cluster_id : cluster_ids) {
        const auto& idxs = members[cluster_id];
        double intra_sum_norm = 0.0;
        double intra_sum_abs = 0.0;
        std::size_t intra_pairs = 0;
        double intra_max_norm = 0.0;
        double intra_max_abs = 0.0;
        for (std::size_t a = 0; a < idxs.size(); ++a) {
            for (std::size_t b = a + 1; b < idxs.size(); ++b) {
                intra_sum_norm += dist_norm[idxs[a]][idxs[b]];
                intra_sum_abs += static_cast<double>(dist_abs[idxs[a]][idxs[b]]);
                intra_max_norm = std::max(intra_max_norm, dist_norm[idxs[a]][idxs[b]]);
                intra_max_abs = std::max(intra_max_abs, static_cast<double>(dist_abs[idxs[a]][idxs[b]]));
                ++intra_pairs;
            }
        }
        const double intra_mean_norm = safe_div(intra_sum_norm, static_cast<double>(intra_pairs));
        const double intra_mean_abs = safe_div(intra_sum_abs, static_cast<double>(intra_pairs));

        double between_cluster_sum_norm = 0.0;
        double between_cluster_sum_abs = 0.0;
        std::size_t between_cluster_count = 0;
        double between_cluster_min_norm = std::numeric_limits<double>::infinity();
        double between_cluster_max_norm = 0.0;
        double between_cluster_min_abs = std::numeric_limits<double>::infinity();
        double between_cluster_max_abs = 0.0;
        std::size_t nearest_other_cluster_id = 0;
        double nearest_other_min_norm = std::numeric_limits<double>::infinity();
        for (const std::size_t other_cluster : cluster_ids) {
            if (other_cluster == cluster_id) {
                continue;
            }
            const auto& other = members[other_cluster];
            double sum_norm = 0.0;
            double sum_abs = 0.0;
            std::size_t pairs = 0;
            double pair_min_norm = std::numeric_limits<double>::infinity();
            for (const std::size_t i : idxs) {
                for (const std::size_t j : other) {
                    sum_norm += dist_norm[i][j];
                    sum_abs += static_cast<double>(dist_abs[i][j]);
                    pair_min_norm = std::min(pair_min_norm, dist_norm[i][j]);
                    ++pairs;
                }
            }
            if (pairs == 0) {
                continue;
            }
            const double m_norm = sum_norm / static_cast<double>(pairs);
            const double m_abs = sum_abs / static_cast<double>(pairs);
            between_cluster_sum_norm += m_norm;
            between_cluster_sum_abs += m_abs;
            ++between_cluster_count;
            between_cluster_min_norm = std::min(between_cluster_min_norm, m_norm);
            between_cluster_max_norm = std::max(between_cluster_max_norm, m_norm);
            between_cluster_min_abs = std::min(between_cluster_min_abs, m_abs);
            between_cluster_max_abs = std::max(between_cluster_max_abs, m_abs);
            if (pair_min_norm < nearest_other_min_norm) {
                nearest_other_min_norm = pair_min_norm;
                nearest_other_cluster_id = other_cluster;
            }
        }
        const double between_cluster_mean_norm = safe_div(
            between_cluster_sum_norm,
            static_cast<double>(between_cluster_count));
        const double between_cluster_mean_abs = safe_div(
            between_cluster_sum_abs,
            static_cast<double>(between_cluster_count));
        if (!std::isfinite(between_cluster_min_norm)) {
            between_cluster_min_norm = 0.0;
        }
        if (!std::isfinite(between_cluster_min_abs)) {
            between_cluster_min_abs = 0.0;
        }
        if (!std::isfinite(nearest_other_min_norm)) {
            nearest_other_min_norm = 0.0;
            nearest_other_cluster_id = 0;
        }

        const double max_norm_dist = std::clamp(1.0 - min_similarity, 0.0, 1.0);
        const double separation_margin_norm = nearest_other_min_norm - intra_max_norm;
        const double threshold_margin_norm = nearest_other_min_norm - max_norm_dist;
        std::string separation_quality = "single_cluster";
        if (cluster_ids.size() > 1) {
            if (threshold_margin_norm < 0.0 || separation_margin_norm < 0.0) {
                separation_quality = "poor";
            } else if (threshold_margin_norm < 0.02 || separation_margin_norm < 0.02) {
                separation_quality = "borderline";
            } else {
                separation_quality = "good";
            }
        }

        std::size_t path_support_sum = 0;
        double silhouette_sum = 0.0;
        double weighted_silhouette_sum = 0.0;
        for (const std::size_t idx : idxs) {
            path_support_sum += unique_alleles[idx].path_support;
            silhouette_sum += silhouette_scores[idx];
            weighted_silhouette_sum +=
                silhouette_scores[idx] * static_cast<double>(unique_alleles[idx].path_support);
        }

        out << cluster_id << '\t'
            << idxs.size() << '\t'
            << path_support_sum << '\t'
            << safe_div(silhouette_sum, static_cast<double>(idxs.size())) << '\t'
            << safe_div(weighted_silhouette_sum, static_cast<double>(path_support_sum)) << '\t'
            << intra_mean_norm << '\t'
            << (1.0 - intra_mean_norm) << '\t'
            << intra_mean_abs << '\t'
            << intra_max_norm << '\t'
            << intra_max_abs << '\t'
            << between_cluster_mean_norm << '\t'
            << (1.0 - between_cluster_mean_norm) << '\t'
            << between_cluster_mean_abs << '\t'
            << between_cluster_min_norm << '\t'
            << between_cluster_max_norm << '\t'
            << between_cluster_min_abs << '\t'
            << between_cluster_max_abs << '\t'
            << nearest_other_cluster_id << '\t'
            << nearest_other_min_norm << '\t'
            << (1.0 - nearest_other_min_norm) << '\t'
            << separation_margin_norm << '\t'
            << threshold_margin_norm << '\t'
            << separation_quality << '\n';
    }
}

std::size_t node_length_bp(
    const std::unordered_map<std::string, Node>& nodes,
    const std::string& node_id) {

    const auto it = nodes.find(node_id);
    if (it == nodes.end() || it->second.sequence.empty()) {
        return 1;
    }
    return std::max<std::size_t>(1, it->second.sequence.size());
}

std::string allele_node_orientation_cell(
    const UniqueAllele& allele,
    const std::string& node_id) {

    std::string cell;
    for (const PathStep& step : allele.steps) {
        if (step.node_id == node_id) {
            cell.push_back(step.reverse ? '-' : '+');
        }
    }
    return cell.empty() ? std::string("|") : cell;
}

std::vector<std::string> cluster_signature_node_order(
    const Cluster& cluster,
    const std::vector<UniqueAllele>& unique_alleles) {

    std::vector<std::string> order;
    std::unordered_set<std::string> seen;

    auto add_nodes = [&](std::size_t allele_idx) {
        if (allele_idx >= unique_alleles.size()) {
            return;
        }
        for (const PathStep& step : unique_alleles[allele_idx].steps) {
            if (seen.insert(step.node_id).second) {
                order.push_back(step.node_id);
            }
        }
    };

    add_nodes(cluster.representative_idx);
    for (const std::size_t idx : cluster.member_unique_idxs) {
        add_nodes(idx);
    }
    return order;
}

void write_similarity_cluster_signature_tables(
    const std::string& output_dir,
    const std::unordered_map<std::string, Node>& nodes,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<Cluster>& clusters) {

    std::filesystem::create_directories(output_dir);
    for (const Cluster& cluster : clusters) {
        const std::string output_path = output_dir + "/cluster_" + std::to_string(cluster.cluster_id) + ".tsv";
        std::ofstream out(output_path);
        if (!out) {
            throw std::runtime_error("Failed to write cluster signature table: " + output_path);
        }

        out << "node_id";
        for (const std::size_t idx : cluster.member_unique_idxs) {
            out << "\tallele_" << unique_alleles[idx].allele_id;
        }
        out << "\tnode_length_bp\n";

        const std::vector<std::string> node_order = cluster_signature_node_order(cluster, unique_alleles);
        for (const std::string& node_id : node_order) {
            out << node_id;
            for (const std::size_t idx : cluster.member_unique_idxs) {
                out << '\t' << allele_node_orientation_cell(unique_alleles[idx], node_id);
            }
            out << '\t' << node_length_bp(nodes, node_id) << '\n';
        }
    }
}

} // namespace

SimilarityReportPaths similarity_report_paths_for_bubble(
    const std::string& output_dir,
    std::size_t bubble_id) {

    SimilarityReportPaths out;
    out.bubble_dir = output_dir + "/bubble_" + std::to_string(bubble_id);
    out.alleles_tsv = out.bubble_dir + "/alleles.tsv";
    out.matrix_norm_tsv = out.bubble_dir + "/distance_matrix_norm.tsv";
    out.matrix_abs_tsv = out.bubble_dir + "/distance_matrix_abs.tsv";
    out.stats_tsv = out.bubble_dir + "/cluster_stats.tsv";
    out.cluster_signatures_dir = out.bubble_dir + "/cluster_signatures";
    return out;
}

void write_similarity_reports_for_bubble(
    const SimilarityReportPaths& paths,
    const Bubble& bubble,
    const Graph& graph,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<Cluster>& clusters,
    const std::vector<std::size_t>& cluster_of_unique,
    const DistanceMatrices& dists,
    double min_similarity) {

    if (unique_alleles.empty()) {
        return;
    }

    std::filesystem::create_directories(paths.bubble_dir);
    const std::vector<double> silhouette_scores =
        compute_silhouette_scores(cluster_of_unique, dists.norm);

    write_similarity_alleles_tsv(paths.alleles_tsv, bubble, unique_alleles, cluster_of_unique);
    write_similarity_distance_matrix_tsv(paths.matrix_norm_tsv, unique_alleles, dists.norm);
    write_similarity_distance_matrix_abs_tsv(paths.matrix_abs_tsv, unique_alleles, dists.abs);
    write_similarity_cluster_stats_tsv(
        paths.stats_tsv,
        unique_alleles,
        cluster_of_unique,
        dists.norm,
        silhouette_scores,
        dists.abs,
        min_similarity);
    write_similarity_cluster_signature_tables(
        paths.cluster_signatures_dir,
        graph.nodes,
        unique_alleles,
        clusters);
}

} // namespace panvar
