#include "panvar/genotype_model.hpp"
#include "panvar/ls_hmm.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void require_close(double observed, double expected, double tolerance, const std::string& message) {
    if (!std::isfinite(observed) || !std::isfinite(expected) ||
        std::abs(observed - expected) > tolerance) {
        throw std::runtime_error(message + ": observed=" + std::to_string(observed) +
                                 " expected=" + std::to_string(expected));
    }
}

void require_throws(const std::function<void()>& operation, const std::string& message) {
    try {
        operation();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(message);
}

long double brute_haploid_prior(
    const panvar::PanelAlleleMatrix& panel,
    const panvar::MosaicPath& mosaic,
    const std::vector<double>& switch_probability) {
    const std::size_t n = panel.template_count();
    long double total = 0.0L;

    std::function<void(std::size_t, std::size_t, long double)> extend;
    extend = [&](std::size_t block, std::size_t previous, long double probability) {
        if (block == panel.block_count()) {
            total += probability;
            return;
        }
        for (std::size_t h = 0; h < n; ++h) {
            if (panel.template_alleles[h][block] != mosaic.alleles[block]) continue;
            const long double r = switch_probability[block - 1];
            const long double transition =
                (h == previous ? 1.0L - r : 0.0L) + r / static_cast<long double>(n);
            extend(block + 1, h, probability * transition);
        }
    };

    for (std::size_t h = 0; h < n; ++h) {
        if (panel.template_alleles[h][0] != mosaic.alleles[0]) continue;
        extend(1, h, 1.0L / static_cast<long double>(n));
    }
    return total;
}

long double brute_best_haploid_history(
    const panvar::PanelAlleleMatrix& panel,
    const panvar::MosaicPath& mosaic,
    const std::vector<double>& switch_probability) {
    const std::size_t n = panel.template_count();
    long double best = 0.0L;
    std::function<void(std::size_t, std::size_t, long double)> extend;
    extend = [&](std::size_t block, std::size_t previous, long double probability) {
        if (block == panel.block_count()) {
            best = std::max(best, probability);
            return;
        }
        for (std::size_t h = 0; h < n; ++h) {
            if (panel.template_alleles[h][block] != mosaic.alleles[block]) continue;
            const long double r = switch_probability[block - 1];
            const long double transition =
                (h == previous ? 1.0L - r : 0.0L) + r / static_cast<long double>(n);
            extend(block + 1, h, probability * transition);
        }
    };
    for (std::size_t h = 0; h < n; ++h) {
        if (panel.template_alleles[h][0] != mosaic.alleles[0]) continue;
        extend(1, h, 1.0L / static_cast<long double>(n));
    }
    return best;
}

panvar::SequencingRead make_read(const std::string& name, std::size_t length) {
    return {name, std::string(length, 'A'), std::string(length, 'I')};
}

void test_library_model() {
    panvar::LibraryModel model;
    model.insert_mean_bp = 350.0;
    model.insert_sd_bp = 50.0;
    model.support_sigmas = 6.0;
    model.provenance = "fixture: declared";

    const panvar::ReadPair overlapping{make_read("r/1", 150), make_read("r/2", 150)};
    const panvar::InsertSupport support = panvar::fragment_insert_support(overlapping, model);
    require(support.min_bp == 150 && support.max_bp == 650,
            "overlapping mates must use max(read lengths), giving support [150,650]");

    const panvar::DiscreteInsertDistribution distribution =
        panvar::build_fragment_insert_distribution(overlapping, model);
    double probability_sum = 0.0;
    for (std::size_t length = support.min_bp; length <= support.max_bp; ++length) {
        probability_sum += std::exp(distribution.log_probability(length));
    }
    require_close(probability_sum, 1.0, 2e-14,
                  "discrete insert probabilities must sum to one on declared support");
    require(std::isinf(distribution.log_probability(149)) &&
                distribution.log_probability(149) < 0.0,
            "insert below physical support must have zero probability");

    const panvar::ReadPair variable{make_read("v/1", 100), make_read("v/2", 250)};
    const panvar::InsertSupport variable_support = panvar::fragment_insert_support(variable, model);
    require(variable_support.min_bp == 250 && variable_support.max_bp == 650,
            "insert floor must be fragment-specific for variable read lengths");

    const panvar::ReadPair unsupported{make_read("u/1", 700), make_read("u/2", 150)};
    require(panvar::fragment_insert_support(unsupported, model).empty(),
            "a read longer than declared insert support must be an explicit unsupported state");
    require(panvar::build_fragment_insert_distribution(unsupported, model)
                .log_probability_by_offset.empty(),
            "unsupported insert distribution must not expose a partial table");

    panvar::ReadPair bad_quality = overlapping;
    bad_quality.first.qualities.pop_back();
    require_throws([&] { panvar::validate_read_pair(bad_quality); },
                   "base/quality length mismatch must refuse");

    panvar::LibraryModel no_provenance = model;
    no_provenance.provenance.clear();
    require_throws([&] { panvar::validate_library_model(no_provenance); },
                   "an unprovenanced library model must refuse");
}

void test_li_stephens_prior() {
    panvar::PanelAlleleMatrix panel;
    panel.template_names = {"h0", "h1", "h2"};
    panel.template_alleles = {
        {0, 0, 1, 0},
        {0, 1, 1, 1},
        {1, 0, 0, 1},
    };
    panvar::MosaicPath recombinant{{0, 0, 1, 1}};

    for (const std::vector<double>& recombination :
         std::vector<std::vector<double>>{{0.0, 0.0, 0.0},
                                          {1.0, 1.0, 1.0},
                                          {0.05, 0.3, 0.9}}) {
        const long double brute = brute_haploid_prior(panel, recombinant, recombination);
        const double exact =
            panvar::log_haploid_mosaic_prior(panel, recombinant, recombination);
        if (brute == 0.0L) {
            require(std::isinf(exact) && exact < 0.0,
                    "an unreachable exhaustive mosaic must have exactly zero prior");
            continue;
        }
        require_close(exact, std::log(static_cast<double>(brute)), 2e-14,
                      "O(blocks*templates) mosaic prior must equal exhaustive template histories");
    }

    const std::vector<double> intermediate{0.05, 0.3, 0.9};
    const long double summed_history = brute_haploid_prior(panel, recombinant, intermediate);
    const long double best_history = brute_best_haploid_history(panel, recombinant, intermediate);
    require(summed_history > best_history * 1.01L,
            "fixture must distinguish exact marginalization from the best template history");
    const panvar::MosaicPath second{{1, 0, 0, 1}};
    const double first_prior =
        panvar::log_haploid_mosaic_prior(panel, recombinant, intermediate);
    const double second_prior =
        panvar::log_haploid_mosaic_prior(panel, second, intermediate);

    panvar::MosaicDiplotype heterozygous{recombinant, second};
    require_close(
        panvar::log_unordered_diplotype_prior(panel, heterozygous, intermediate),
        first_prior + second_prior + std::log(2.0),
        2e-14,
        "distinct unordered mosaics must aggregate both homologue orders exactly once");

    panvar::MosaicDiplotype homozygous{recombinant, recombinant};
    require_close(
        panvar::log_unordered_diplotype_prior(panel, homozygous, intermediate),
        2.0 * first_prior,
        2e-14,
        "identical mosaics must not receive a spurious log(2) multiplicity");

    panvar::MosaicPath absent{{9, 9, 9, 9}};
    const double absent_prior = panvar::log_haploid_mosaic_prior(panel, absent, intermediate);
    require(std::isinf(absent_prior) && absent_prior < 0.0,
            "a mosaic with an unavailable allele must have zero prior");

    require_throws(
        [&] { panvar::log_haploid_mosaic_prior(panel, recombinant, {0.2, 0.2}); },
        "missing boundary probability must refuse");
    require_throws(
        [&] { panvar::log_haploid_mosaic_prior(panel, recombinant, {0.2, -0.1, 0.2}); },
        "invalid boundary probability must refuse");
}

} // namespace

int main() {
    try {
        test_library_model();
        test_li_stephens_prior();
        std::cout << "genotype mosaic core: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype mosaic core: FAIL: " << error.what() << '\n';
        return 1;
    }
}
