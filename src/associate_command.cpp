#include "panvar/associate_command.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "panvar/gz_reader.hpp"

namespace panvar {
namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> out;
    std::size_t start = 0;
    for (std::size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == delim) {
            out.push_back(s.substr(start, i - start));
            start = i + 1;
        }
    }
    return out;
}

std::string trim(const std::string& s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

bool is_na_token(const std::string& t) {
    std::string u = trim(t);
    for (char& c : u) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return u.empty() || u == "na" || u == "nan" || u == "." || u == "-9" || u == "null";
}

// Parse a numeric cell; returns NaN for missing tokens.
double parse_num(const std::string& t) {
    if (is_na_token(t)) return kNaN;
    try {
        return std::stod(trim(t));
    } catch (const std::exception&) {
        return kNaN;
    }
}

// Two-sided Wald p-value from a z statistic. Computed as the upper tail directly via erfc
// (p = erfc(|z|/sqrt2)) to avoid the catastrophic cancellation of 1 - Phi(|z|) for large |z|.
// Floored at 1e-300 so -log10(p) stays finite for plotting at extreme significance.
double wald_p(double z) {
    if (!std::isfinite(z)) return kNaN;
    const double p = std::erfc(std::fabs(z) / std::sqrt(2.0));
    return p < 1e-300 ? 1e-300 : p;
}

// Solve A x = b and also return A^{-1} (Gauss-Jordan, partial pivot). A is p x p (row-major),
// small (p <= ~12). Returns false if singular.
bool solve_and_invert(std::vector<double> A, const std::vector<double>& b,
                      std::size_t p, std::vector<double>& x, std::vector<double>& inv) {
    // augmented [A | I]
    std::vector<double> M(p * 2 * p, 0.0);
    for (std::size_t i = 0; i < p; ++i) {
        for (std::size_t j = 0; j < p; ++j) M[i * 2 * p + j] = A[i * p + j];
        M[i * 2 * p + p + i] = 1.0;
    }
    for (std::size_t col = 0; col < p; ++col) {
        std::size_t piv = col;
        double best = std::fabs(M[col * 2 * p + col]);
        for (std::size_t r = col + 1; r < p; ++r) {
            const double v = std::fabs(M[r * 2 * p + col]);
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-12) return false;
        if (piv != col) for (std::size_t j = 0; j < 2 * p; ++j) std::swap(M[col * 2 * p + j], M[piv * 2 * p + j]);
        const double d = M[col * 2 * p + col];
        for (std::size_t j = 0; j < 2 * p; ++j) M[col * 2 * p + j] /= d;
        for (std::size_t r = 0; r < p; ++r) {
            if (r == col) continue;
            const double f = M[r * 2 * p + col];
            if (f == 0.0) continue;
            for (std::size_t j = 0; j < 2 * p; ++j) M[r * 2 * p + j] -= f * M[col * 2 * p + j];
        }
    }
    inv.assign(p * p, 0.0);
    for (std::size_t i = 0; i < p; ++i)
        for (std::size_t j = 0; j < p; ++j) inv[i * p + j] = M[i * 2 * p + p + j];
    x.assign(p, 0.0);
    for (std::size_t i = 0; i < p; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < p; ++j) s += inv[i * p + j] * b[j];
        x[i] = s;
    }
    return true;
}

struct FitResult {
    bool ok = false;
    double beta = kNaN;  // genotype coefficient (column 1)
    double se = kNaN;
    double z = kNaN;
    double p = kNaN;
};

// X is n x p row-major (col 0 = intercept, col 1 = genotype, cols 2.. = covariates). Returns the
// Wald test on the genotype coefficient (index 1).
FitResult fit_linear(const std::vector<double>& X, const std::vector<double>& y, std::size_t n, std::size_t p) {
    FitResult r;
    if (n <= p + 1) return r;
    std::vector<double> XtX(p * p, 0.0), Xty(p, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t a = 0; a < p; ++a) {
            const double xa = X[i * p + a];
            Xty[a] += xa * y[i];
            for (std::size_t b = a; b < p; ++b) XtX[a * p + b] += xa * X[i * p + b];
        }
    }
    for (std::size_t a = 0; a < p; ++a) for (std::size_t b = 0; b < a; ++b) XtX[a * p + b] = XtX[b * p + a];
    std::vector<double> beta, inv;
    if (!solve_and_invert(XtX, Xty, p, beta, inv)) return r;
    double rss = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double yhat = 0.0;
        for (std::size_t a = 0; a < p; ++a) yhat += X[i * p + a] * beta[a];
        const double e = y[i] - yhat;
        rss += e * e;
    }
    const double sigma2 = rss / static_cast<double>(n - p);
    const double var_b1 = sigma2 * inv[1 * p + 1];
    if (!(var_b1 > 0.0)) return r;
    r.beta = beta[1];
    r.se = std::sqrt(var_b1);
    r.z = r.beta / r.se;
    r.p = wald_p(r.z);
    r.ok = std::isfinite(r.p);
    return r;
}

FitResult fit_logistic(const std::vector<double>& X, const std::vector<double>& y, std::size_t n, std::size_t p) {
    FitResult r;
    if (n <= p + 1) return r;
    std::vector<double> beta(p, 0.0), inv;
    for (int iter = 0; iter < 50; ++iter) {
        std::vector<double> XtWX(p * p, 0.0), XtWz(p, 0.0);
        bool bad = false;
        for (std::size_t i = 0; i < n; ++i) {
            double eta = 0.0;
            for (std::size_t a = 0; a < p; ++a) eta += X[i * p + a] * beta[a];
            const double mu = 1.0 / (1.0 + std::exp(-eta));
            const double w = std::max(mu * (1.0 - mu), 1e-9);
            const double z = eta + (y[i] - mu) / w;
            for (std::size_t a = 0; a < p; ++a) {
                const double xa = X[i * p + a];
                XtWz[a] += xa * w * z;
                for (std::size_t b = a; b < p; ++b) XtWX[a * p + b] += xa * w * X[i * p + b];
            }
        }
        for (std::size_t a = 0; a < p; ++a) for (std::size_t b = 0; b < a; ++b) XtWX[a * p + b] = XtWX[b * p + a];
        std::vector<double> nb;
        if (!solve_and_invert(XtWX, XtWz, p, nb, inv)) { bad = true; }
        if (bad) return r;
        double delta = 0.0;
        for (std::size_t a = 0; a < p; ++a) { delta += std::fabs(nb[a] - beta[a]); beta[a] = nb[a]; }
        if (delta < 1e-8) break;
    }
    const double var_b1 = inv.empty() ? kNaN : inv[1 * p + 1];
    if (!(var_b1 > 0.0)) return r;
    r.beta = beta[1];
    r.se = std::sqrt(var_b1);
    r.z = r.beta / r.se;
    r.p = wald_p(r.z);
    r.ok = std::isfinite(r.p);
    return r;
}

// Benjamini-Hochberg q-values for a vector of p-values (NaN entries get NaN q).
std::vector<double> bh_qvalues(const std::vector<double>& pv) {
    std::vector<std::size_t> idx;
    for (std::size_t i = 0; i < pv.size(); ++i) if (std::isfinite(pv[i])) idx.push_back(i);
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return pv[a] < pv[b]; });
    const std::size_t m = idx.size();
    std::vector<double> q(pv.size(), kNaN);
    double prev = 1.0;
    for (std::size_t k = m; k-- > 0;) {
        const std::size_t i = idx[k];
        double val = pv[i] * static_cast<double>(m) / static_cast<double>(k + 1);
        if (val > 1.0) val = 1.0;
        prev = std::min(prev, val);
        q[i] = prev;
    }
    return q;
}

struct Options {
    std::string genotypes;
    std::string samples;
    std::string feature_annot;
    std::string phenotype;
    std::string out_prefix = "associate";
    double min_maf = 0.01;
    std::string model = "auto";  // auto|linear|logistic
    bool quiet = false;
};

void print_help() {
    std::cout
        << "Usage: panvar associate --genotypes <bimbam.gz> --samples <samples.txt[.gz]> \\\n"
        << "                        --phenotype <table.tsv> -o <prefix> [options]\n\n"
        << "GWAS on a BIMBAM dosage matrix (from `describe`) vs a phenotype/covariate table.\n\n"
        << "Required:\n"
        << "  --genotypes <path>     BIMBAM mean-genotype dosage (describe bimbam_{kmers,graph}.bimbam.gz)\n"
        << "  --samples <path>       sample (column) order, one per line (describe bimbam.samples.txt.gz)\n"
        << "  --phenotype <path>     TSV: sample <tab> phenotype [<tab> covariate1 ...]; header required.\n"
        << "                         Phenotype/covariate cells may be NA (NA-phenotype samples are dropped).\n"
        << "  -o, --out-prefix <p>   output prefix\n\n"
        << "Options:\n"
        << "  --feature-annot <path> describe feature_annot.tsv.gz (adds layer/bubbles/nodes to output)\n"
        << "  --min-maf <X>          drop features whose minor (non-modal) genotype frequency < X (default 0.01)\n"
        << "  --model <m>            auto|linear|logistic (default auto: binary phenotype -> logistic)\n"
        << "  -q, --quiet            less logging\n\n"
        << "Outputs: <prefix>.assoc.tsv (per-feature beta/se/p/p_bonf/q_bh) and <prefix>.summary.tsv\n"
        << "  (n_tests, Bonferroni threshold 0.05/n_tests, significant counts). Plot with scripts/plot_gwas.R.\n";
}

struct FeatureAnnot {
    std::string layer = ".";
    std::string bubbles = ".";
    std::string nodes = ".";
};

}  // namespace

int run_associate_command(const std::vector<std::string>& args) {
    if (args.empty() || args[0] == "-h" || args[0] == "--help") { print_help(); return args.empty() ? 1 : 0; }
    Options opt;
    auto need = [&](std::size_t& i) -> std::string {
        if (i + 1 >= args.size()) throw std::runtime_error("missing value for " + args[i]);
        return args[++i];
    };
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& a = args[i];
        if (a == "--genotypes") opt.genotypes = need(i);
        else if (a == "--samples") opt.samples = need(i);
        else if (a == "--feature-annot") opt.feature_annot = need(i);
        else if (a == "--phenotype") opt.phenotype = need(i);
        else if (a == "-o" || a == "--out-prefix") opt.out_prefix = need(i);
        else if (a == "--min-maf") opt.min_maf = std::stod(need(i));
        else if (a == "--model") opt.model = need(i);
        else if (a == "-q" || a == "--quiet") opt.quiet = true;
        else throw std::runtime_error("Unknown option for associate: " + a);
    }
    if (opt.genotypes.empty() || opt.samples.empty() || opt.phenotype.empty())
        throw std::runtime_error("associate requires --genotypes, --samples, --phenotype");
    if (opt.model != "auto" && opt.model != "linear" && opt.model != "logistic")
        throw std::runtime_error("--model must be auto|linear|logistic");

    // ---- sample (genotype column) order ----
    std::vector<std::string> geno_samples;
    {
        GzLineReader r(opt.samples);
        if (!r.ok()) throw std::runtime_error("cannot open --samples: " + opt.samples);
        std::string line;
        while (r.getline(line)) { line = trim(line); if (!line.empty()) geno_samples.push_back(line); }
    }
    if (geno_samples.empty()) throw std::runtime_error("no samples in " + opt.samples);

    // ---- phenotype + covariates: sample -> (phenotype, covariates[]) ----
    std::vector<std::string> covar_names;
    std::unordered_map<std::string, double> pheno;
    std::unordered_map<std::string, std::vector<double>> covars;
    {
        GzLineReader r(opt.phenotype);
        if (!r.ok()) throw std::runtime_error("cannot open --phenotype: " + opt.phenotype);
        std::string line;
        bool header = true;
        std::size_t ncov = 0;
        while (r.getline(line)) {
            if (line.empty()) continue;
            std::vector<std::string> f = split(line, '\t');
            if (f.size() < 2) f = split(line, ',');
            if (header) {
                header = false;
                for (std::size_t j = 2; j < f.size(); ++j) covar_names.push_back(trim(f[j]));
                ncov = covar_names.size();
                continue;
            }
            const std::string s = trim(f[0]);
            if (s.empty()) continue;
            pheno[s] = (f.size() > 1) ? parse_num(f[1]) : kNaN;
            std::vector<double> cv(ncov, kNaN);
            for (std::size_t j = 0; j < ncov; ++j) if (j + 2 < f.size()) cv[j] = parse_num(f[j + 2]);
            covars[s] = std::move(cv);
        }
    }

    // ---- usable samples: present in genotype order AND with non-NA phenotype + all covariates ----
    // sample_used[g] = index into the compact phenotype/covariate arrays, or -1.
    const std::size_t ncov = covar_names.size();
    std::vector<long long> col_to_used(geno_samples.size(), -1);
    std::vector<double> y;          // phenotype for used samples
    std::vector<std::vector<double>> Z;  // covariates for used samples
    for (std::size_t c = 0; c < geno_samples.size(); ++c) {
        const auto pit = pheno.find(geno_samples[c]);
        if (pit == pheno.end() || !std::isfinite(pit->second)) continue;
        const auto cit = covars.find(geno_samples[c]);
        std::vector<double> cv = (cit != covars.end()) ? cit->second : std::vector<double>(ncov, kNaN);
        bool ok = true;
        for (double v : cv) if (!std::isfinite(v)) { ok = false; break; }
        if (!ok) continue;
        col_to_used[c] = static_cast<long long>(y.size());
        y.push_back(pit->second);
        Z.push_back(std::move(cv));
    }
    const std::size_t n_used = y.size();
    if (n_used < ncov + 3) throw std::runtime_error("too few usable samples after NA filtering (" +
                                                    std::to_string(n_used) + ")");

    // ---- phenotype type / model ----
    bool binary = true;
    for (double v : y) if (v != 0.0 && v != 1.0) { binary = false; break; }
    std::string model = opt.model;
    if (model == "auto") model = binary ? "logistic" : "linear";
    if (model == "logistic" && !binary)
        throw std::runtime_error("--model logistic but phenotype is not binary (0/1)");

    // ---- optional feature annotation ----
    std::unordered_map<std::string, FeatureAnnot> annot;
    if (!opt.feature_annot.empty()) {
        GzLineReader r(opt.feature_annot);
        if (!r.ok()) throw std::runtime_error("cannot open --feature-annot: " + opt.feature_annot);
        std::string line; bool header = true;
        while (r.getline(line)) {
            if (line.empty()) continue;
            if (header) { header = false; continue; }
            std::vector<std::string> f = split(line, '\t');
            if (f.size() < 5) continue;
            FeatureAnnot a; a.layer = f[1]; a.bubbles = f[3]; a.nodes = f[4];
            annot[f[0]] = std::move(a);
        }
    }

    // ---- per-feature association over the BIMBAM rows ----
    struct Row {
        std::string id, layer, bubbles, nodes;
        std::size_t n = 0;
        double minor_freq = kNaN, beta = kNaN, se = kNaN, z = kNaN, p = kNaN;
    };
    std::vector<Row> rows;
    std::size_t n_geno_rows = 0, n_dropped_maf = 0, n_dropped_fit = 0;
    const std::size_t p_dim = 2 + ncov;  // intercept + genotype + covariates

    GzLineReader gr(opt.genotypes);
    if (!gr.ok()) throw std::runtime_error("cannot open --genotypes: " + opt.genotypes);
    std::string line;
    while (gr.getline(line)) {
        if (line.empty()) continue;
        // BIMBAM mean genotype: id, A, B, dose1, dose2, ...   (comma-separated, possibly space-padded)
        std::vector<std::string> f = split(line, ',');
        if (f.size() < 3 + geno_samples.size()) continue;  // malformed / wrong sample count
        ++n_geno_rows;
        const std::string id = trim(f[0]);

        // collect genotype dosage for used samples (skip NA genotypes)
        std::vector<double> g; g.reserve(n_used);
        std::vector<double> yy; yy.reserve(n_used);
        std::vector<std::vector<double>> zz; zz.reserve(n_used);
        for (std::size_t c = 0; c < geno_samples.size(); ++c) {
            if (col_to_used[c] < 0) continue;
            const double dose = parse_num(f[3 + c]);
            if (!std::isfinite(dose)) continue;  // genuinely missing (NA = non-traversing)
            const std::size_t u = static_cast<std::size_t>(col_to_used[c]);
            g.push_back(dose);
            yy.push_back(y[u]);
            zz.push_back(Z[u]);
        }
        const std::size_t n = g.size();
        if (n < p_dim + 1) { ++n_dropped_fit; continue; }

        // minor (non-modal) genotype frequency: 1 - freq(most common rounded dosage). Works for both
        // presence/absence (0/1) and CN dosage (drops invariant / cohort-rare features).
        std::unordered_map<long long, std::size_t> cat;
        for (double v : g) ++cat[std::llround(v)];
        std::size_t mode = 0;
        for (const auto& kv : cat) mode = std::max(mode, kv.second);
        const double minor_freq = 1.0 - static_cast<double>(mode) / static_cast<double>(n);
        if (minor_freq < opt.min_maf) { ++n_dropped_maf; continue; }

        // design matrix
        std::vector<double> X(n * p_dim);
        for (std::size_t i = 0; i < n; ++i) {
            X[i * p_dim + 0] = 1.0;
            X[i * p_dim + 1] = g[i];
            for (std::size_t j = 0; j < ncov; ++j) X[i * p_dim + 2 + j] = zz[i][j];
        }
        const FitResult fr = (model == "logistic")
            ? fit_logistic(X, yy, n, p_dim)
            : fit_linear(X, yy, n, p_dim);
        if (!fr.ok) { ++n_dropped_fit; continue; }

        Row row;
        row.id = id;
        if (auto it = annot.find(id); it != annot.end()) {
            row.layer = it->second.layer; row.bubbles = it->second.bubbles; row.nodes = it->second.nodes;
        } else { row.layer = "."; row.bubbles = "."; row.nodes = "."; }
        row.n = n; row.minor_freq = minor_freq;
        row.beta = fr.beta; row.se = fr.se; row.z = fr.z; row.p = fr.p;
        rows.push_back(std::move(row));
    }

    // ---- multiple testing over the ACTUALLY TESTED features (not genome-wide) ----
    const std::size_t n_tests = rows.size();
    std::vector<double> pv(n_tests);
    for (std::size_t i = 0; i < n_tests; ++i) pv[i] = rows[i].p;
    const std::vector<double> qv = bh_qvalues(pv);
    const double bonf_threshold = n_tests > 0 ? 0.05 / static_cast<double>(n_tests) : kNaN;
    std::size_t n_sig_bonf = 0, n_sig_fdr = 0;
    for (std::size_t i = 0; i < n_tests; ++i) {
        if (std::isfinite(rows[i].p) && rows[i].p < bonf_threshold) ++n_sig_bonf;
        if (std::isfinite(qv[i]) && qv[i] < 0.05) ++n_sig_fdr;
    }

    // ---- write assoc.tsv (sorted by p) ----
    std::vector<std::size_t> order(n_tests);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        if (!std::isfinite(rows[a].p)) return false;
        if (!std::isfinite(rows[b].p)) return true;
        return rows[a].p < rows[b].p;
    });
    auto fmt = [](double v) {
        if (!std::isfinite(v)) return std::string("NA");
        std::ostringstream o; o << v; return o.str();
    };
    const std::string effect = (model == "logistic") ? "log_or" : "beta";
    {
        std::ofstream out(opt.out_prefix + ".assoc.tsv");
        if (!out) throw std::runtime_error("cannot write " + opt.out_prefix + ".assoc.tsv");
        out << "feature_id\tlayer\tbubbles\tnodes\tn\tminor_freq\t" << effect
            << "\tse\tz\tp\tp_bonf\tq_bh\n";
        for (std::size_t k = 0; k < n_tests; ++k) {
            const std::size_t i = order[k];
            const Row& r = rows[i];
            const double p_bonf = std::isfinite(r.p) ? std::min(1.0, r.p * static_cast<double>(n_tests)) : kNaN;
            out << r.id << '\t' << r.layer << '\t' << r.bubbles << '\t' << r.nodes << '\t'
                << r.n << '\t' << fmt(r.minor_freq) << '\t' << fmt(r.beta) << '\t' << fmt(r.se)
                << '\t' << fmt(r.z) << '\t' << fmt(r.p) << '\t' << fmt(p_bonf) << '\t' << fmt(qv[i]) << '\n';
        }
    }
    // ---- summary (also the per-plot threshold lines) ----
    {
        std::ofstream out(opt.out_prefix + ".summary.tsv");
        out << "key\tvalue\n"
            << "model\t" << model << '\n'
            << "phenotype_type\t" << (binary ? "binary" : "quantitative") << '\n'
            << "covariates\t" << (ncov ? std::to_string(ncov) : std::string("0")) << '\n'
            << "samples_used\t" << n_used << '\n'
            << "genotype_rows\t" << n_geno_rows << '\n'
            << "features_tested\t" << n_tests << '\n'
            << "dropped_min_maf\t" << n_dropped_maf << '\n'
            << "dropped_fit\t" << n_dropped_fit << '\n'
            << "bonferroni_threshold\t" << fmt(bonf_threshold) << '\n'
            << "nominal_threshold\t0.05\n"
            << "significant_bonferroni\t" << n_sig_bonf << '\n'
            << "significant_fdr05\t" << n_sig_fdr << '\n';
    }

    if (!opt.quiet) {
        std::cerr << "[associate] model=" << model << " (" << (binary ? "binary" : "quantitative")
                  << "), samples=" << n_used << ", covariates=" << ncov << "\n"
                  << "[associate] features: tested=" << n_tests << " dropped(maf)=" << n_dropped_maf
                  << " dropped(fit)=" << n_dropped_fit << "\n"
                  << "[associate] Bonferroni threshold (0.05/" << n_tests << ") = " << fmt(bonf_threshold)
                  << "; significant: Bonferroni=" << n_sig_bonf << " FDR<0.05=" << n_sig_fdr << "\n"
                  << "[associate] wrote " << opt.out_prefix << ".assoc.tsv + .summary.tsv\n";
    }
    return 0;
}

}  // namespace panvar
