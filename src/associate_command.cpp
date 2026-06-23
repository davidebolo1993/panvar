#include "panvar/associate_command.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

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

// Benjamini-Hochberg FDR q-values (NaN entries get NaN q). Benjamini & Hochberg 1995, JRSS-B.
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

// Genomic-inflation factor: median of the chi-square stats (z^2) divided by the median of a chi-square_1
// (0.4549364). lambda ~ 1 means a well-calibrated test; lambda > 1 flags residual structure/relatedness.
double genomic_lambda(const std::vector<double>& zs) {
    std::vector<double> chi;
    for (double z : zs) if (std::isfinite(z)) chi.push_back(z * z);
    if (chi.empty()) return kNaN;
    std::sort(chi.begin(), chi.end());
    const std::size_t m = chi.size();
    const double med = (m % 2) ? chi[m / 2] : 0.5 * (chi[m / 2 - 1] + chi[m / 2]);
    return med / 0.4549364;
}

// ---- LMM (EMMAX-style) ---------------------------------------------------------------------------
// Eigendecompose the kinship K once, rotate phenotype + covariates, estimate delta = sigma_e^2/sigma_g^2
// once under the null, then test each feature by GLS in the rotated space with weights 1/(d_i + delta).
// EMMAX: Kang et al. 2010, https://doi.org/10.1038/ng.548 ; cf. GEMMA, Zhou & Stephens 2012.
struct LmmNull {
    Eigen::VectorXd d;     // kinship eigenvalues (n)
    Eigen::MatrixXd U;     // kinship eigenvectors (n x n)
    Eigen::MatrixXd UtX;   // U^T [intercept | covariates]  (n x q)
    Eigen::VectorXd Uty;   // U^T phenotype  (n)
    double delta = 1.0;    // sigma_e^2 / sigma_g^2 (estimated once under the null)
};

// REML objective to MINIMIZE over delta (profile likelihood, additive constants dropped):
//   (n-q) ln(sigma_g^2) + sum_i ln(d_i + delta) + ln det(X^T W X),  W = diag(1/(d_i+delta)).
double lmm_reml_obj(const LmmNull& m, double delta) {
    const Eigen::Index n = m.d.size();
    const Eigen::Index q = m.UtX.cols();
    const Eigen::ArrayXd w = 1.0 / (m.d.array() + delta);
    const Eigen::MatrixXd Xw = m.UtX.array().colwise() * w;   // row i scaled by w_i
    const Eigen::MatrixXd A = m.UtX.transpose() * Xw;         // q x q
    const Eigen::VectorXd b = Xw.transpose() * m.Uty;         // q
    const Eigen::VectorXd beta = A.ldlt().solve(b);
    const Eigen::VectorXd resid = m.Uty - m.UtX * beta;
    const double rss = (resid.array().square() * w).sum();
    const double sigma2 = rss / static_cast<double>(n - q);
    if (!(sigma2 > 0.0)) return std::numeric_limits<double>::infinity();
    const double logdetV = (m.d.array() + delta).log().sum();
    const double logdetA = std::log(std::max(A.determinant(), 1e-300));
    return static_cast<double>(n - q) * std::log(sigma2) + logdetV + logdetA;
}

// Coarse log-grid + golden-section refine for the delta minimizing the REML objective.
double estimate_delta(const LmmNull& m) {
    double best_l = -5.0, best_f = std::numeric_limits<double>::infinity();
    for (int k = -50; k <= 50; ++k) {
        const double l = k / 10.0;            // log10 delta in [-5, 5]
        const double f = lmm_reml_obj(m, std::pow(10.0, l));
        if (f < best_f) { best_f = f; best_l = l; }
    }
    double a = best_l - 0.1, b = best_l + 0.1;
    const double gr = 0.6180339887;
    double c = b - gr * (b - a), e = a + gr * (b - a);
    double fc = lmm_reml_obj(m, std::pow(10.0, c)), fe = lmm_reml_obj(m, std::pow(10.0, e));
    for (int it = 0; it < 40; ++it) {
        if (fc < fe) { b = e; e = c; fe = fc; c = b - gr * (b - a); fc = lmm_reml_obj(m, std::pow(10.0, c)); }
        else         { a = c; c = e; fc = fe; e = a + gr * (b - a); fe = lmm_reml_obj(m, std::pow(10.0, e)); }
    }
    return std::pow(10.0, 0.5 * (a + b));
}

// Per-feature GLS Wald test in the rotated space (genotype is the appended last column).
FitResult lmm_test(const LmmNull& m, const Eigen::VectorXd& g) {
    FitResult r;
    const Eigen::Index n = m.d.size();
    const Eigen::Index q = m.UtX.cols();
    if (n <= q + 2) return r;
    const Eigen::VectorXd Utg = m.U.transpose() * g;
    Eigen::MatrixXd X(n, q + 1);
    X.leftCols(q) = m.UtX;
    X.col(q) = Utg;
    const Eigen::ArrayXd w = 1.0 / (m.d.array() + m.delta);
    const Eigen::MatrixXd Xw = X.array().colwise() * w;
    const Eigen::MatrixXd A = X.transpose() * Xw;            // (q+1) x (q+1)
    const Eigen::VectorXd b = Xw.transpose() * m.Uty;
    const Eigen::MatrixXd Ainv = A.inverse();
    const Eigen::VectorXd beta = Ainv * b;
    const Eigen::VectorXd resid = m.Uty - X * beta;
    const double rss = (resid.array().square() * w).sum();
    const double sigma2 = rss / static_cast<double>(n - (q + 1));
    const double var_b = sigma2 * Ainv(q, q);
    if (!(var_b > 0.0)) return r;
    r.beta = beta(q);
    r.se = std::sqrt(var_b);
    r.z = r.beta / r.se;
    r.p = wald_p(r.z);
    r.ok = std::isfinite(r.p);
    return r;
}

// Top-N eigenvectors of a symmetric PSD matrix by orthogonal (subspace) iteration -- cheap relative
// to a full decomposition; the span is what we need for PC covariates.
Eigen::MatrixXd top_eigvecs(const Eigen::MatrixXd& K, int N) {
    const Eigen::Index n = K.rows();
    const Eigen::Index k = std::min<Eigen::Index>(N, n);
    Eigen::MatrixXd Q = Eigen::MatrixXd::Random(n, k);
    Q = Eigen::HouseholderQR<Eigen::MatrixXd>(Q).householderQ() * Eigen::MatrixXd::Identity(n, k);
    for (int it = 0; it < 60; ++it) {
        Eigen::MatrixXd Z = K * Q;
        Q = Eigen::HouseholderQR<Eigen::MatrixXd>(Z).householderQ() * Eigen::MatrixXd::Identity(n, k);
    }
    return Q;
}

struct Options {
    std::string genotypes;
    std::string samples;
    std::string feature_annot;
    std::string node_genes;
    std::string phenotype;
    std::string out_prefix = "associate";
    std::string kinship;         // precomputed GRM file (--kinship)
    bool make_kinship = false;   // build GRM from the genotype matrix (--make-kinship)
    int pca = 0;                 // top-N kinship PCs as GLM covariates (--pca)
    double min_maf = 0.01;
    std::string model = "auto";  // auto|linear|logistic|lmm
    bool quiet = false;
};

// Extract every maximal run of digits from a provenance string and map each to a node id. Handles
// node ids ("4789"), edge keys ("4789+>4789+"), and ";"-joined k-mer node lists ("123;456").
std::vector<std::string> node_ids_in(const std::string& s) {
    std::vector<std::string> out;
    std::size_t i = 0;
    while (i < s.size()) {
        if (std::isdigit(static_cast<unsigned char>(s[i]))) {
            std::size_t j = i;
            while (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j]))) ++j;
            out.push_back(s.substr(i, j - i));
            i = j;
        } else {
            ++i;
        }
    }
    return out;
}

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
        << "  --node-genes <path>    call node_genes.tsv (call --gtf): adds a `gene` column by node\n"
        << "  --min-maf <X>          drop features whose minor (non-modal) genotype frequency < X (default 0.01)\n"
        << "  --model <m>            auto|linear|logistic|lmm (default auto: binary->logistic, else linear)\n"
        << "                         lmm = linear mixed model (EMMAX); needs --kinship or --make-kinship.\n"
        << "  --kinship <path>       precomputed n x n GRM (rows/cols in --samples order) for --model lmm / --pca\n"
        << "  --make-kinship         build the GRM from the genotype matrix itself (region-proximal; see docs)\n"
        << "  --pca <N>              add the top-N kinship PCs as covariates to the GLM (cheap structure control)\n"
        << "  -q, --quiet            less logging\n\n"
        << "Outputs: <prefix>.assoc.tsv (per-feature beta/se/p/p_bonf/q_bh[/gene]) and <prefix>.summary.tsv\n"
        << "  (n_tests, Bonferroni threshold 0.05/n_tests, significant counts). Plot with scripts/plot_associate.R.\n";
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
        else if (a == "--node-genes") opt.node_genes = need(i);
        else if (a == "--phenotype") opt.phenotype = need(i);
        else if (a == "-o" || a == "--out-prefix") opt.out_prefix = need(i);
        else if (a == "--min-maf") opt.min_maf = std::stod(need(i));
        else if (a == "--model") opt.model = need(i);
        else if (a == "--kinship") opt.kinship = need(i);
        else if (a == "--make-kinship") opt.make_kinship = true;
        else if (a == "--pca") opt.pca = std::stoi(need(i));
        else if (a == "-q" || a == "--quiet") opt.quiet = true;
        else throw std::runtime_error("Unknown option for associate: " + a);
    }
    if (opt.genotypes.empty() || opt.samples.empty() || opt.phenotype.empty())
        throw std::runtime_error("associate requires --genotypes, --samples, --phenotype");
    if (opt.model != "auto" && opt.model != "linear" && opt.model != "logistic" && opt.model != "lmm")
        throw std::runtime_error("--model must be auto|linear|logistic|lmm");
    const bool need_kinship = (opt.model == "lmm") || (opt.pca > 0);
    if (need_kinship && opt.kinship.empty() && !opt.make_kinship)
        throw std::runtime_error("--model lmm / --pca require --kinship <file> or --make-kinship");

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
    if (model == "lmm" && binary)
        throw std::runtime_error("--model lmm needs a quantitative phenotype (binary->logistic-LMM not yet implemented)");

    // ---- kinship (for --model lmm and/or --pca): build over the used samples ----
    // K is a GRM in the used-sample order. Either read from --kinship (rows/cols in --samples order,
    // subset to used) or built from the genotype matrix (--make-kinship: standardize each feature,
    // K = Z Z^T / m; note this is region-proximal -- see docs). PCs (if --pca) are appended as fixed
    // GLM covariates; the LMM uses the full K as a random effect.
    Eigen::MatrixXd K;
    if (need_kinship) {
        std::vector<std::size_t> used_cols;
        for (std::size_t c = 0; c < geno_samples.size(); ++c) if (col_to_used[c] >= 0) used_cols.push_back(c);
        K = Eigen::MatrixXd::Zero(n_used, n_used);
        if (!opt.kinship.empty()) {
            // dense n_geno x n_geno matrix (whitespace/comma), subset to used rows+cols
            std::vector<std::vector<double>> M;
            GzLineReader r(opt.kinship);
            if (!r.ok()) throw std::runtime_error("cannot open --kinship: " + opt.kinship);
            std::string line;
            while (r.getline(line)) {
                if (trim(line).empty()) continue;
                std::vector<std::string> f = split(line, '\t');
                if (f.size() < 2) f = split(line, ',');
                if (f.size() < 2) f = split(line, ' ');
                std::vector<double> row; row.reserve(f.size());
                for (const std::string& t : f) { std::string u = trim(t); if (!u.empty()) row.push_back(parse_num(u)); }
                if (!row.empty()) M.push_back(std::move(row));
            }
            if (M.size() < geno_samples.size())
                throw std::runtime_error("--kinship matrix has fewer rows than samples");
            for (std::size_t a = 0; a < n_used; ++a)
                for (std::size_t b = 0; b < n_used; ++b)
                    K(a, b) = M[used_cols[a]][used_cols[b]];
        } else {  // --make-kinship: one pass over the genotype matrix
            GzLineReader gr(opt.genotypes);
            if (!gr.ok()) throw std::runtime_error("cannot open --genotypes: " + opt.genotypes);
            std::string line; std::size_t m_feat = 0;
            while (gr.getline(line)) {
                if (line.empty()) continue;
                std::vector<std::string> f = split(line, ',');
                if (f.size() < 3 + geno_samples.size()) continue;
                Eigen::VectorXd z(n_used);
                double sum = 0.0; std::size_t nf = 0;
                for (std::size_t a = 0; a < n_used; ++a) {
                    const double v = parse_num(f[3 + used_cols[a]]);
                    z(a) = v; if (std::isfinite(v)) { sum += v; ++nf; }
                }
                if (nf == 0) continue;
                const double mean = sum / static_cast<double>(nf);
                for (std::size_t a = 0; a < n_used; ++a) if (!std::isfinite(z(a))) z(a) = mean;  // mean-impute
                z.array() -= mean;
                const double sd = std::sqrt(z.squaredNorm() / static_cast<double>(n_used));
                if (!(sd > 0.0)) continue;
                z /= sd;
                K.selfadjointView<Eigen::Lower>().rankUpdate(z);  // K += z z^T
                ++m_feat;
            }
            if (m_feat == 0) throw std::runtime_error("--make-kinship: no usable features to build the GRM");
            K = K.selfadjointView<Eigen::Lower>();
            K /= static_cast<double>(m_feat);
        }
    }

    // PCA: append the top-N kinship eigenvectors as fixed covariates (cheap structure control).
    if (opt.pca > 0) {
        const Eigen::MatrixXd V = top_eigvecs(K, opt.pca);
        for (Eigen::Index j = 0; j < V.cols(); ++j) covar_names.push_back("kPC" + std::to_string(j + 1));
        for (std::size_t u = 0; u < n_used; ++u)
            for (Eigen::Index j = 0; j < V.cols(); ++j) Z[u].push_back(V(static_cast<Eigen::Index>(u), j));
    }
    const std::size_t ncov_eff = covar_names.size();  // base covariates (+ PCs)

    // LMM null: eigendecompose K once, rotate phenotype + covariates, estimate the variance ratio.
    LmmNull lmm;
    if (model == "lmm") {
        if (!opt.quiet) std::cerr << "[associate] LMM: eigendecomposing kinship (" << n_used << " samples)...\n";
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K);
        if (es.info() != Eigen::Success) throw std::runtime_error("kinship eigendecomposition failed");
        lmm.d = es.eigenvalues();
        lmm.U = es.eigenvectors();
        Eigen::MatrixXd Xcov(n_used, 1 + ncov_eff);
        Eigen::VectorXd yv(n_used);
        for (std::size_t u = 0; u < n_used; ++u) {
            Xcov(static_cast<Eigen::Index>(u), 0) = 1.0;
            for (std::size_t j = 0; j < ncov_eff; ++j) Xcov(static_cast<Eigen::Index>(u), 1 + j) = Z[u][j];
            yv(static_cast<Eigen::Index>(u)) = y[u];
        }
        lmm.UtX = lmm.U.transpose() * Xcov;
        lmm.Uty = lmm.U.transpose() * yv;
        lmm.delta = estimate_delta(lmm);
        if (!opt.quiet) std::cerr << "[associate] LMM: variance ratio delta=sigma_e^2/sigma_g^2 = " << lmm.delta << "\n";
    }

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

    // ---- optional node -> gene(s) map (call --gtf node_genes.tsv) ----
    std::unordered_map<std::string, std::vector<std::string>> node_to_genes;
    if (!opt.node_genes.empty()) {
        GzLineReader r(opt.node_genes);
        if (!r.ok()) throw std::runtime_error("cannot open --node-genes: " + opt.node_genes);
        std::string line; bool header = true;
        while (r.getline(line)) {
            if (line.empty()) continue;
            if (header) { header = false; continue; }
            std::vector<std::string> f = split(line, '\t');
            if (f.size() < 2) continue;
            node_to_genes[trim(f[0])] = split(f[1], ';');  // genes are ';'-separated
        }
    }
    // Join the node->gene map onto a feature's node provenance string -> sorted, ';'-joined gene set.
    auto genes_for = [&](const std::string& nodes) -> std::string {
        if (node_to_genes.empty() || nodes == "." || nodes.empty()) return ".";
        std::set<std::string> g;
        for (const std::string& nid : node_ids_in(nodes)) {
            const auto it = node_to_genes.find(nid);
            if (it == node_to_genes.end()) continue;
            for (const std::string& gene : it->second) { std::string t = trim(gene); if (!t.empty()) g.insert(t); }
        }
        if (g.empty()) return ".";
        std::string out;
        for (const std::string& gene : g) { if (!out.empty()) out += ';'; out += gene; }
        return out;
    };

    // ---- per-feature association over the BIMBAM rows ----
    struct Row {
        std::string id, layer, bubbles, nodes, gene = ".";
        std::size_t n = 0;
        double minor_freq = kNaN, beta = kNaN, se = kNaN, z = kNaN, p = kNaN;
    };
    std::vector<Row> rows;
    std::size_t n_geno_rows = 0, n_dropped_maf = 0, n_dropped_fit = 0;
    const std::size_t p_dim = 2 + ncov_eff;  // intercept + genotype + covariates(+PCs)
    const bool is_lmm = (model == "lmm");

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

        // collect genotype dosage for used samples. GLM/logistic drop NA samples per feature; LMM keeps
        // every used sample (the rotation is over a fixed set) and mean-imputes NA.
        std::vector<double> g; g.reserve(n_used);
        std::vector<double> yy; yy.reserve(n_used);
        std::vector<std::vector<double>> zz; zz.reserve(n_used);
        Eigen::VectorXd glmm;            // full-length imputed genotype (LMM only)
        if (is_lmm) {
            glmm.resize(n_used);
            double sum = 0.0; std::size_t nf = 0;
            for (std::size_t c = 0; c < geno_samples.size(); ++c) {
                if (col_to_used[c] < 0) continue;
                const std::size_t u = static_cast<std::size_t>(col_to_used[c]);
                const double dose = parse_num(f[3 + c]);
                glmm(static_cast<Eigen::Index>(u)) = dose;
                if (std::isfinite(dose)) { sum += dose; ++nf; g.push_back(dose); }
            }
            if (nf == 0) { ++n_dropped_fit; continue; }
            const double mean = sum / static_cast<double>(nf);
            for (std::size_t u = 0; u < n_used; ++u)
                if (!std::isfinite(glmm(static_cast<Eigen::Index>(u)))) glmm(static_cast<Eigen::Index>(u)) = mean;
        } else {
            for (std::size_t c = 0; c < geno_samples.size(); ++c) {
                if (col_to_used[c] < 0) continue;
                const double dose = parse_num(f[3 + c]);
                if (!std::isfinite(dose)) continue;  // genuinely missing (NA = non-traversing)
                const std::size_t u = static_cast<std::size_t>(col_to_used[c]);
                g.push_back(dose);
                yy.push_back(y[u]);
                zz.push_back(Z[u]);
            }
        }
        const std::size_t n = is_lmm ? n_used : g.size();
        if (n < p_dim + 1) { ++n_dropped_fit; continue; }

        // minor (non-modal) genotype frequency: 1 - freq(most common rounded dosage). Works for both
        // presence/absence (0/1) and CN dosage (drops invariant / cohort-rare features). Computed on the
        // observed (non-imputed) dosages collected in g.
        std::unordered_map<long long, std::size_t> cat;
        for (double v : g) ++cat[std::llround(v)];
        std::size_t mode = 0;
        for (const auto& kv : cat) mode = std::max(mode, kv.second);
        const double minor_freq = 1.0 - static_cast<double>(mode) / static_cast<double>(g.size());
        if (minor_freq < opt.min_maf) { ++n_dropped_maf; continue; }

        FitResult fr;
        if (is_lmm) {
            fr = lmm_test(lmm, glmm);
        } else {
            // design matrix: intercept, genotype, covariates(+PCs)
            std::vector<double> X(n * p_dim);
            for (std::size_t i = 0; i < n; ++i) {
                X[i * p_dim + 0] = 1.0;
                X[i * p_dim + 1] = g[i];
                for (std::size_t j = 0; j < ncov_eff; ++j) X[i * p_dim + 2 + j] = zz[i][j];
            }
            fr = (model == "logistic") ? fit_logistic(X, yy, n, p_dim) : fit_linear(X, yy, n, p_dim);
        }
        if (!fr.ok) { ++n_dropped_fit; continue; }

        Row row;
        row.id = id;
        if (auto it = annot.find(id); it != annot.end()) {
            row.layer = it->second.layer; row.bubbles = it->second.bubbles; row.nodes = it->second.nodes;
        } else { row.layer = "."; row.bubbles = "."; row.nodes = "."; }
        row.gene = genes_for(row.nodes);
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
    std::vector<double> zs(n_tests);
    for (std::size_t i = 0; i < n_tests; ++i) zs[i] = rows[i].z;
    const double lambda_gc = genomic_lambda(zs);

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
            << "\tse\tz\tp\tp_bonf\tq_bh\tgene\n";
        for (std::size_t k = 0; k < n_tests; ++k) {
            const std::size_t i = order[k];
            const Row& r = rows[i];
            const double p_bonf = std::isfinite(r.p) ? std::min(1.0, r.p * static_cast<double>(n_tests)) : kNaN;
            out << r.id << '\t' << r.layer << '\t' << r.bubbles << '\t' << r.nodes << '\t'
                << r.n << '\t' << fmt(r.minor_freq) << '\t' << fmt(r.beta) << '\t' << fmt(r.se)
                << '\t' << fmt(r.z) << '\t' << fmt(r.p) << '\t' << fmt(p_bonf) << '\t' << fmt(qv[i])
                << '\t' << r.gene << '\n';
        }
    }
    // ---- summary (also the per-plot threshold lines) ----
    {
        std::ofstream out(opt.out_prefix + ".summary.tsv");
        out << "key\tvalue\n"
            << "model\t" << model << '\n'
            << "phenotype_type\t" << (binary ? "binary" : "quantitative") << '\n'
            << "covariates\t" << ncov_eff << '\n'
            << "pca_covariates\t" << opt.pca << '\n'
            << "samples_used\t" << n_used << '\n'
            << "genotype_rows\t" << n_geno_rows << '\n'
            << "features_tested\t" << n_tests << '\n'
            << "dropped_min_maf\t" << n_dropped_maf << '\n'
            << "dropped_fit\t" << n_dropped_fit << '\n'
            << "bonferroni_threshold\t" << fmt(bonf_threshold) << '\n'
            << "nominal_threshold\t0.05\n"
            << "significant_bonferroni\t" << n_sig_bonf << '\n'
            << "significant_fdr05\t" << n_sig_fdr << '\n'
            << "lambda_gc\t" << fmt(lambda_gc) << '\n';
        if (model == "lmm") out << "lmm_delta\t" << fmt(lmm.delta) << '\n';
    }

    if (!opt.quiet) {
        std::cerr << "[associate] model=" << model << " (" << (binary ? "binary" : "quantitative")
                  << "), samples=" << n_used << ", covariates=" << ncov_eff
                  << (opt.pca > 0 ? " (incl. " + std::to_string(opt.pca) + " PCs)" : "") << "\n"
                  << "[associate] features: tested=" << n_tests << " dropped(maf)=" << n_dropped_maf
                  << " dropped(fit)=" << n_dropped_fit << "\n"
                  << "[associate] Bonferroni threshold (0.05/" << n_tests << ") = " << fmt(bonf_threshold)
                  << "; significant: Bonferroni=" << n_sig_bonf << " FDR<0.05=" << n_sig_fdr << "\n"
                  << "[associate] genomic inflation lambda_GC = " << fmt(lambda_gc) << "\n"
                  << "[associate] wrote " << opt.out_prefix << ".assoc.tsv + .summary.tsv\n";
    }
    return 0;
}

}  // namespace panvar
