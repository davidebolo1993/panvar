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
#include <unordered_set>
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

// Squared Pearson correlation between two equal-length dosage vectors: LD r^2 between variants (clumping)
// and the collinearity guard for feature-tier conditioning. A (near-)constant vector returns 0.
double r2_vec(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    const double ma = a.mean(), mb = b.mean();
    const Eigen::VectorXd da = a.array() - ma, db = b.array() - mb;
    const double na = da.norm(), nb = db.norm();
    if (na < 1e-12 || nb < 1e-12) return 0.0;
    const double r = da.dot(db) / (na * nb);
    return r * r;
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
    std::string kinship;         // precomputed (genome-wide) GRM file (--kinship)
    int pca = 0;                 // top-N kinship PCs as GLM covariates (--pca)
    double min_maf = 0.01;
    std::string model = "auto";  // auto|linear|logistic|lmm
    std::string unit = "auto";   // auto|variant|feature -- the multiple-testing unit (see below)
    double ld_r2 = 0.8;          // variant tier: genotype r^2 above which a variant is an LD shadow of a lead
    long min_ac = 3;             // variant tier: minor-allele-count floor below which a call is flagged low_af
    double cojo_p = -1.0;        // variant tier: forward-stepwise (COJO) entry p; <0 -> use 0.05/Meff
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
        << "  --unit <u>             multiple-testing unit: auto|variant|feature (default auto). variant = one\n"
        << "                         test per SV call (from describe --variant-vcf): honest n_tests + LD-clumping.\n"
        << "                         feature = k-mer/node/edge tests with an effective-tests (Meff) Bonferroni.\n"
        << "                         auto = variant when the feature_annot layer is `variant`, else feature.\n"
        << "  --ld-r2 <X>            variant tier: r^2 above which a variant is an LD shadow of a lead (default 0.8)\n"
        << "  --min-ac <N>           variant tier: flag low_af when the minor-allele count < N (default 3)\n"
        << "  --cojo-p <X>           variant tier: forward-stepwise conditional (COJO) entry p (default 0.05/Meff).\n"
        << "                         selects independent signals; adds p_conditional + cond_role columns.\n"
        << "  --model <m>            auto|linear|logistic|lmm (default auto: binary->logistic, else linear)\n"
        << "                         lmm = linear mixed model (EMMAX); needs an external --kinship GRM.\n"
        << "  --kinship <path>       precomputed (genome-wide) n x n GRM (rows/cols in --samples order)\n"
        << "                         for --model lmm / --pca. panvar is local, so it does not build a GRM itself;\n"
        << "                         supply a genome-wide one, or use PC columns in the phenotype table.\n"
        << "  --pca <N>              add the top-N kinship PCs as covariates to the GLM (needs --kinship)\n"
        << "  -q, --quiet            less logging\n\n"
        << "Outputs: <prefix>.assoc.tsv (per-feature beta/se/p/p_bonf/p_bonf_meff/q_bh, af/an/low_af,\n"
        << "  clump/is_lead[/gene]) and <prefix>.summary.tsv (n_tests, meff, Bonferroni thresholds, significant\n"
        << "  counts, lambda_gc). Plot with scripts/plot_associate.R.\n";
}

struct FeatureAnnot {
    std::string layer = ".";
    std::string bubbles = ".";
    std::string nodes = ".";
    // Extra columns present only in the variant sidecar (feature_annot.variant.tsv.gz):
    // feature_id, layer, encoding, bubbles, nodes, svtype, gene, AF, AN
    std::string svtype = ".";
    std::string gene = ".";
    std::string af = ".";
    std::string an = ".";
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
        else if (a == "--unit") opt.unit = need(i);
        else if (a == "--ld-r2") opt.ld_r2 = std::stod(need(i));
        else if (a == "--min-ac") opt.min_ac = std::stol(need(i));
        else if (a == "--cojo-p") opt.cojo_p = std::stod(need(i));
        else if (a == "--model") opt.model = need(i);
        else if (a == "--kinship") opt.kinship = need(i);
        else if (a == "--pca") opt.pca = std::stoi(need(i));
        else if (a == "-q" || a == "--quiet") opt.quiet = true;
        else throw std::runtime_error("Unknown option for associate: " + a);
    }
    if (opt.genotypes.empty() || opt.samples.empty() || opt.phenotype.empty())
        throw std::runtime_error("associate requires --genotypes, --samples, --phenotype");
    if (opt.model != "auto" && opt.model != "linear" && opt.model != "logistic" && opt.model != "lmm")
        throw std::runtime_error("--model must be auto|linear|logistic|lmm");
    if (opt.unit != "auto" && opt.unit != "variant" && opt.unit != "feature")
        throw std::runtime_error("--unit must be auto|variant|feature");
    const bool need_kinship = (opt.model == "lmm") || (opt.pca > 0);
    if (need_kinship && opt.kinship.empty())
        throw std::runtime_error("--model lmm / --pca require an external --kinship <file> "
                                 "(panvar is local and does not build a GRM; or use PC covariates)");

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

    // ---- kinship (for --model lmm and/or --pca): read an external GRM over the used samples ----
    // K is the GRM in used-sample order, read from --kinship (rows/cols in --samples order, subset to
    // used). panvar is local, so it does not build a GRM from its own region genotypes (that would be
    // proximally contaminated); supply a genome-wide GRM. PCs (--pca) are appended as fixed GLM
    // covariates; the LMM uses the full K as a random effect.
    Eigen::MatrixXd K;
    if (need_kinship) {
        std::vector<std::size_t> used_cols;
        for (std::size_t c = 0; c < geno_samples.size(); ++c) if (col_to_used[c] >= 0) used_cols.push_back(c);
        K = Eigen::MatrixXd::Zero(n_used, n_used);
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
            if (f.size() >= 9) {  // variant sidecar carries svtype/gene/AF/AN
                a.svtype = f[5]; a.gene = f[6]; a.af = f[7]; a.an = f[8];
            }
            annot[f[0]] = std::move(a);
        }
    }

    // ---- multiple-testing unit: variant (one test per SV call) vs feature (k-mer/node/edge) ----
    // Variant mode tests the SV calls directly, so the tests are weakly correlated and n_tests is an
    // honest Bonferroni denominator (refined by LD-clumping). Feature mode keeps the fine-grained tests
    // but corrects with an effective-tests count Meff (the features within one variant are correlated).
    bool variant_mode;
    if (opt.unit == "variant") variant_mode = true;
    else if (opt.unit == "feature") variant_mode = false;
    else {  // auto: variant when the feature_annot is the variant sidecar (its layer is `variant`)
        std::size_t n_var = 0;
        for (const auto& kv : annot) if (kv.second.layer == "variant") ++n_var;
        variant_mode = (!annot.empty() && n_var * 2 >= annot.size());
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
        std::string af = ".", an = ".", svtype = ".";
        int low_af = -1;          // 1/0 in variant mode (AF/AN known), -1 = unknown -> "." in output
        int clump = -1;           // LD-clump id (variant mode), -1 -> "."
        int is_lead = -1;         // 1/0 lead-of-clump (variant mode), -1 -> "."
        std::size_t n = 0;
        double minor_freq = kNaN, beta = kNaN, se = kNaN, z = kNaN, p = kNaN;
    };
    std::vector<Row> rows;
    // Variant tier only: full-length mean-imputed dosage per retained row, for LD r^2 between variants.
    // Kept only in variant mode (few variants), so the k-mer substrate never holds a giant matrix.
    std::vector<Eigen::VectorXd> var_dose;
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
        std::string ann_gene;
        if (auto it = annot.find(id); it != annot.end()) {
            row.layer = it->second.layer; row.bubbles = it->second.bubbles; row.nodes = it->second.nodes;
            row.svtype = it->second.svtype; row.af = it->second.af; row.an = it->second.an;
            ann_gene = it->second.gene;
        } else { row.layer = "."; row.bubbles = "."; row.nodes = "."; }
        // gene: the variant sidecar's GENES if present, else the node->gene join.
        row.gene = (ann_gene != "." && !ann_gene.empty()) ? ann_gene : genes_for(row.nodes);
        // low_af: too few minority observations -> underpowered / asymptotically unstable p. Use the
        // OBSERVED non-modal genotype count (minor_freq * n), not the VCF carrier-AF: for a DUP the
        // dosage is the continuous CN gradient, so a carrier-AF near 1 still carries ample variance and
        // must not be flagged. This metric is uniform for binary (0/1) and copy-number genotypes.
        if (variant_mode)
            row.low_af = (minor_freq * static_cast<double>(n) < static_cast<double>(opt.min_ac)) ? 1 : 0;
        row.n = n; row.minor_freq = minor_freq;
        row.beta = fr.beta; row.se = fr.se; row.z = fr.z; row.p = fr.p;
        rows.push_back(std::move(row));

        // Retain a full-length mean-imputed dosage for variant-tier LD-clumping (cheap: few variants).
        if (variant_mode) {
            Eigen::VectorXd vd(n_used);
            if (is_lmm) {
                vd = glmm;  // already full-length, mean-imputed
            } else {
                for (std::size_t u = 0; u < n_used; ++u) vd(static_cast<Eigen::Index>(u)) =
                    std::numeric_limits<double>::quiet_NaN();
                double sum = 0.0; std::size_t nf = 0;
                for (std::size_t c = 0; c < geno_samples.size(); ++c) {
                    if (col_to_used[c] < 0) continue;
                    const double dose = parse_num(f[3 + c]);
                    if (!std::isfinite(dose)) continue;
                    vd(static_cast<Eigen::Index>(col_to_used[c])) = dose; sum += dose; ++nf;
                }
                const double mean = nf ? sum / static_cast<double>(nf) : 0.0;
                for (std::size_t u = 0; u < n_used; ++u)
                    if (!std::isfinite(vd(static_cast<Eigen::Index>(u)))) vd(static_cast<Eigen::Index>(u)) = mean;
            }
            var_dose.push_back(std::move(vd));
        }
    }

    // ---- multiple testing over the ACTUALLY TESTED features (not genome-wide) ----
    const std::size_t n_tests = rows.size();
    std::vector<double> pv(n_tests);
    for (std::size_t i = 0; i < n_tests; ++i) pv[i] = rows[i].p;
    const std::vector<double> qv = bh_qvalues(pv);

    // Effective number of independent tests (Meff). Raw n_tests over-counts because many tests are
    // correlated -- features within one variant, or variants in LD. Variant mode: greedy LD-clumping
    // (a lead variant, lowest p, claims neighbours with genotype r^2 > --ld-r2; Meff = #leads). Feature
    // mode: Meff = number of distinct bubbles the features map to (the correlated block). Bonferroni then
    // uses Meff; BH-FDR stays the primary control. See docs/algorithms/associate.md.
    std::size_t meff = n_tests;
    if (variant_mode && n_tests > 0) {
        std::vector<std::size_t> byp(n_tests);
        std::iota(byp.begin(), byp.end(), 0);
        std::sort(byp.begin(), byp.end(), [&](std::size_t a, std::size_t b) {
            if (!std::isfinite(rows[a].p)) return false;
            if (!std::isfinite(rows[b].p)) return true;
            return rows[a].p < rows[b].p;
        });
        int next_clump = 0;
        for (std::size_t oi = 0; oi < n_tests; ++oi) {
            const std::size_t i = byp[oi];
            if (rows[i].clump != -1 || !std::isfinite(rows[i].p)) continue;  // shadow or unfittable
            // AF/MAC floor: a low-AF variant has an unstable r^2 and a fragile p, so it must not ANCHOR a
            // clump (and so cannot inflate Meff). It can still be claimed as a shadow by a genuine lead.
            if (rows[i].low_af == 1) continue;
            rows[i].clump = next_clump; rows[i].is_lead = 1;
            for (std::size_t oj = oi + 1; oj < n_tests; ++oj) {
                const std::size_t j = byp[oj];
                if (rows[j].clump != -1 || !std::isfinite(rows[j].p)) continue;
                if (r2_vec(var_dose[i], var_dose[j]) > opt.ld_r2) { rows[j].clump = next_clump; rows[j].is_lead = 0; }
            }
            ++next_clump;
        }
        meff = static_cast<std::size_t>(std::max(1, next_clump));
    } else if (n_tests > 0) {
        std::unordered_set<std::string> blocks;
        for (const Row& r : rows)
            for (const std::string& b : split(r.bubbles, ';')) { const std::string t = trim(b); if (!t.empty() && t != ".") blocks.insert(t); }
        if (!blocks.empty()) meff = blocks.size();
    }

    // --- variant-tier conditional / joint analysis (forward-stepwise, COJO-style) ---
    // Marginal r^2-clumping cannot establish independence: a weak genotypic tag of an extremely strong
    // locus stays genome-wide significant even at r^2 well below --ld-r2. We instead select a set of
    // jointly-independent signals by forward selection -- repeatedly add the variant with the smallest
    // p conditioned on the already-selected set, until none clears the entry threshold (--cojo-p, default
    // 0.05/Meff). Then every variant is reported conditioned on the selected set MINUS itself: a shadow
    // collapses (large p_conditional), a true independent signal survives and is flagged cond_role=signal.
    // Variant tier only (one test per event -> no within-event collinearity); linear/logistic (LMM needs
    // rotation, not done here). Each fit re-uses fit_linear/fit_logistic with the conditioning dosages as
    // extra covariates; the genotype of interest stays at column 1, so FitResult is its conditional Wald.
    std::vector<double> p_cond(rows.size(), kNaN);
    // cond_role: variant tier -> "signal" (COJO-selected) / "shadow"; feature tier -> "lead" /
    // "collinear" (same-event redundancy, not scored) / "conditioned"; "." when not computed.
    std::vector<std::string> cond_role(rows.size(), ".");
    std::size_t cojo_n_signals = 0;
    if (variant_mode && !is_lmm && rows.size() > 1) {
        // conditional p of variant i given a set of conditioning dosage vectors (empty -> marginal-style).
        auto cond_p = [&](std::size_t i, const std::vector<std::size_t>& cond) -> double {
            const std::size_t pc = p_dim + cond.size();        // intercept, g_i, covariates, |cond|
            std::vector<double> X(n_used * pc);
            for (std::size_t u = 0; u < n_used; ++u) {
                std::size_t col = 0;
                X[u * pc + col++] = 1.0;
                X[u * pc + col++] = var_dose[i](static_cast<Eigen::Index>(u));   // target -> column 1
                for (std::size_t j = 0; j < ncov_eff; ++j) X[u * pc + col++] = Z[u][j];
                for (std::size_t c : cond) X[u * pc + col++] = var_dose[c](static_cast<Eigen::Index>(u));
            }
            const FitResult fc = (model == "logistic") ? fit_logistic(X, y, n_used, pc)
                                                       : fit_linear(X, y, n_used, pc);
            return fc.ok ? fc.p : kNaN;
        };
        const double entry = opt.cojo_p > 0.0 ? opt.cojo_p
                                              : (meff > 0 ? 0.05 / static_cast<double>(meff) : 0.05);
        // forward selection of jointly-independent signals
        std::vector<std::size_t> selected;
        std::vector<char> in_sel(rows.size(), 0);
        for (std::size_t step = 0; step < rows.size(); ++step) {
            std::size_t best = rows.size(); double best_p = std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < rows.size(); ++i) {
                if (in_sel[i] || !std::isfinite(rows[i].p)) continue;
                const double pc_i = cond_p(i, selected);
                if (std::isfinite(pc_i) && pc_i < best_p) { best_p = pc_i; best = i; }
            }
            if (best == rows.size() || !(best_p < entry)) break;   // nothing new clears the bar
            selected.push_back(best); in_sel[best] = 1;
        }
        cojo_n_signals = selected.size();
        // report each variant conditioned on the selected set minus itself
        for (std::size_t i = 0; i < rows.size(); ++i) {
            if (!std::isfinite(rows[i].p)) continue;
            cond_role[i] = in_sel[i] ? "signal" : "shadow";
            std::vector<std::size_t> cond;
            for (std::size_t c : selected) if (c != i) cond.push_back(c);
            if (!cond.empty()) p_cond[i] = cond_p(i, cond);   // empty cond (the sole signal) -> NA
        }
    } else if (!variant_mode && !is_lmm && rows.size() > 1) {
        // --- feature-tier single-lead conditional p with a within-bubble collinearity guard ---
        // k-mer / node / edge dosages inside one bubble are ~collinear (they measure the same event), so
        // conditioning them on the top feature is degenerate -> flag features with r^2 > 0.95 vs the lead
        // ("collinear") instead of scoring them. Cross-bubble features get a real conditional p: a mere tag
        // of the lead's variant collapses. Two bounded streaming passes (no feature x sample matrix kept).
        auto parse_full = [&](const std::vector<std::string>& f, Eigen::VectorXd& vd) -> bool {
            vd = Eigen::VectorXd::Constant(static_cast<Eigen::Index>(n_used), kNaN);
            double sum = 0.0; std::size_t nf = 0;
            for (std::size_t c = 0; c < geno_samples.size(); ++c) {
                if (col_to_used[c] < 0) continue;
                const double dose = parse_num(f[3 + c]);
                if (!std::isfinite(dose)) continue;
                vd(static_cast<Eigen::Index>(col_to_used[c])) = dose; sum += dose; ++nf;
            }
            if (nf == 0) return false;
            const double mean = sum / static_cast<double>(nf);
            for (std::size_t u = 0; u < n_used; ++u)
                if (!std::isfinite(vd(static_cast<Eigen::Index>(u)))) vd(static_cast<Eigen::Index>(u)) = mean;
            return true;
        };
        std::size_t lead_idx = 0; double best_p = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < rows.size(); ++i)
            if (std::isfinite(rows[i].p) && rows[i].p < best_p) { best_p = rows[i].p; lead_idx = i; }
        const std::string lead_id = rows[lead_idx].id;
        std::unordered_map<std::string, std::size_t> id_to_row;
        for (std::size_t i = 0; i < rows.size(); ++i) id_to_row[rows[i].id] = i;
        const std::size_t need_cols = 3 + geno_samples.size();
        Eigen::VectorXd gl;  // pass 1: capture the lead feature's dosage
        { GzLineReader g2(opt.genotypes); std::string line;
          while (g2.getline(line)) {
              if (line.empty()) continue;
              std::vector<std::string> f = split(line, ',');
              if (f.size() < need_cols) continue;
              if (trim(f[0]) == lead_id) { parse_full(f, gl); break; }
          } }
        if (gl.size() == static_cast<Eigen::Index>(n_used)) {
            const std::size_t pc = p_dim + 1;  // intercept, g_i, covariates, g_lead
            GzLineReader g3(opt.genotypes); std::string line;
            while (g3.getline(line)) {
                if (line.empty()) continue;
                std::vector<std::string> f = split(line, ',');
                if (f.size() < need_cols) continue;
                const auto it = id_to_row.find(trim(f[0]));
                if (it == id_to_row.end()) continue;
                const std::size_t i = it->second;
                if (i == lead_idx) { cond_role[i] = "lead"; continue; }
                Eigen::VectorXd vd;
                if (!parse_full(f, vd)) continue;
                if (r2_vec(vd, gl) > 0.95) { cond_role[i] = "collinear"; continue; }  // same-event redundancy
                std::vector<double> X(n_used * pc);
                for (std::size_t u = 0; u < n_used; ++u) {
                    X[u * pc + 0] = 1.0;
                    X[u * pc + 1] = vd(static_cast<Eigen::Index>(u));
                    for (std::size_t j = 0; j < ncov_eff; ++j) X[u * pc + 2 + j] = Z[u][j];
                    X[u * pc + 2 + ncov_eff] = gl(static_cast<Eigen::Index>(u));
                }
                const FitResult fc = (model == "logistic") ? fit_logistic(X, y, n_used, pc)
                                                           : fit_linear(X, y, n_used, pc);
                if (fc.ok) { p_cond[i] = fc.p; cond_role[i] = "conditioned"; }
            }
        }
    }

    const double bonf_threshold = n_tests > 0 ? 0.05 / static_cast<double>(n_tests) : kNaN;
    const double bonf_threshold_meff = meff > 0 ? 0.05 / static_cast<double>(meff) : kNaN;
    std::size_t n_sig_bonf = 0, n_sig_fdr = 0, n_sig_meff = 0;
    for (std::size_t i = 0; i < n_tests; ++i) {
        if (std::isfinite(rows[i].p) && rows[i].p < bonf_threshold) ++n_sig_bonf;
        if (std::isfinite(rows[i].p) && rows[i].p < bonf_threshold_meff) ++n_sig_meff;
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
    auto flag = [](int v) { return v < 0 ? std::string(".") : std::to_string(v); };
    {
        std::ofstream out(opt.out_prefix + ".assoc.tsv");
        if (!out) throw std::runtime_error("cannot write " + opt.out_prefix + ".assoc.tsv");
        // p_bonf = raw 0.05/n_tests scaling; p_bonf_meff = the honest effective-tests scaling (Meff).
        out << "feature_id\tlayer\tbubbles\tnodes\tn\tminor_freq\t" << effect
            << "\tse\tz\tp\tp_bonf\tp_bonf_meff\tq_bh\taf\tan\tlow_af\tclump\tis_lead\tgene"
               "\tp_conditional\tcond_role\n";
        for (std::size_t k = 0; k < n_tests; ++k) {
            const std::size_t i = order[k];
            const Row& r = rows[i];
            const double p_bonf = std::isfinite(r.p) ? std::min(1.0, r.p * static_cast<double>(n_tests)) : kNaN;
            const double p_bonf_meff = std::isfinite(r.p) ? std::min(1.0, r.p * static_cast<double>(meff)) : kNaN;
            out << r.id << '\t' << r.layer << '\t' << r.bubbles << '\t' << r.nodes << '\t'
                << r.n << '\t' << fmt(r.minor_freq) << '\t' << fmt(r.beta) << '\t' << fmt(r.se)
                << '\t' << fmt(r.z) << '\t' << fmt(r.p) << '\t' << fmt(p_bonf) << '\t' << fmt(p_bonf_meff)
                << '\t' << fmt(qv[i]) << '\t' << r.af << '\t' << r.an << '\t' << flag(r.low_af)
                << '\t' << flag(r.clump) << '\t' << flag(r.is_lead) << '\t' << r.gene
                << '\t' << fmt(p_cond[i]) << '\t' << cond_role[i] << '\n';
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
            << "unit\t" << (variant_mode ? "variant" : "feature") << '\n'
            << "meff\t" << meff << '\n'
            << (variant_mode ? "independent_variants\t" : "distinct_bubbles\t") << meff << '\n'
            << "dropped_min_maf\t" << n_dropped_maf << '\n'
            << "dropped_fit\t" << n_dropped_fit << '\n'
            << "bonferroni_threshold\t" << fmt(bonf_threshold) << '\n'
            << "bonferroni_threshold_meff\t" << fmt(bonf_threshold_meff) << '\n'
            << "nominal_threshold\t0.05\n"
            << "significant_bonferroni\t" << n_sig_bonf << '\n'
            << "significant_bonferroni_meff\t" << n_sig_meff << '\n'
            << "significant_fdr05\t" << n_sig_fdr << '\n'
            << "lambda_gc\t" << fmt(lambda_gc) << '\n';
        if (variant_mode) out << "cojo_independent_signals\t" << cojo_n_signals << '\n';
        if (model == "lmm") out << "lmm_delta\t" << fmt(lmm.delta) << '\n';
    }

    if (!opt.quiet) {
        std::cerr << "[associate] model=" << model << " (" << (binary ? "binary" : "quantitative")
                  << "), samples=" << n_used << ", covariates=" << ncov_eff
                  << (opt.pca > 0 ? " (incl. " + std::to_string(opt.pca) + " PCs)" : "") << "\n"
                  << "[associate] unit=" << (variant_mode ? "variant" : "feature")
                  << " tested=" << n_tests << " Meff=" << meff << " dropped(maf)=" << n_dropped_maf
                  << " dropped(fit)=" << n_dropped_fit << "\n"
                  << "[associate] Bonferroni (0.05/Meff=" << meff << ") = " << fmt(bonf_threshold_meff)
                  << "; significant: Bonferroni(Meff)=" << n_sig_meff << " FDR<0.05=" << n_sig_fdr << "\n"
                  << "[associate] genomic inflation lambda_GC = " << fmt(lambda_gc) << "\n";
        if (variant_mode)
            std::cerr << "[associate] COJO (forward-stepwise conditional): "
                      << cojo_n_signals << " jointly-independent signal(s)\n";
        std::cerr << "[associate] wrote " << opt.out_prefix << ".assoc.tsv + .summary.tsv\n";
    }
    return 0;
}

}  // namespace panvar
