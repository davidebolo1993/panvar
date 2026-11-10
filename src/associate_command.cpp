#include "panvar/associate_command.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/cosigt.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>
#include <cstdint>
#include <cctype>

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

// The missing-value tokens: empty, NA, NaN, "." and null. PLINK's -9 is missing only as a phenotype,
// where the phenotype reader treats the number -9 so; a covariate of -9 is a value.
bool is_na_token(const std::string& t) {
    std::string u = trim(t);
    for (char& c : u) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return u.empty() || u == "na" || u == "nan" || u == "." || u == "null";
}

// A numeric cell; NaN for a missing token. Anything else must be a complete, finite number (1.5, not
// 1.5abc). The callers add the file and line.
double parse_num(const std::string& t) {
    if (is_na_token(t)) return kNaN;
    const std::string u = trim(t);
    char* end = nullptr;
    errno = 0;
    const double v = std::strtod(u.c_str(), &end);
    if (end == u.c_str() || *end != '\0' || errno == ERANGE || !std::isfinite(v))
        throw std::runtime_error("'" + u + "' is not a number (missing values are NA, NaN, ., null or empty)");
    return v;
}

// A value as the input would write it (5, not 5.000000).
std::string fmt_value(double v) {
    std::ostringstream o;
    o << v;
    return o.str();
}

// Two-sided normal p from z, as erfc(|z|/sqrt2) so a large |z| keeps its precision; floored at 1e-300
// so -log10(p) stays finite.
double wald_p(double z) {
    if (!std::isfinite(z)) return kNaN;
    const double p = std::erfc(std::fabs(z) / std::sqrt(2.0));
    return p < 1e-300 ? 1e-300 : p;
}

// Regularized incomplete beta I_x(a,b), by the continued fraction of Numerical Recipes 6.4, for the
// Student-t tail below.
double betacf(double a, double b, double x) {
    const int kMaxIt = 300;
    const double kEps = 3e-14, kTiny = 1e-300;
    const double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0, d = 1.0 - qab * x / qap;
    if (std::fabs(d) < kTiny) d = kTiny;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= kMaxIt; ++m) {
        const int m2 = 2 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d; if (std::fabs(d) < kTiny) d = kTiny;
        c = 1.0 + aa / c; if (std::fabs(c) < kTiny) c = kTiny;
        d = 1.0 / d; h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d; if (std::fabs(d) < kTiny) d = kTiny;
        c = 1.0 + aa / c; if (std::fabs(c) < kTiny) c = kTiny;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < kEps) break;
    }
    return h;
}

double betai(double a, double b, double x) {
    if (!(x > 0.0)) return 0.0;
    if (!(x < 1.0)) return 1.0;
    const double bt = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
                               a * std::log(x) + b * std::log1p(-x));
    return (x < (a + 1.0) / (a + b + 2.0)) ? bt * betacf(a, b, x) / a
                                           : 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

// Natural log of the two-sided normal tail, unfloored: rows past the 1e-300 floor are still ranked by
// their real p. Past the double range of erfc, its asymptotic series.
double log_wald_p(double z) {
    if (!std::isfinite(z)) return kNaN;
    const double x = std::fabs(z) / std::sqrt(2.0);
    const double p = std::erfc(x);
    if (p > 1e-290) return std::log(p);
    const double x2 = x * x;
    return -x2 - std::log(x) - 0.5 * std::log(M_PI) + std::log1p(-1.0 / (2.0 * x2) + 3.0 / (4.0 * x2 * x2));
}

// Two-sided Student-t tail on `df` degrees of freedom. Falls back to the normal above 1e6 df, where the
// two differ by under 1% for |t| < 10 (about 2% at t = 20, p ~ 1e-88).
double student_p(double t, double df) {
    if (!std::isfinite(t) || !(df > 0.0)) return kNaN;
    if (df > 1e6) return wald_p(t);
    const double p = betai(0.5 * df, 0.5, df / (df + t * t));
    return p < 1e-300 ? 1e-300 : (p > 1.0 ? 1.0 : p);
}

// Natural log of student_p, unfloored: the incomplete beta's prefactor taken in logs, so a t far past
// the floor keeps its rank, and a t on fewer degrees of freedom ranks after one with the same |t| on more.
double log_student_p(double t, double df) {
    if (!std::isfinite(t) || !(df > 0.0)) return kNaN;
    if (df > 1e6) return log_wald_p(t);
    const double a = 0.5 * df, b = 0.5, x = df / (df + t * t);
    if (!(x > 0.0)) return -std::numeric_limits<double>::infinity();
    if (x < (a + 1.0) / (a + b + 2.0)) {
        const double lbt = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) + a * std::log(x) + b * std::log1p(-x);
        return std::min(0.0, lbt + std::log(betacf(a, b, x)) - std::log(a));
    }
    return std::log(std::min(1.0, betai(a, b, x)));
}

// Solve A x = b and also return A^{-1} (Gauss-Jordan, partial pivot). A is p x p (row-major),
// symmetric positive semi-definite (a cross-product X'WX), small (p <= ~12). Returns false if singular.
//
// A is first scaled to a unit diagonal, D^-1/2 A D^-1/2, so the singularity test is on correlations
// and not on the units of the columns.
bool solve_and_invert(std::vector<double> A, const std::vector<double>& b,
                      std::size_t p, std::vector<double>& x, std::vector<double>& inv) {
    std::vector<double> d(p);
    for (std::size_t i = 0; i < p; ++i) {
        const double a = A[i * p + i];
        if (!(a > 0.0) || !std::isfinite(a)) return false;   // a column of zeros, or worse
        d[i] = std::sqrt(a);
    }
    // augmented [As | I], As = D^-1/2 A D^-1/2
    std::vector<double> M(p * 2 * p, 0.0);
    for (std::size_t i = 0; i < p; ++i) {
        for (std::size_t j = 0; j < p; ++j) M[i * 2 * p + j] = A[i * p + j] / (d[i] * d[j]);
        M[i * 2 * p + p + i] = 1.0;
    }
    for (std::size_t col = 0; col < p; ++col) {
        std::size_t piv = col;
        double best = std::fabs(M[col * 2 * p + col]);
        for (std::size_t r = col + 1; r < p; ++r) {
            const double v = std::fabs(M[r * 2 * p + col]);
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-10) return false;
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
    // A^-1 = D^-1/2 As^-1 D^-1/2
    inv.assign(p * p, 0.0);
    for (std::size_t i = 0; i < p; ++i)
        for (std::size_t j = 0; j < p; ++j) inv[i * p + j] = M[i * 2 * p + p + j] / (d[i] * d[j]);
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
    bool spa = false;      // p came from the saddlepoint rather than the normal tail
    bool exact_tail = false;  // p is the exact boundary probability, not any approximation
    bool carriers_exact = false;  // p enumerates the carriers exactly (score_exact_carriers)
    bool normal_fallback = false; // past the cutoff, no tail correction could be made: normal tail
    double beta = kNaN;  // genotype coefficient (column 1)
    double se = kNaN;
    double z = kNaN;
    double p = kNaN;
    double lp = kNaN;    // ln p, unfloored where the tail allows: the rank of a row
};

// X is n x p row-major (col 0 = intercept, col 1 = genotype, cols 2.. = covariates). Returns the
// Wald test on the genotype coefficient (index 1).
FitResult fit_linear(const std::vector<double>& X, const std::vector<double>& y, std::size_t n, std::size_t p) {
    FitResult r;
    if (n <= p) return r;   // one residual degree of freedom is enough for a t test
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
    // A phenotype constant on these samples has nothing to explain; a fit that leaves residuals is
    // tested, however small they are (the covariates alone reproducing the phenotype is refused up front).
    bool constant = true;
    for (std::size_t i = 1; i < n && constant; ++i) if (y[i] != y[0]) constant = false;
    if (constant || !(rss > 0.0)) return r;
    const double sigma2 = rss / static_cast<double>(n - p);
    const double var_b1 = sigma2 * inv[1 * p + 1];
    if (!(var_b1 > 0.0)) return r;
    r.beta = beta[1];
    r.se = std::sqrt(var_b1);
    r.z = r.beta / r.se;
    // Student-t on n - p degrees of freedom: sigma^2 comes from the same residuals, and with a few dozen
    // samples the normal tail would be too small.
    r.p = student_p(r.z, static_cast<double>(n - p));
    r.lp = log_student_p(r.z, static_cast<double>(n - p));
    r.ok = std::isfinite(r.p);
    return r;
}

// The residual sum of squares of y on the n x p design X (row-major); NaN when X is singular.
double rss_of(const std::vector<double>& X, const std::vector<double>& y, std::size_t n, std::size_t p) {
    std::vector<double> XtX(p * p, 0.0), Xty(p, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t a = 0; a < p; ++a) {
            Xty[a] += X[i * p + a] * y[i];
            for (std::size_t b = 0; b < p; ++b) XtX[a * p + b] += X[i * p + a] * X[i * p + b];
        }
    std::vector<double> beta, inv;
    if (!solve_and_invert(XtX, Xty, p, beta, inv)) return kNaN;
    double rss = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double f = 0.0;
        for (std::size_t a = 0; a < p; ++a) f += X[i * p + a] * beta[a];
        rss += (y[i] - f) * (y[i] - f);
    }
    return rss;
}

// P(F(d1, d2) > f).
double f_upper(double f, double d1, double d2) {
    if (!(f > 0.0)) return 1.0;
    return betai(0.5 * d2, 0.5 * d1, d2 / (d2 + d1 * f));
}

// A dosage's classes as a general (not additive) effect, for a measured phenotype: each class of samples
// at one dosage gets its own mean, a class of fewer than 5 samples joined to its nearest. p_genotypic tests
// the classes against the covariates alone (k - 1 df); p_nonadditive tests them against the additive line
// (k - 2 df), so it flags a threshold, a dominant or a U-shaped dose-response. NaN below 3 classes for
// p_nonadditive, below 2 for p_genotypic.
void genotypic_test(const std::vector<double>& g, const std::vector<double>& yk,
                    const std::vector<std::vector<double>>& zk, double* p_geno, double* p_nonadd) {
    *p_geno = kNaN; *p_nonadd = kNaN;
    const std::size_t n = g.size();
    std::map<double, std::size_t> cnt;
    for (double v : g) ++cnt[v];
    std::vector<std::pair<double, std::size_t>> cls(cnt.begin(), cnt.end());
    // Merge small classes into their nearest neighbour, smallest first.
    for (;;) {
        if (cls.size() < 2) break;
        std::size_t small = cls.size();
        for (std::size_t k = 0; k < cls.size(); ++k)
            if (cls[k].second < 5 && (small == cls.size() || cls[k].second < cls[small].second)) small = k;
        if (small == cls.size()) break;
        const std::size_t to = small == 0 ? 1 : small + 1 == cls.size() ? small - 1
                             : (cls[small].first - cls[small - 1].first <= cls[small + 1].first - cls[small].first ? small - 1 : small + 1);
        cls[to].second += cls[small].second;
        cls.erase(cls.begin() + static_cast<long>(small));
    }
    const std::size_t k = cls.size();
    if (k < 2) return;
    // A class that took in another holds several dosages: the additive line is then no longer inside the
    // class model, so the dosage itself joins the class model to keep the two nested.
    const bool merged = cls.size() < cnt.size();
    // A sample's class: the merged class whose original values cover its dosage (nearest kept value).
    auto class_of = [&](double v) {
        std::size_t best = 0;
        for (std::size_t c = 1; c < k; ++c)
            if (std::fabs(cls[c].first - v) < std::fabs(cls[best].first - v)) best = c;
        return best;
    };
    const std::size_t nc = zk.empty() ? 0 : zk[0].size();
    const std::size_t p0 = 1 + nc, pa = p0 + 1, pg = p0 + (k - 1) + (merged ? 1 : 0);
    if (n <= pg + 1) return;
    std::vector<double> X0(n * p0), Xa(n * pa), Xg(n * pg);
    for (std::size_t i = 0; i < n; ++i) {
        X0[i * p0] = 1.0; Xa[i * pa] = 1.0; Xg[i * pg] = 1.0;
        for (std::size_t j = 0; j < nc; ++j) { X0[i * p0 + 1 + j] = zk[i][j]; Xa[i * pa + 1 + j] = zk[i][j]; Xg[i * pg + 1 + j] = zk[i][j]; }
        Xa[i * pa + p0] = g[i];
        const std::size_t c = class_of(g[i]);
        for (std::size_t d = 1; d < k; ++d) Xg[i * pg + p0 + d - 1] = (c == d) ? 1.0 : 0.0;
        if (merged) Xg[i * pg + p0 + k - 1] = g[i];
    }
    const double r0 = rss_of(X0, yk, n, p0), ra = rss_of(Xa, yk, n, pa), rg = rss_of(Xg, yk, n, pg);
    if (!std::isfinite(r0) || !std::isfinite(rg) || !(rg > 0.0)) return;
    const double d2 = static_cast<double>(n - pg);
    const double dg = static_cast<double>(pg - p0), dn = static_cast<double>(pg - pa);
    // Floored at 1e-300, as every p in the output.
    const auto floored = [](double v) { return std::isfinite(v) ? std::max(v, 1e-300) : v; };
    *p_geno = floored(f_upper(((r0 - rg) / dg) / (rg / d2), dg, d2));
    if (dn >= 1.0 && std::isfinite(ra))
        *p_nonadd = floored(f_upper(((ra - rg) / dn) / (rg / d2), dn, d2));
}

// IRLS for logistic regression on an arbitrary design. Returns the coefficients, the inverse
// information matrix and the fitted probabilities, so callers can build either a Wald or a score test
// from one fit.
bool logistic_irls(const std::vector<double>& X, const std::vector<double>& y, std::size_t n, std::size_t p,
                   std::vector<double>& beta, std::vector<double>& inv, std::vector<double>* mu_out) {
    beta.assign(p, 0.0);
    inv.clear();
    bool converged = false;
    for (int iter = 0; iter < 50; ++iter) {
        std::vector<double> XtWX(p * p, 0.0), XtWz(p, 0.0);
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
        if (!solve_and_invert(XtWX, XtWz, p, nb, inv)) return false;
        bool small = true;
        for (std::size_t a = 0; a < p; ++a) {
            if (std::fabs(nb[a] - beta[a]) > 1e-8 * (1.0 + std::fabs(nb[a]))) small = false;
            beta[a] = nb[a];
        }
        if (small) { converged = true; break; }
    }
    // A fit that hit the iteration cap has not found the maximum (under separation it diverges).
    if (!converged) return false;
    if (mu_out != nullptr) {
        mu_out->assign(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double eta = 0.0;
            for (std::size_t a = 0; a < p; ++a) eta += X[i * p + a] * beta[a];
            (*mu_out)[i] = 1.0 / (1.0 + std::exp(-eta));
        }
    }
    return true;
}

// Firth-penalised logistic regression (Firth 1993; Heinze & Schemper 2002).
//
// Under separation the maximum likelihood estimate diverges. Firth adds Jeffreys' invariant prior,
// |I(beta)|^{1/2}, as a penalty, giving the penalised score
//
//   U*_j = sum_i x_ij [ y_i - mu_i + h_i (0.5 - mu_i) ],    h = diag( W^{1/2} X (X'WX)^-1 X' W^{1/2} )
//
// which stays finite because h_i(0.5 - mu_i) pulls the fitted probabilities back off 0 and 1, and is
// first-order unbiased where ordinary ML is not.
//
// Used only for the reported effect size when ML fails; the p-value stays the score test's.
bool firth_logistic(const std::vector<double>& X, const std::vector<double>& y, std::size_t n,
                    std::size_t p, std::vector<double>& beta, std::vector<double>& inv) {
    beta.assign(p, 0.0);
    inv.clear();
    std::vector<double> mu(n), w(n);
    // The penalised log-likelihood l(beta) + 0.5 log|X'WX| at `b`, filling mu, w and X'WX on the way;
    // -inf when X'WX is not positive definite. log mu and log(1-mu) as softplus, so a fitted
    // probability of 1e-300 does not round to log(0).
    std::vector<double> XtWX(p * p);
    auto penalised = [&](const std::vector<double>& b) {
        std::fill(XtWX.begin(), XtWX.end(), 0.0);
        double ll = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double eta = 0.0;
            for (std::size_t a = 0; a < p; ++a) eta += X[i * p + a] * b[a];
            const auto softplus = [](double t) { return t > 30.0 ? t : std::log1p(std::exp(t)); };
            ll -= y[i] * softplus(-eta) + (1.0 - y[i]) * softplus(eta);
            mu[i] = 1.0 / (1.0 + std::exp(-eta));
            w[i] = std::max(mu[i] * (1.0 - mu[i]), 1e-10);
            for (std::size_t a = 0; a < p; ++a)
                for (std::size_t c = a; c < p; ++c) XtWX[a * p + c] += X[i * p + a] * w[i] * X[i * p + c];
        }
        Eigen::MatrixXd F(p, p);
        for (std::size_t a = 0; a < p; ++a)
            for (std::size_t c = 0; c < p; ++c) F(a, c) = XtWX[std::min(a, c) * p + std::max(a, c)];
        const Eigen::LLT<Eigen::MatrixXd> llt(F);
        if (llt.info() != Eigen::Success) return -std::numeric_limits<double>::infinity();
        const Eigen::ArrayXd diag = llt.matrixL().toDenseMatrix().diagonal().array();
        if ((diag <= 0.0).any()) return -std::numeric_limits<double>::infinity();
        return ll + diag.log().sum();   // 0.5 log|F| = sum log diag(L)
    };
    double cur = penalised(beta);
    if (!std::isfinite(cur)) return false;
    for (int iter = 0; iter < 200; ++iter) {
        // penalised() left mu, w and X'WX at the current beta
        for (std::size_t a = 0; a < p; ++a) for (std::size_t c = 0; c < a; ++c) XtWX[a * p + c] = XtWX[c * p + a];
        std::vector<double> unit(p, 0.0), sol;
        if (!solve_and_invert(XtWX, unit, p, sol, inv)) return false;
        if (inv.size() != p * p) return false;
        // leverages h_i = w_i * x_i' (X'WX)^-1 x_i
        std::vector<double> Ustar(p, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double h = 0.0;
            for (std::size_t a = 0; a < p; ++a) {
                double t = 0.0;
                for (std::size_t b = 0; b < p; ++b) t += inv[a * p + b] * X[i * p + b];
                h += X[i * p + a] * t;
            }
            h *= w[i];
            const double resid = (y[i] - mu[i]) + h * (0.5 - mu[i]);
            for (std::size_t a = 0; a < p; ++a) Ustar[a] += X[i * p + a] * resid;
        }
        std::vector<double> d(p, 0.0);
        double largest = 0.0;
        for (std::size_t a = 0; a < p; ++a) {
            for (std::size_t b = 0; b < p; ++b) d[a] += inv[a * p + b] * Ustar[b];
            largest = std::max(largest, std::fabs(d[a]));
        }
        if (largest < 1e-8) return true;   // inv is the inverse information at the estimate
        // The Newton step, at most 5 in any coordinate, then halved until the penalised likelihood
        // rises (Heinze & Schemper; logistf).
        const double shrink = largest > 5.0 ? 5.0 / largest : 1.0;
        for (double& v : d) v *= shrink;
        std::vector<double> next(p);
        double val = -std::numeric_limits<double>::infinity();
        bool rose = false;
        for (int half = 0; half < 40; ++half) {
            for (std::size_t a = 0; a < p; ++a) next[a] = beta[a] + d[a];
            val = penalised(next);
            // Strictly: two points either side of the maximum have equal value.
            if (std::isfinite(val) && val > cur) { rose = true; break; }
            for (double& v : d) v *= 0.5;
        }
        if (!rose) {
            penalised(beta);   // restore mu, w, X'WX at beta for the final inverse
            for (std::size_t a = 0; a < p; ++a) for (std::size_t c = 0; c < a; ++c) XtWX[a * p + c] = XtWX[c * p + a];
            std::vector<double> sol2;
            return solve_and_invert(XtWX, unit, p, sol2, inv) && largest < 1e-4;
        }
        beta = next;
        cur = val;
    }
    return false;   // did not settle even with the penalty
}

// Rao score test on the genotype column, evaluated under the NULL fit (covariates only, genotype
// coefficient fixed at zero).
//
// Not the Wald test: for a rare variant in an unbalanced study the fit nears separation, the standard
// error grows faster than the estimate, and the Wald statistic collapses. The score test never fits the
// alternative. Each feature gets its own null fit, on its own complete-case samples.

// ---- Saddlepoint approximation for the binary score statistic (Dey et al. fastSPA; SAIGE) --------
//
// The score S = sum_i Gt_i (y_i - mu_i) is a sum of INDEPENDENT bounded terms, so its cumulant
// generating function is available in closed form -- no approximation needed:
//
//   K(t)   = sum_i [ log(1 - mu_i + mu_i e^{t Gt_i}) - t Gt_i mu_i ]
//   K'(t)  = sum_i [ Gt_i mu_i e^{t Gt_i} / (1 - mu_i + mu_i e^{t Gt_i}) - Gt_i mu_i ]
//   K''(t) = sum_i [ Gt_i^2 mu_i (1-mu_i) e^{t Gt_i} / (1 - mu_i + mu_i e^{t Gt_i})^2 ]
//
// The normal approximation matches only the first two cumulants and drifts in the far tail when the
// terms are skewed (a rare variant, unbalanced cases). The saddlepoint expands about the point where the
// tilted distribution is centred on the observed value, and stays accurate there.
struct SpaCgf {
    const std::vector<double>* gt;
    const std::vector<double>* mu;
};

// log(a) + log(1 + e^{b-a}) without overflowing either branch.
double log_add(double a, double b) {
    if (!std::isfinite(a)) return b;
    if (!std::isfinite(b)) return a;
    const double m = std::max(a, b);
    return m + std::log1p(std::exp(std::min(a, b) - m));
}

// log(1 - mu + mu e^x), as a log-sum-exp so a large saddlepoint cannot overflow exp(x).
double log_mix(double mu, double x) {
    if (!(mu > 0.0)) return 0.0;
    if (!(mu < 1.0)) return x;
    return log_add(std::log1p(-mu), std::log(mu) + x);
}

// The exponentially tilted mean of y_i, in the numerically stable logistic form
//   mu e^x / (1 - mu + mu e^x) = 1 / (1 + ((1-mu)/mu) e^{-x})
// which is bounded for either sign of x.
double tilted_mean(double mu, double x) {
    if (!(mu > 0.0)) return 0.0;
    if (!(mu < 1.0)) return 1.0;
    const double z = std::log1p(-mu) - std::log(mu) - x;   // log((1-mu)/mu) - x
    if (z > 40.0) return 0.0;
    if (z < -40.0) return 1.0;
    return 1.0 / (1.0 + std::exp(z));
}

double spa_K(double t, const SpaCgf& c) {
    double k = 0.0;
    for (std::size_t i = 0; i < c.gt->size(); ++i) {
        const double g = (*c.gt)[i], m = (*c.mu)[i];
        if (g == 0.0) continue;
        k += log_mix(m, t * g) - t * g * m;
    }
    return k;
}

// First and second derivatives together: the saddlepoint solve needs both at every step.
void spa_K12(double t, const SpaCgf& c, double& k1, double& k2) {
    k1 = 0.0; k2 = 0.0;
    for (std::size_t i = 0; i < c.gt->size(); ++i) {
        const double g = (*c.gt)[i], m = (*c.mu)[i];
        if (g == 0.0) continue;
        const double q = tilted_mean(m, t * g);
        k1 += g * (q - m);
        k2 += g * g * q * (1.0 - q);
    }
}

// Solve K'(zeta) = s. K' is strictly increasing: bracket, then bisect with a Newton step where safe.
bool spa_root(double s, const SpaCgf& c, double& zeta) {
    double k1 = 0.0, k2 = 0.0;
    spa_K12(0.0, c, k1, k2);
    if (!(k2 > 0.0)) return false;
    double lo = 0.0, hi = 0.0;
    const double step = (s > k1) ? 1.0 : -1.0;
    double t = 0.0;
    for (int it = 0; it < 200; ++it) {
        t += step * (0.05 + 0.5 * std::fabs(t));
        spa_K12(t, c, k1, k2);
        if (!std::isfinite(k1)) return false;
        if ((step > 0.0 && k1 >= s) || (step < 0.0 && k1 <= s)) { lo = std::min(0.0, t); hi = std::max(0.0, t); break; }
        if (it == 199) return false;
    }
    if (lo == hi) return false;
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        spa_K12(mid, c, k1, k2);
        if (!std::isfinite(k1)) return false;
        if (k1 < s) lo = mid; else hi = mid;
        if (hi - lo < 1e-12 * std::max(1.0, std::fabs(mid))) break;
    }
    zeta = 0.5 * (lo + hi);
    // A saddlepoint this far out means the search ran to the support boundary rather than a genuine
    // interior root; the expansion is not valid there.
    return std::isfinite(zeta) && std::fabs(zeta) < 50.0;
}

// The FAR tail at s, by Lugannani-Rice: P(S >= s) when the saddlepoint is positive, P(S <= s) when it
// is negative. Each side is evaluated in the form that is numerically stable there -- the survival
// function above the mean, the distribution function below it -- rather than one of them by subtraction.
double spa_far_tail(double s, const SpaCgf& c) {
    double zeta = 0.0;
    if (!spa_root(s, c, zeta)) return kNaN;
    if (std::fabs(zeta) < 1e-10) return 0.5;
    double k1 = 0.0, k2 = 0.0;
    spa_K12(zeta, c, k1, k2);
    if (!(k2 > 0.0)) return kNaN;
    const double arg = 2.0 * (zeta * s - spa_K(zeta, c));
    if (!(arg > 0.0)) return kNaN;
    const double w = (zeta > 0.0 ? 1.0 : -1.0) * std::sqrt(arg);
    const double v = zeta * std::sqrt(k2);                  // shares the sign of zeta
    if (!(std::fabs(w) > 1e-8) || !(std::fabs(v) > 1e-12)) return kNaN;
    const double pdf = std::exp(-0.5 * w * w) / std::sqrt(2.0 * M_PI);
    const double corr = pdf * (1.0 / w - 1.0 / v);
    // Phi(w) via erfc keeps both tails accurate without cancellation.
    const double t = (zeta > 0.0) ? (0.5 * std::erfc(w / std::sqrt(2.0)) - corr)   // P(S >= s)
                                  : (0.5 * std::erfc(-w / std::sqrt(2.0)) + corr); // P(S <= s)
    // Negative on a lattice-like score (a handful of carriers): the expansion failed.
    if (!std::isfinite(t) || !(t > 0.0)) return kNaN;
    return std::min(1.0, t);
}

// Chernoff bound on a tail: P(S >= s) <= exp(K(t) - t s) for every t > 0 (and the mirror for the lower
// tail, t < 0), minimised over a grid of t. An upper bound, so using it where the saddlepoint cannot be
// formed errs conservative. It serves the tail opposite the observation, whose saddlepoint can lie past
// the edge of the support.
double chernoff_tail(double s, const SpaCgf& c, bool upper) {
    double best = 0.0;   // log bound; t = 0 gives exp(0) = 1
    for (double t = 0.01; t <= 1e4; t *= 1.5) {
        const double tt = upper ? t : -t;
        const double lb = spa_K(tt, c) - tt * s;
        if (std::isfinite(lb)) best = std::min(best, lb);
    }
    return std::exp(best);
}

// Two-sided p for the observed score. The score's distribution is NOT symmetric under case/control
// imbalance, so both tails are computed rather than one doubled.
double spa_two_sided(double s, const std::vector<double>& gt, const std::vector<double>& mu,
                     const std::vector<double>* y = nullptr, bool* exact_out = nullptr) {
    // S is a bounded sum, and each tail is evaluated on its own: a threshold can be outside the support
    // (probability zero), on its boundary (one configuration), or inside it.
    double smax = 0.0, smin = 0.0;
    for (std::size_t i = 0; i < gt.size(); ++i) {
        const double a = gt[i] * (1.0 - mu[i]), b = gt[i] * (0.0 - mu[i]);
        smax += std::max(a, b);
        smin += std::min(a, b);
    }
    const double span = smax - smin;
    if (!(span > 0.0)) return kNaN;

    // P(S = smax) = prod_{gt>0} mu_i * prod_{gt<0} (1 - mu_i): a single Bernoulli configuration, so it
    // is exact. In logs, so it does not underflow.
    auto log_extreme = [&](bool upper) {
        double lg = 0.0;
        for (std::size_t i = 0; i < gt.size(); ++i) {
            if (gt[i] == 0.0) continue;
            const bool want_one = upper ? (gt[i] > 0.0) : (gt[i] < 0.0);
            const double q = want_one ? mu[i] : (1.0 - mu[i]);
            if (!(q > 0.0)) return -std::numeric_limits<double>::infinity();
            lg += std::log(q);
        }
        return lg;
    };
    // Whether the observed outcome is the extreme configuration, read from y rather than from a tolerance
    // on the statistic; it applies only to the side the observation is on.
    auto observed_is_extreme = [&](bool upper) {
        if (y == nullptr || y->size() != gt.size()) return true;
        for (std::size_t i = 0; i < gt.size(); ++i) {
            if (gt[i] == 0.0) continue;
            const double want = (upper ? (gt[i] > 0.0) : (gt[i] < 0.0)) ? 1.0 : 0.0;
            if (std::fabs((*y)[i] - want) > 1e-12) return false;
        }
        return true;
    };

    const double tol = 1e-8 * span;
    const SpaCgf c{&gt, &mu};
    bool up_exact = false, dn_exact = false;

    // Upper tail P(S >= |s|)
    double up = kNaN;
    const double u = std::fabs(s);
    if (u > smax + tol) { up = 0.0; up_exact = true; }                       // outside the support
    else if (u >= smax - tol && (s < 0.0 || observed_is_extreme(true))) {
        const double lu = log_extreme(true);
        up = std::isfinite(lu) ? std::exp(lu) : kNaN;
        up_exact = std::isfinite(up);
    } else {
        up = spa_far_tail(u, c);
        if (!std::isfinite(up) && s < 0.0) up = chernoff_tail(u, c, true);    // the opposite tail
    }

    // Lower tail P(S <= -|s|)
    double dn = kNaN;
    const double l = -std::fabs(s);
    if (l < smin - tol) { dn = 0.0; dn_exact = true; }                       // outside the support
    else if (l <= smin + tol && (s > 0.0 || observed_is_extreme(false))) {
        const double ll = log_extreme(false);
        dn = std::isfinite(ll) ? std::exp(ll) : kNaN;
        dn_exact = std::isfinite(dn);
    } else {
        dn = spa_far_tail(l, c);
        if (!std::isfinite(dn) && s >= 0.0) dn = chernoff_tail(l, c, false);  // the opposite tail
    }

    if (!std::isfinite(up) || !std::isfinite(dn)) return kNaN;
    const double p = up + dn;
    if (!std::isfinite(p) || p <= 0.0) return kNaN;
    // Exact only when both tails were exact or zero; a mix is reported as the saddlepoint.
    if (exact_out != nullptr) *exact_out = (up_exact && dn_exact);
    return std::min(1.0, p);
}

// The two-sided tail counted over the outcomes, for few carriers (samples off the modal genotype), whose
// outcomes form a lattice of at most 2^m points where both the normal tail and the saddlepoint fail. Exact
// without covariates; with covariates the fitted chances are taken as known.
//   - With at most `max_enum` samples in all, the whole score S = sum_i gt_i (y_i - mu_i) is enumerated,
//     the y_i independent Bernoulli(mu_i): P(S >= |s|) + P(S <= -|s|).
//   - Otherwise the test is on the carriers alone: D = sum over carriers of (g_i - modal)(y_i - mu_i), the
//     carriers' outcomes enumerated at their fitted chances, or, without covariates, given the number of
//     cases (each set of s carriers being the cases with chance C(n - m, K - s) / C(n, K), Fisher's
//     hypergeometric, two-sided by distance from the expected count rather than by probability). Adding the non-carriers' many small terms to the score would smear each lattice point,
//     counting the observed one by half (a mid-p): with ten carriers, three of them among a tenth of
//     cases, 0.038 where the exact tail is 0.069.
// NaN when the carriers number more than `max_enum`, or `max_enum` is 0.
double carriers_exact_p(double s, const std::vector<double>& g, const std::vector<double>& gt,
                        const std::vector<double>& mu, std::size_t max_enum, const std::vector<double>& y,
                        bool intercept_only = false) {
    const std::size_t n = g.size();
    if (max_enum == 0 || n == 0) return kNaN;
    std::map<double, std::size_t> count;
    for (double v : g) ++count[v];
    double modal = 0.0; std::size_t best = 0;
    for (const auto& [v, k] : count) if (k > best) { best = k; modal = v; }   // ties: the smaller value
    std::vector<std::size_t> enumerated, rest;
    double span = 0.0, gt_scale = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        span += std::fabs(gt[i]);
        gt_scale = std::max(gt_scale, std::fabs(gt[i]));
        if (n <= max_enum || g[i] != modal) enumerated.push_back(i); else rest.push_back(i);
    }
    if (enumerated.size() > max_enum) return kNaN;
    const double tol = 1e-9 * std::max(span, 1e-300);
    // The enumerated sum, configuration by configuration: (value, probability). Each sample's term: the
    // whole adjusted score when every sample is enumerated, else the carrier's distance from the mode.
    // Without covariates the carriers are always tested given the number of cases, however few the samples.
    if (intercept_only) {
        std::vector<std::size_t> car;
        for (std::size_t i = 0; i < n; ++i) if (g[i] != modal) car.push_back(i);
        enumerated.swap(car);
        double K = 0.0;
        for (double v : y) K += v;
        const double m = static_cast<double>(enumerated.size()), nn = static_cast<double>(n);
        auto lchoose = [](double a, double b) {
            return (b < 0.0 || b > a) ? -std::numeric_limits<double>::infinity()
                                      : std::lgamma(a + 1.0) - std::lgamma(b + 1.0) - std::lgamma(a - b + 1.0);
        };
        const double mu_c = K / nn;
        double obs = 0.0, scale = 0.0;
        for (const std::size_t i : enumerated) { obs += (g[i] - modal) * (y[i] - mu_c); scale += std::fabs(g[i] - modal); }
        // Carriers in classes of equal weight: a class of m_c carriers with k_c cases counts C(m_c, k_c)
        // subsets, so the classes' case counts are enumerated, not the subsets (a 0/1 row is one class).
        std::map<double, int> wcls;
        for (const std::size_t i : enumerated) ++wcls[g[i] - modal];
        std::vector<std::pair<double, int>> cls_v(wcls.begin(), wcls.end());
        // (value, cases among carriers, log multiplicity)
        std::vector<std::tuple<double, int, double>> dist{{0.0, 0, 0.0}};
        for (const auto& [w, mc] : cls_v) {
            std::vector<std::tuple<double, int, double>> next;
            next.reserve(dist.size() * static_cast<std::size_t>(mc + 1));
            for (const auto& [v, k, lm] : dist)
                for (int c = 0; c <= mc; ++c)
                    next.emplace_back(v + w * (static_cast<double>(c) - static_cast<double>(mc) * mu_c), k + c,
                                      lm + lchoose(static_cast<double>(mc), static_cast<double>(c)));
            dist.swap(next);
        }
        const double a = std::fabs(obs), tol_h = 1e-9 * std::max(scale, 1e-300);
        const double lall = lchoose(nn, K);
        double p = 0.0;
        for (const auto& [v, k, lm] : dist)
            if (v >= a - tol_h || v <= -a + tol_h) p += std::exp(lm + lchoose(nn - m, K - k) - lall);
        if (!(p > 0.0) || !std::isfinite(p)) return kNaN;
        return std::min(1.0, p);
    }
    const bool whole = rest.empty();
    double obs = whole ? s : 0.0, scale = 0.0;
    for (const std::size_t i : enumerated) {
        const double w = whole ? gt[i] : g[i] - modal;
        if (!whole) obs += w * (y[i] - mu[i]);
        scale += std::fabs(w);
    }
    const double tol_d = whole ? tol : 1e-9 * std::max(scale, 1e-300);
    // Carriers with whole-number weights (dosage steps): D >= a is a bound on T = sum w_i y_i, whose
    // distribution is a polynomial product over the carriers, O(m * sum w), not 2^m.
    if (!whole) {
        bool integral = true;
        long total = 0;
        for (const std::size_t i : enumerated) {
            const double w = g[i] - modal;
            // A weight that rounds to 0 is no step of the lattice, however close to 0 it is.
            if (std::fabs(w) > 1000.0 || std::fabs(w - std::round(w)) > 1e-9 || std::lround(w) == 0) { integral = false; break; }
            total += std::labs(std::lround(w));
        }
        if (integral && total < 100000) {
            // T counts from the most negative possible sum; offset by the negative weights.
            long neg = 0;
            for (const std::size_t i : enumerated) { const long w = std::lround(g[i] - modal); if (w < 0) neg += -w; }
            std::vector<double> pt(static_cast<std::size_t>(total + 1), 0.0);
            pt[static_cast<std::size_t>(neg)] = 1.0;   // T + neg, with no carrier a case
            double mean_shift = 0.0;
            for (const std::size_t i : enumerated) {
                const long w = std::lround(g[i] - modal);
                mean_shift += static_cast<double>(w) * mu[i];
                std::vector<double> nx(pt.size(), 0.0);
                for (std::size_t t = 0; t < pt.size(); ++t) {
                    if (pt[t] == 0.0) continue;
                    nx[t] += pt[t] * (1.0 - mu[i]);
                    const long t2 = static_cast<long>(t) + w;
                    if (t2 >= 0 && t2 < static_cast<long>(nx.size())) nx[static_cast<std::size_t>(t2)] += pt[t] * mu[i];
                }
                pt.swap(nx);
            }
            const double a = std::fabs(obs);
            double p = 0.0;
            for (std::size_t t = 0; t < pt.size(); ++t) {
                const double v = static_cast<double>(static_cast<long>(t) - neg) - mean_shift;
                if (v >= a - tol_d || v <= -a + tol_d) p += pt[t];
            }
            if (!(p > 0.0) || !std::isfinite(p)) return kNaN;
            return std::min(1.0, p);
        }
    }
    std::vector<std::pair<double, double>> dist{{0.0, 1.0}};
    for (const std::size_t i : enumerated) {
        const double w = whole ? gt[i] : g[i] - modal;
        std::vector<std::pair<double, double>> next;
        next.reserve(dist.size() * 2);
        for (const auto& [v, pr] : dist) {
            next.emplace_back(v + w * (1.0 - mu[i]), pr * mu[i]);
            next.emplace_back(v - w * mu[i], pr * (1.0 - mu[i]));
        }
        dist.swap(next);
    }
    const double a = std::fabs(obs);
    double p = 0.0;
    for (const auto& [v, pr] : dist)
        if (v >= a - tol_d || v <= -a + tol_d) p += pr;
    if (!(p > 0.0) || !std::isfinite(p)) return kNaN;
    return std::min(1.0, p);
}

FitResult score_logistic(const std::vector<double>& X, const std::vector<double>& y,
                         std::size_t n, std::size_t p, std::size_t exact_carriers, double spa_cutoff = 2.0) {
    FitResult r;
    if (p < 2 || n <= p + 1) return r;
    const std::size_t p0 = p - 1;
    std::vector<double> X0(n * p0);
    std::vector<double> g(n);
    for (std::size_t i = 0; i < n; ++i) {
        g[i] = X[i * p + 1];
        X0[i * p0 + 0] = X[i * p + 0];                       // intercept
        for (std::size_t a = 2; a < p; ++a) X0[i * p0 + (a - 1)] = X[i * p + a];
    }
    std::vector<double> b0, inv0, mu;
    if (!logistic_irls(X0, y, n, p0, b0, inv0, &mu)) return r;
    if (inv0.size() != p0 * p0) return r;

    double U = 0.0, gWg = 0.0;
    std::vector<double> gWX(p0, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double w = std::max(mu[i] * (1.0 - mu[i]), 1e-12);
        U += g[i] * (y[i] - mu[i]);
        gWg += g[i] * w * g[i];
        for (std::size_t a = 0; a < p0; ++a) gWX[a] += g[i] * w * X0[i * p0 + a];
    }
    double proj = 0.0;
    for (std::size_t a = 0; a < p0; ++a)
        for (std::size_t b = 0; b < p0; ++b) proj += gWX[a] * inv0[a * p0 + b] * gWX[b];
    const double V = gWg - proj;                             // variance of the score, covariate-adjusted
    if (!(V > 0.0) || !std::isfinite(U)) return r;
    const double chi2 = (U * U) / V;
    r.z = (U >= 0.0 ? 1.0 : -1.0) * std::sqrt(chi2);
    r.p = wald_p(r.z);                                       // chi2 on 1 df is the squared normal
    r.lp = log_wald_p(r.z);
    r.ok = std::isfinite(r.p);

    // With few carriers the score is a coarse lattice and the normal tail is wrong at any z (two carriers,
    // one a case: 0.047 where the exact tail is 0.18), so it is counted exactly, or corrected. Otherwise, past
    // |z| = 2 the tail is corrected; below it the normal agrees and is far cheaper (as in SAIGE).
    std::size_t carriers = 0;
    {
        std::map<double, std::size_t> count;
        for (double v : g) ++count[v];
        std::size_t best = 0;
        for (const auto& kv : count) best = std::max(best, kv.second);
        carriers = n - best;
    }
    // With covariates the carriers' count holds the non-carriers fixed, which is conservative once the
    // carriers are a sizeable share of the samples: there the saddlepoint on the adjusted genotype is used,
    // at any z. With n within the limit every sample's outcome is enumerated, always: no lattice is smoothed,
    // though with covariates the fitted chances are taken as known, so the tail is a close approximation.
    const bool small = exact_carriers > 0 && carriers <= exact_carriers;
    const bool few = exact_carriers > 0 &&
                     (n <= exact_carriers ||
                      (small && (p0 == 1 || static_cast<double>(carriers) <= 0.05 * static_cast<double>(n))));
    if (few || small || (spa_cutoff > 0.0 && std::fabs(r.z) > spa_cutoff)) {
        // Covariate-adjusted genotype: Gt = g - X0 (X0'W X0)^-1 X0'W g. The score is Gt'(y - mu), and
        // its terms are independent, which is what makes the exact CGF above available.
        std::vector<double> bproj(p0, 0.0);
        for (std::size_t a = 0; a < p0; ++a)
            for (std::size_t b = 0; b < p0; ++b) bproj[a] += inv0[a * p0 + b] * gWX[b];
        std::vector<double> gt(n);
        for (std::size_t i = 0; i < n; ++i) {
            double fit = 0.0;
            for (std::size_t a = 0; a < p0; ++a) fit += X0[i * p0 + a] * bproj[a];
            gt[i] = g[i] - fit;
        }
        // Scaled to a largest |gt| of 1: the tails are the same, and the saddlepoint's checks then do not
        // depend on the dosage's units.
        double gmax = 0.0;
        for (double v : gt) gmax = std::max(gmax, std::fabs(v));
        if (gmax > 0.0) for (double& v : gt) v /= gmax;
        double S = 0.0;
        for (std::size_t i = 0; i < n; ++i) S += gt[i] * (y[i] - mu[i]);
        const double pc = few ? carriers_exact_p(S, g, gt, mu, exact_carriers, y, p0 == 1) : kNaN;
        if (std::isfinite(pc)) {
            r.p = std::max(pc, 1e-300);
            r.lp = pc > 0.0 ? std::log(pc) : std::log(r.p);   // ranked below the output floor too
            r.carriers_exact = true;
            return r;
        }
        bool exact = false;
        const double ps = spa_two_sided(S, gt, mu, &y, &exact);
        if (std::isfinite(ps) && ps > 0.0) {
            r.p = std::max(ps, 1e-300);
            r.lp = std::log(ps);
            r.spa = !exact;
            r.exact_tail = exact;
        } else {
            // No correction could be formed: the normal tail stands, labelled.
            r.normal_fallback = true;
        }
    }
    return r;
}

// Case/control, a conditioning column (columns `first`.. of the n x p row-major design X) on which no case
// lies below any control, or none above, separates them: the null fit has no maximum, and as it diverges the
// fitted chances of the samples off the boundary value tend to their own outcomes, so that those samples add
// nothing to the score or to its variance. The conditional test is that limit: the samples at the boundary
// value are kept, and a column then constant on them, the intercept again, is left out. A 0/1 column beside
// it (a missingness indicator, the column constant where it is 1) lets those samples take any value, so they
// are kept whatever side they would fall on. Repeated until no column separates. On the samples left, a
// covariate or conditioning column that repeats the intercept and the columns before it (constant there, or
// a column mean-filled where its indicator is 1, now the indicator scaled) is left out too: it adds nothing,
// and kept it would make the null fit singular. `first` then counts the columns before the conditioning
// ones that remain. False when no case or no control is left: the columns reproduce the phenotype.
std::vector<char> spanning_columns(const Eigen::MatrixXd& B, Eigen::MatrixXd& E);

bool drop_separated(std::vector<double>& X, std::vector<double>& y, std::size_t& n, std::size_t& p,
                    std::size_t& first) {
    const double inf = std::numeric_limits<double>::infinity();
    const std::size_t n_in = n;
    std::vector<char> keep;
    for (bool changed = true; changed;) {
        changed = false;
        std::size_t cases = 0;
        for (std::size_t i = 0; i < n; ++i) if (y[i] > 0.5) ++cases;
        if (cases == 0 || cases == n) return false;   // one outcome only
        for (std::size_t j = first; j < p && !changed; ++j) {
            // k == p: no free samples; else the 0/1 column k whose 1s are free
            for (std::size_t k = first; k <= p && !changed; ++k) {
                if (k == j) continue;
                if (k < p) {
                    bool binary = true, any = false, same = true; double c = 0.0;
                    for (std::size_t i = 0; i < n && binary && same; ++i) {
                        const double v = X[i * p + k];
                        if (v != 0.0 && v != 1.0) binary = false;
                        else if (v == 1.0) {
                            if (!any) { c = X[i * p + j]; any = true; } else if (X[i * p + j] != c) same = false;
                        }
                    }
                    if (!binary || !any || !same) continue;
                }
                auto is_free = [&](std::size_t i) { return k < p && X[i * p + k] == 1.0; };
                double lo1 = inf, hi1 = -inf, lo0 = inf, hi0 = -inf;
                for (std::size_t i = 0; i < n; ++i) {
                    if (is_free(i)) continue;
                    const double v = X[i * p + j];
                    if (y[i] > 0.5) { lo1 = std::min(lo1, v); hi1 = std::max(hi1, v); }
                    else { lo0 = std::min(lo0, v); hi0 = std::max(hi0, v); }
                }
                keep.assign(n, 0);
                if (hi0 < lo1 || hi1 < lo0) {   // every other sample separated
                    for (std::size_t i = 0; i < n; ++i) keep[i] = is_free(i);
                } else if (hi0 == lo1 || hi1 == lo0) {
                    const double t = hi0 == lo1 ? hi0 : hi1;
                    for (std::size_t i = 0; i < n; ++i) keep[i] = is_free(i) || X[i * p + j] == t;
                } else {
                    continue;
                }
                std::size_t m = 0;
                for (std::size_t i = 0; i < n; ++i) {
                    if (!keep[i]) continue;
                    if (m != i) {
                        std::copy(X.begin() + static_cast<std::ptrdiff_t>(i * p), X.begin() + static_cast<std::ptrdiff_t>((i + 1) * p),
                                  X.begin() + static_cast<std::ptrdiff_t>(m * p));
                        y[m] = y[i];
                    }
                    ++m;
                }
                if (m < n) { n = m; changed = true; }
            }
        }
    }
    X.resize(n * p);
    y.resize(n);
    std::vector<std::size_t> cols;
    if (n < n_in && p > 2) {
        // Covariates first, then the conditioning columns, each against the intercept and those kept before
        // it; the target (column 1) stays. The columns kept enter centred, as in spanning_columns.
        const auto ni = static_cast<Eigen::Index>(n);
        Eigen::MatrixXd B(ni, 1), E(ni, static_cast<Eigen::Index>(p - 2));
        for (std::size_t i = 0; i < n; ++i) {
            const auto r = static_cast<Eigen::Index>(i);
            B(r, 0) = X[i * p];
            for (std::size_t j = 2; j < p; ++j) E(r, static_cast<Eigen::Index>(j - 2)) = X[i * p + j];
        }
        const std::vector<char> kept = spanning_columns(B, E);
        std::size_t first_left = 2;
        cols = {0, 1};
        for (std::size_t j = 2; j < p; ++j) {
            if (!kept[j - 2]) continue;
            cols.push_back(j);
            if (j < first) ++first_left;
            for (std::size_t i = 0; i < n; ++i) X[i * p + j] = E(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j - 2));
        }
        first = first_left;
    } else {
        for (std::size_t j = 0; j < p; ++j) {
            bool varies = j < first;
            for (std::size_t i = 1; i < n && !varies; ++i) if (X[i * p + j] != X[j]) varies = true;
            if (varies) cols.push_back(j);
        }
    }
    if (cols.size() < p) {
        std::vector<double> Xc(n * cols.size());
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t k = 0; k < cols.size(); ++k) Xc[i * cols.size() + k] = X[i * p + cols[k]];
        X.swap(Xc);
        p = cols.size();
    }
    return true;
}

FitResult fit_logistic(const std::vector<double>& X, const std::vector<double>& y, std::size_t n, std::size_t p) {
    FitResult r;
    if (n <= p + 1) return r;
    std::vector<double> beta, inv;
    if (!logistic_irls(X, y, n, p, beta, inv, nullptr)) return r;   // includes non-convergence
    if (inv.size() != p * p) return r;
    const double var_b1 = inv[1 * p + 1];
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

// The 1-df chi-square value whose upper tail is p: z^2 with z = Phi^-1(1 - p/2). Acklam's rational
// approximation of the normal quantile (relative error ~1e-9), then one Halley step on erfc.
double chisq1_upper_quantile(double p) {
    if (!(p > 0.0)) return std::numeric_limits<double>::infinity();
    if (p >= 1.0) return 0.0;
    const double t = 0.5 * p;   // upper normal tail
    double z;
    if (t < 0.02425) {
        static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                                   -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
        static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                                   3.754408661907416e+00};
        const double q = std::sqrt(-2.0 * std::log(t));
        z = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    } else {
        static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                                   1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
        static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                                   6.680131188771972e+01, -1.328068155288572e+01};
        const double q = t - 0.5, r = q * q;
        z = -(((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    }
    const double e = 0.5 * std::erfc(z / std::sqrt(2.0)) - t;   // upper tail at z, minus the target
    const double u = e * std::sqrt(2.0 * M_PI) * std::exp(0.5 * z * z);
    if (std::isfinite(u)) z = z + u / (1.0 - z * u / 2.0);   // not past exp's range (p < ~1e-308)
    return z * z;
}

// Genomic-inflation factor: the median 1-df chi-square over its expected median (0.4549364), each
// chi-square read back from the reported p (which may be a t, exact or saddlepoint p).
double genomic_lambda(const std::vector<double>& ps) {
    std::vector<double> chi;
    for (double p : ps) if (std::isfinite(p)) chi.push_back(chisq1_upper_quantile(p));
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

// Which columns of E add to the span of B's columns and of the columns of E kept before them (Gram-Schmidt
// in order). B's first column is the intercept, and E's columns are first centred about it, in place: the
// span with B is the same, but a constant added to a column then changes neither this check nor the fits on
// them, which see a column of small spread on a large mean as nearly the intercept. A column left with under
// 1e-8 of its centred squared norm repeats them: kept, it would make the design singular, while dropping it
// leaves the span, and so every fit on it, unchanged.
std::vector<char> spanning_columns(const Eigen::MatrixXd& B, Eigen::MatrixXd& E) {
    const double one2 = B.cols() > 0 ? B.col(0).squaredNorm() : 0.0;
    auto centre = [&](Eigen::VectorXd& x) {
        if (!(one2 > 0.0)) return;
        const double raw = x.squaredNorm();
        x -= B.col(0) * (B.col(0).dot(x) / one2);
        if (!(x.squaredNorm() > 1e-24 * raw)) x.setZero();   // constant but for rounding
    };
    for (Eigen::Index j = 0; j < E.cols(); ++j) {
        Eigen::VectorXd x = E.col(j);
        centre(x);
        E.col(j) = x;
    }
    std::vector<Eigen::VectorXd> basis;
    auto add = [&](const Eigen::VectorXd& x) -> bool {
        const double n2 = x.squaredNorm();
        if (!(n2 > 0.0) || !std::isfinite(n2)) return false;
        Eigen::VectorXd r = x;
        for (int pass = 0; pass < 2; ++pass)   // twice, so rounding in the first pass leaves no residue
            for (const Eigen::VectorXd& q : basis) r -= q.dot(r) * q;
        const double r2 = r.squaredNorm();
        if (!(r2 > 1e-8 * n2)) return false;
        basis.push_back(r / std::sqrt(r2));
        return true;
    };
    for (Eigen::Index j = 0; j < B.cols(); ++j) {
        Eigen::VectorXd x = B.col(j);
        if (j > 0) centre(x);   // a covariate is judged about the intercept too
        add(x);
    }
    std::vector<char> keep(static_cast<std::size_t>(E.cols()), 0);
    for (Eigen::Index j = 0; j < E.cols(); ++j) keep[static_cast<std::size_t>(j)] = add(E.col(j)) ? 1 : 0;
    return keep;
}

// ---- LMM ------------------------------------------------------------------------------------------
// Eigendecompose the kinship K once, rotate phenotype + covariates, estimate delta = sigma_e^2/sigma_g^2
// under the null, then test each feature by GLS in the rotated space with weights 1/(d_i + delta), delta
// refitted with the feature in the model (GEMMA, Zhou & Stephens 2012) or kept at the null's (EMMAX, Kang
// et al. 2010).
struct LmmNull {
    Eigen::VectorXd d;     // kinship eigenvalues (n)
    Eigen::MatrixXd U;     // kinship eigenvectors (n x n)
    Eigen::MatrixXd UtX;   // U^T [intercept | covariates]  (n x q)
    Eigen::VectorXd Uty;   // U^T phenotype  (n)
    double delta = 1.0;    // sigma_e^2 / sigma_g^2 (estimated once under the null)
    bool at_boundary = false;   // delta at the edge of the searched range: no interior optimum
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
    // One checked factorisation gives the solve and the log-determinant (from its diagonal, which does
    // not underflow), and reports a singular design.
    const Eigen::LDLT<Eigen::MatrixXd> ldlt(A);
    if (ldlt.info() != Eigen::Success || !ldlt.isPositive())
        return std::numeric_limits<double>::infinity();
    const Eigen::VectorXd beta = ldlt.solve(b);
    const Eigen::VectorXd resid = m.Uty - m.UtX * beta;
    const double rss = (resid.array().square() * w).sum();
    const double sigma2 = rss / static_cast<double>(n - q);
    if (!(sigma2 > 0.0)) return std::numeric_limits<double>::infinity();
    const double logdetV = (m.d.array() + delta).log().sum();
    const Eigen::ArrayXd dA = ldlt.vectorD().array();
    if ((dA <= 0.0).any()) return std::numeric_limits<double>::infinity();
    const double logdetA = dA.log().sum();
    return static_cast<double>(n - q) * std::log(sigma2) + logdetV + logdetA;
}

// Coarse log-grid + golden-section refine for the delta minimizing the REML objective. delta is in the
// units of K, so the grid is centred on K's mean eigenvalue: K and 1e6 K give the same fit.
double estimate_delta(const LmmNull& m, bool* at_boundary = nullptr) {
    const double kscale = m.d.size() > 0 && m.d.mean() > 0.0 ? m.d.mean() : 1.0;
    const double l0 = std::log10(kscale);
    double best_l = l0, best_f = std::numeric_limits<double>::infinity();
    int best_k = 0;
    for (int k = -50; k <= 50; ++k) {
        const double l = l0 + k / 10.0;       // log10 delta in mean(d) * [1e-5, 1e5]
        const double f = lmm_reml_obj(m, std::pow(10.0, l));
        if (f < best_f) { best_f = f; best_l = l; best_k = k; }
    }
    if (at_boundary != nullptr) *at_boundary = (best_k == -50 || best_k == 50);
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

// Per-feature GLS Wald test in the rotated space (genotype is the appended last column), at variance
// ratio `delta`: the null fit's (EMMAX), or the row's own REML fit with the row among the fixed effects
// (exact, as GEMMA's -lmm 1). At one locus the exact fit costs little: a search about the null ratio.
FitResult lmm_test(const LmmNull& m, const Eigen::VectorXd& g, bool exact = true) {
    FitResult r;
    const Eigen::Index n = m.d.size();
    const Eigen::Index q = m.UtX.cols();
    if (n <= q + 2) return r;
    const Eigen::VectorXd Utg = m.U.transpose() * g;
    Eigen::MatrixXd X(n, q + 1);
    X.leftCols(q) = m.UtX;
    X.col(q) = Utg;
    double delta = m.delta;
    if (exact) {
        LmmNull mr;
        mr.d = m.d; mr.UtX = X; mr.Uty = m.Uty;
        // The whole range the null fit searches (a strong row can move the ratio by decades), then a
        // golden-section refinement about the best grid point.
        const double kscale = m.d.size() > 0 && m.d.mean() > 0.0 ? m.d.mean() : 1.0;
        const double l0 = std::log10(kscale);
        double best_l = std::log10(m.delta), best_f = std::numeric_limits<double>::infinity();
        for (int k = -50; k <= 50; ++k) {
            const double l = l0 + k / 10.0;
            const double f = lmm_reml_obj(mr, std::pow(10.0, l));
            if (f < best_f) { best_f = f; best_l = l; }
        }
        double a = best_l - 0.1, b = best_l + 0.1;
        const double gr = 0.6180339887;
        double c = b - gr * (b - a), e = a + gr * (b - a);
        double fc = lmm_reml_obj(mr, std::pow(10.0, c)), fe = lmm_reml_obj(mr, std::pow(10.0, e));
        for (int it = 0; it < 30; ++it) {
            if (fc < fe) { b = e; e = c; fe = fc; c = b - gr * (b - a); fc = lmm_reml_obj(mr, std::pow(10.0, c)); }
            else         { a = c; c = e; fc = fe; e = a + gr * (b - a); fe = lmm_reml_obj(mr, std::pow(10.0, e)); }
        }
        if (std::isfinite(std::min(fc, fe))) delta = std::pow(10.0, 0.5 * (a + b));
    }
    const Eigen::ArrayXd w = 1.0 / (m.d.array() + delta);
    const Eigen::MatrixXd Xw = X.array().colwise() * w;
    const Eigen::MatrixXd A = X.transpose() * Xw;            // (q+1) x (q+1)
    const Eigen::VectorXd b = Xw.transpose() * m.Uty;
    // A checked LDLT reports a singular design (a collinear covariate, a constant genotype). Only the
    // (q,q) entry of A^-1 is needed, solved against the last unit vector.
    const Eigen::LDLT<Eigen::MatrixXd> ldlt(A);
    if (ldlt.info() != Eigen::Success || !ldlt.isPositive()) return r;
    const Eigen::ArrayXd dA = ldlt.vectorD().array();
    if ((dA <= 0.0).any()) return r;
    const Eigen::VectorXd beta = ldlt.solve(b);
    Eigen::VectorXd e = Eigen::VectorXd::Zero(q + 1);
    e(q) = 1.0;
    const double Ainv_qq = ldlt.solve(e)(q);
    const Eigen::VectorXd resid = m.Uty - X * beta;
    const double rss = (resid.array().square() * w).sum();
    if (!(rss > 0.0)) return r;
    const double sigma2 = rss / static_cast<double>(n - (q + 1));
    const double var_b = sigma2 * Ainv_qq;
    if (!(var_b > 0.0) || !std::isfinite(var_b)) return r;
    r.beta = beta(q);
    r.se = std::sqrt(var_b);
    r.z = r.beta / r.se;
    // Student-t on n - (q+1) degrees of freedom, as the linear model.
    r.p = student_p(r.z, static_cast<double>(n - (q + 1)));
    r.lp = log_student_p(r.z, static_cast<double>(n - (q + 1)));
    r.ok = std::isfinite(r.p);
    return r;
}

// ---------------------------------------------------------------------------------------------------
// Fine-mapping: the sum of single effects (SuSiE; Wang et al. 2020, JRSS-B), fitted by iterative Bayesian
// stepwise selection on the samples' own data, with susieR's refinement. The phenotype is up to L effects,
// each one row with a normal effect; each effect's probabilities over the rows (alpha) say which row it
// is. A row's PIP is the chance that at least one effect is it; a credible set is the fewest rows holding
// --cs-coverage of one effect's probability, kept when its rows are correlated (|r| >= 0.5), so it marks
// one signal.
struct CredibleSet {
    std::vector<std::size_t> rows;   // by decreasing probability
    std::vector<double> prob;        // each row's probability within the set
    double coverage = 0.0, purity = 0.0;
};
struct SusieFit {
    std::vector<double> pip;
    std::vector<CredibleSet> sets;
    std::size_t iterations = 0;
    bool converged = false;
};

double log_sum_exp(const Eigen::ArrayXd& v) {
    const double m = v.maxCoeff();
    return m + std::log((v - m).exp().sum());
}

struct SusieState {
    Eigen::MatrixXd alpha, mu, mu2;   // L x p
    Eigen::VectorXd V;                // each effect's prior variance; 0 = no effect
    double elbo = -std::numeric_limits<double>::infinity();
    std::size_t iterations = 0;
    bool converged = false;
};

// One fit. X: n x p, columns centred and scaled to unit variance; y centred (covariates projected out of
// both, so no intercept). log_prior: each row's prior log probability of being an effect (-inf = never).
// init: a starting alpha, mu, mu2; else every effect starts spread over the rows.
SusieState susie_ibss(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, std::size_t L,
                      const Eigen::ArrayXd& log_prior, const SusieState* init,
                      std::size_t max_iter = 100, double tol = 1e-3) {
    const Eigen::Index n = X.rows(), p = X.cols();
    const double nd = static_cast<double>(n);
    const Eigen::ArrayXd d = X.colwise().squaredNorm().transpose().array();
    const double var_y = y.squaredNorm() / (nd - 1.0);
    const Eigen::ArrayXd prior = log_prior.exp();
    const auto Li = static_cast<Eigen::Index>(L);
    SusieState s;
    if (init != nullptr) {
        s.alpha = init->alpha; s.mu = init->mu; s.mu2 = init->mu2;
    } else {
        s.alpha = prior.matrix().transpose().replicate(Li, 1);
        s.mu = Eigen::MatrixXd::Zero(Li, p);
        s.mu2 = Eigen::MatrixXd::Zero(Li, p);
    }
    s.V = Eigen::VectorXd::Constant(Li, 0.2 * var_y);
    Eigen::MatrixXd Xb(n, Li);    // each effect's fitted values
    for (Eigen::Index l = 0; l < Li; ++l)
        Xb.col(l) = X * (s.alpha.row(l).array() * s.mu.row(l).array()).matrix().transpose();
    Eigen::VectorXd Xr = Xb.rowwise().sum();
    Eigen::VectorXd KL = Eigen::VectorXd::Zero(Li);
    double sigma2 = var_y;
    const double ln2pi = std::log(2.0 * M_PI);
    for (std::size_t it = 0; it < max_iter; ++it) {
        for (Eigen::Index l = 0; l < Li; ++l) {
            const Eigen::VectorXd R = y - Xr + Xb.col(l);
            const Eigen::ArrayXd XtR = (X.transpose() * R).array();
            const Eigen::ArrayXd bhat = XtR / d;
            const Eigen::ArrayXd s2 = sigma2 / d;
            // log Bayes factor of each row being the effect, under prior variance v
            auto lbf = [&](double v) {
                return (0.5 * (s2 / (v + s2)).log() + 0.5 * bhat.square() * v / (s2 * (v + s2))).eval();
            };
            auto objective = [&](double lv) { return log_sum_exp(lbf(std::exp(lv)) + log_prior); };
            // The prior variance by empirical Bayes: golden-section search on log v; none when no v beats 0.
            double a = std::log(var_y) - 20.0, b = std::log(var_y) + 5.0;
            const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
            double c1 = b - gr * (b - a), c2 = a + gr * (b - a);
            double f1 = objective(c1), f2 = objective(c2);
            for (int k = 0; k < 60; ++k) {
                if (f1 > f2) { b = c2; c2 = c1; f2 = f1; c1 = b - gr * (b - a); f1 = objective(c1); }
                else { a = c1; c1 = c2; f1 = f2; c2 = a + gr * (b - a); f2 = objective(c2); }
            }
            s.V(l) = std::max(f1, f2) > 0.0 ? std::exp(f1 > f2 ? c1 : c2) : 0.0;
            const double Rss = R.squaredNorm();
            double lbf_model = 0.0;
            if (s.V(l) > 0.0) {
                const Eigen::ArrayXd lb = lbf(s.V(l)) + log_prior;
                lbf_model = log_sum_exp(lb);
                const Eigen::ArrayXd pv = 1.0 / (1.0 / s.V(l) + d / sigma2);
                const Eigen::ArrayXd pm = pv * XtR / sigma2;
                s.alpha.row(l) = (lb - lbf_model).exp().matrix().transpose();
                s.mu.row(l) = pm.matrix().transpose();
                s.mu2.row(l) = (pv + pm.square()).matrix().transpose();
            } else {
                s.alpha.row(l) = prior.matrix().transpose();
                s.mu.row(l).setZero();
                s.mu2.row(l).setZero();
            }
            const Eigen::ArrayXd al = s.alpha.row(l).transpose().array();
            const double loglik = lbf_model - 0.5 * nd * (ln2pi + std::log(sigma2)) - 0.5 * Rss / sigma2;
            const double e_loglik = -0.5 * nd * (ln2pi + std::log(sigma2)) -
                0.5 / sigma2 * (Rss - 2.0 * (al * s.mu.row(l).transpose().array() * XtR).sum() +
                                (al * s.mu2.row(l).transpose().array() * d).sum());
            KL(l) = -loglik + e_loglik;
            const Eigen::VectorXd xb = X * (al * s.mu.row(l).transpose().array()).matrix();
            Xr += xb - Xb.col(l);
            Xb.col(l) = xb;
        }
        double erss = (y - Xr).squaredNorm();
        for (Eigen::Index l = 0; l < Li; ++l)
            erss += (s.alpha.row(l).array() * s.mu2.row(l).array() * d.transpose()).sum() - Xb.col(l).squaredNorm();
        const double elbo = -0.5 * nd * (ln2pi + std::log(sigma2)) - 0.5 * erss / sigma2 - KL.sum();
        s.iterations = it + 1;
        const bool done = elbo - s.elbo < tol;
        s.elbo = elbo;
        if (done) { s.converged = true; break; }
        sigma2 = std::max(erss / nd, 1e-12 * var_y);
    }
    return s;
}

// Credible sets of a fit: per effect, rows by decreasing probability until the coverage is reached.
// Purity is the smallest |r| between its rows: all of them up to 1,000, past that 1,000 spread evenly
// over the set. A set seen from two effects is kept once.
std::vector<CredibleSet> credible_sets(const SusieState& s, const Eigen::MatrixXd& X, double coverage,
                                       double min_purity) {
    const Eigen::Index p = X.cols();
    const Eigen::ArrayXd d = X.colwise().squaredNorm().transpose().array();
    std::vector<CredibleSet> out;
    std::set<std::vector<std::size_t>> seen;
    for (Eigen::Index l = 0; l < s.alpha.rows(); ++l) {
        if (!(s.V(l) > 0.0)) continue;
        std::vector<std::size_t> ord(static_cast<std::size_t>(p));
        std::iota(ord.begin(), ord.end(), 0);
        const auto a_of = [&](std::size_t j) { return s.alpha(l, static_cast<Eigen::Index>(j)); };
        std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) {
            return a_of(a) > a_of(b) || (a_of(a) == a_of(b) && a < b);
        });
        CredibleSet cs;
        for (std::size_t j : ord) {
            cs.rows.push_back(j);
            cs.prob.push_back(a_of(j));
            cs.coverage += a_of(j);
            if (cs.coverage >= coverage) break;
        }
        std::vector<Eigen::Index> pick;
        const std::size_t k = std::min<std::size_t>(cs.rows.size(), 1000);
        for (std::size_t a = 0; a < k; ++a)
            pick.push_back(static_cast<Eigen::Index>(cs.rows[a * cs.rows.size() / k]));
        Eigen::MatrixXd S(X.rows(), static_cast<Eigen::Index>(pick.size()));
        for (std::size_t a = 0; a < pick.size(); ++a)
            S.col(static_cast<Eigen::Index>(a)) = X.col(pick[a]) / std::sqrt(d(pick[a]));
        // Column by column, stopping at the first pair below the bar.
        cs.purity = 1.0;
        for (Eigen::Index a = 0; a + 1 < S.cols() && cs.purity >= min_purity; ++a) {
            const Eigen::VectorXd r = S.rightCols(S.cols() - a - 1).transpose() * S.col(a);
            cs.purity = std::min(cs.purity, r.cwiseAbs().minCoeff());
        }
        if (cs.purity < min_purity) continue;
        std::vector<std::size_t> key = cs.rows;
        std::sort(key.begin(), key.end());
        if (seen.insert(key).second) out.push_back(std::move(cs));
    }
    return out;
}

// The fit, then susieR's refinement: a strong row can absorb a weaker, correlated effect, leaving the
// fit in a local optimum. For each credible set, fit without its rows, refit with every row from there,
// and keep the best of these while it raises the fit's evidence (ELBO).
SusieFit susie(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, std::size_t L, double coverage,
               double min_purity = 0.5) {
    const Eigen::Index n = X.rows(), p = X.cols();
    SusieFit out;
    out.pip.assign(static_cast<std::size_t>(p), 0.0);
    if (p == 0 || n < 3 || L == 0) return out;
    const Eigen::ArrayXd uniform = Eigen::ArrayXd::Constant(p, -std::log(static_cast<double>(p)));
    SusieState best = susie_ibss(X, y, L, uniform, nullptr);
    std::vector<CredibleSet> sets = credible_sets(best, X, coverage, min_purity);
    for (;;) {
        bool improved = false;
        SusieState cand_best;
        for (const CredibleSet& cs : sets) {
            Eigen::ArrayXd w = Eigen::ArrayXd::Ones(p);
            for (std::size_t r : cs.rows) w(static_cast<Eigen::Index>(r)) = 0.0;
            if (!(w.sum() > 0.0)) continue;
            const Eigen::ArrayXd lp = (w / w.sum()).log();
            const SusieState without = susie_ibss(X, y, L, lp, nullptr);
            SusieState refit = susie_ibss(X, y, L, uniform, &without);
            if (refit.elbo > best.elbo && refit.elbo > cand_best.elbo) { cand_best = std::move(refit); improved = true; }
        }
        if (!improved) break;
        best = std::move(cand_best);
        sets = credible_sets(best, X, coverage, min_purity);
    }
    for (Eigen::Index j = 0; j < p; ++j) {
        double none = 1.0;
        for (Eigen::Index l = 0; l < best.alpha.rows(); ++l) if (best.V(l) > 0.0) none *= 1.0 - best.alpha(l, j);
        out.pip[static_cast<std::size_t>(j)] = 1.0 - none;
    }
    out.sets = std::move(sets);
    out.iterations = best.iterations;
    out.converged = best.converged;
    return out;
}

// ---------------------------------------------------------------------------------------------------
// Rare-variant group tests: burden, SKAT and SKAT-O (Wu et al. 2011; Lee et al. 2012), as the SKAT R
// package computes them (method "optimal.adj").

// Regularized upper incomplete gamma Q(a, x), by its series below a + 1 and its continued fraction above.
double gamma_q(double a, double x) {
    if (x <= 0.0) return 1.0;
    const double lg = std::lgamma(a);
    if (x < a + 1.0) {
        double sum = 1.0 / a, term = sum;
        for (int k = 1; k < 100000; ++k) {
            term *= x / (a + k);
            sum += term;
            if (std::fabs(term) < std::fabs(sum) * 1e-16) break;
        }
        return std::max(0.0, 1.0 - sum * std::exp(-x + a * std::log(x) - lg));
    }
    double b = x + 1.0 - a, c = 1.0 / 1e-300, d = 1.0 / b, h = d;
    for (int k = 1; k < 100000; ++k) {
        const double an = -k * (k - a);
        b += 2.0;
        d = an * d + b; if (std::fabs(d) < 1e-300) d = 1e-300;
        c = b + an / c; if (std::fabs(c) < 1e-300) c = 1e-300;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < 1e-16) break;
    }
    return std::exp(-x + a * std::log(x) - lg) * h;
}
double chisq_upper(double x, double df) { return gamma_q(df / 2.0, x / 2.0); }
// The x with P(chi2_df > x) = p, by bisection on log x.
double chisq_upper_quantile(double p, double df) {
    if (!(p > 0.0)) return std::numeric_limits<double>::infinity();
    if (p >= 1.0) return 0.0;
    double lo = std::log(1e-300), hi = std::log(df + 100.0);
    while (chisq_upper(std::exp(hi), df) > p) hi += 1.0;
    for (int k = 0; k < 200; ++k) {
        const double mid = 0.5 * (lo + hi);
        if (chisq_upper(std::exp(mid), df) > p) lo = mid; else hi = mid;
    }
    return std::exp(0.5 * (lo + hi));
}

// Eigenvalues of a symmetric matrix, as SKAT keeps them: the non-negative ones above 1e-5 of their mean.
std::vector<double> mixture_weights(const Eigen::MatrixXd& K) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K, Eigen::EigenvaluesOnly);
    std::vector<double> all, out;
    for (Eigen::Index k = 0; k < es.eigenvalues().size(); ++k)
        if (es.eigenvalues()(k) >= 0.0) all.push_back(es.eigenvalues()(k));
    if (all.empty()) return out;
    const double cut = std::accumulate(all.begin(), all.end(), 0.0) / static_cast<double>(all.size()) / 1e5;
    for (double l : all) if (l > cut) out.push_back(l);
    return out;
}

// P(sum_k lambda_k chi2_1 > q) by the saddlepoint approximation (Kuonen 1999): a few percent off, kept
// for the far tail. At the mean, where the formula is singular, it is interpolated between two points.
double mixture_tail_spa(const std::vector<double>& lam, double q) {
    const double lmax = *std::max_element(lam.begin(), lam.end());
    double mean = 0.0, var = 0.0;
    for (double l : lam) { mean += l; var += 2.0 * l * l; }
    auto K = [&](double t) { double s = 0.0; for (double l : lam) s -= 0.5 * std::log1p(-2.0 * l * t); return s; };
    auto K1 = [&](double t) { double s = 0.0; for (double l : lam) s += l / (1.0 - 2.0 * l * t); return s; };
    auto K2 = [&](double t) { double s = 0.0; for (double l : lam) { const double u = 1.0 - 2.0 * l * t; s += 2.0 * l * l / (u * u); } return s; };
    auto spa = [&](double x) {
        double lo = -1.0, hi = 0.5 / lmax;
        while (K1(lo) > x) lo *= 2.0;
        for (int k = 0; k < 200; ++k) {
            const double mid = 0.5 * (lo + hi);
            if (K1(mid) < x) lo = mid; else hi = mid;
        }
        const double t = 0.5 * (lo + hi);
        const double w = (t > 0.0 ? 1.0 : -1.0) * std::sqrt(std::max(0.0, 2.0 * (t * x - K(t))));
        const double v = t * std::sqrt(K2(t));
        return 0.5 * std::erfc((w + std::log(v / w) / w) / std::sqrt(2.0));
    };
    const double h = 0.02 * std::sqrt(var);
    if (std::fabs(q - mean) < h) {
        const double a = spa(mean - h), b = spa(mean + h);
        return a + (b - a) * (q - (mean - h)) / (2.0 * h);
    }
    return std::min(1.0, spa(q));
}

// 15-point Gauss-Kronrod on [a, b], with the embedded 7-point Gauss estimate for the error.
template <class F>
double gauss_kronrod(const F& f, double a, double b, double* err) {
    static const double xk[8] = {0.991455371120812639, 0.949107912342758525, 0.864864423359769073, 0.741531185599394440,
                                 0.586087235467691130, 0.405845151377397167, 0.207784955007898468, 0.0};
    static const double wk[8] = {0.022935322010529225, 0.063092092629978553, 0.104790010322250184, 0.140653259715525919,
                                 0.169004726639267903, 0.190350578064785410, 0.204432940075298892, 0.209482141084727828};
    static const double wg[4] = {0.129484966168869693, 0.279705391489276668, 0.381830050505118945, 0.417959183673469388};
    const double c = 0.5 * (a + b), h = 0.5 * (b - a);
    const double fc = f(c);
    double k = wk[7] * fc, g = wg[3] * fc;
    for (int j = 0; j < 7; ++j) {
        const double s = f(c - h * xk[j]) + f(c + h * xk[j]);
        k += wk[j] * s;
        if (j % 2 == 1) g += wg[j / 2] * s;
    }
    *err = std::fabs((k - g) * h);
    return k * h;
}
template <class F>
double adaptive_gk(const F& f, double a, double b, double tol, int max_depth = 30, int depth = 0) {
    double err;
    const double v = gauss_kronrod(f, a, b, &err);
    if (err <= tol || depth >= max_depth) return v;
    const double m = 0.5 * (a + b);
    return adaptive_gk(f, a, m, tol / 2.0, max_depth, depth + 1) + adaptive_gk(f, m, b, tol / 2.0, max_depth, depth + 1);
}

// P(sum_k lambda_k chi2_1 > q). One weight: the scaled chi-square. Otherwise Imhof's (1961) inversion,
// 1/2 + (1/pi) int_0^inf sin(theta(u)) / (u rho(u)) du, integrated up to where theta starts falling
// steadily, then between successive zeros of sin(theta), whose alternating sum Wynn's epsilon algorithm
// accelerates; its absolute error is about 1e-11, so below 1e-7 the saddlepoint takes over.
double mixture_tail(const std::vector<double>& lam_in, double q) {
    if (lam_in.empty() || !(q > 0.0)) return 1.0;
    if (lam_in.size() == 1) return std::erfc(std::sqrt(q / (2.0 * lam_in[0])));
    const double lmax = *std::max_element(lam_in.begin(), lam_in.end());
    std::vector<double> lam(lam_in.size());
    for (std::size_t k = 0; k < lam.size(); ++k) lam[k] = lam_in[k] / lmax;
    const double x = q / lmax;
    double lsum = 0.0;
    for (double l : lam) lsum += l;
    auto theta = [&](double u) { double t = 0.0; for (double l : lam) t += std::atan(l * u); return 0.5 * t - 0.5 * x * u; };
    auto f = [&](double u) {
        if (u < 1e-12) return 0.5 * (lsum - x);
        double lr = 0.0;
        for (double l : lam) lr += 0.25 * std::log1p(l * l * u * u);
        return std::sin(theta(u)) / (u * std::exp(lr));
    };
    // u0: past it theta falls at least x/4 per unit, sum_k l/(1 + l^2 u^2) <= x/2.
    auto slope_ok = [&](double u) { double s = 0.0; for (double l : lam) s += l / (1.0 + l * l * u * u); return s <= 0.5 * x; };
    double u0 = 0.0;
    if (!slope_ok(0.0)) {
        double hi = 1.0;
        while (!slope_ok(hi)) hi *= 2.0;
        double lo = 0.0;
        for (int k = 0; k < 60 && hi - lo > 1e-12 * hi; ++k) { const double mid = 0.5 * (lo + hi); if (slope_ok(mid)) hi = mid; else lo = mid; }
        u0 = hi;
    }
    double total = u0 > 0.0 ? adaptive_gk(f, 0.0, u0, 1e-13) : 0.0;
    // Zeros of sin(theta) beyond u0: theta = j pi for the integers j below theta(u0), falling.
    auto zero_at = [&](double target, double from) {
        double lo = from, hi = from + 1.0;
        while (theta(hi) > target) hi = from + 2.0 * (hi - from);
        for (int k = 0; k < 60 && hi - lo > 1e-12 * hi; ++k) { const double mid = 0.5 * (lo + hi); if (theta(mid) > target) lo = mid; else hi = mid; }
        return 0.5 * (lo + hi);
    };
    double j = std::floor(theta(u0) / M_PI);
    if (j * M_PI == theta(u0)) j -= 1.0;
    double a = u0;
    // Partial sums through each half-period, accelerated by Wynn's epsilon table.
    std::vector<double> sums;
    double last = std::numeric_limits<double>::quiet_NaN(), part = total;
    for (int k = 0; k < 400; ++k) {
        const double b = zero_at(j * M_PI, a);
        const double piece = adaptive_gk(f, a, b, 1e-14);
        part += piece;
        sums.push_back(part);
        a = b;
        j -= 1.0;
        // Pieces this small have converged; the epsilon table would only divide by rounding noise.
        if (std::fabs(piece) < 1e-15) { total = part; break; }
        if (sums.size() >= 5 && sums.size() % 2 == 1) {
            std::vector<double> e0(sums.size() + 1, 0.0), e1(sums.begin(), sums.end());
            double best = e1.back();
            for (std::size_t col = 1; e1.size() > 1; ++col) {
                std::vector<double> e2(e1.size() - 1);
                for (std::size_t i = 0; i + 1 < e1.size(); ++i) {
                    const double d = e1[i + 1] - e1[i];
                    e2[i] = e0[i + 1] + (d != 0.0 ? 1.0 / d : 1e300);
                }
                if (col % 2 == 0) best = e2.back();
                e0 = std::move(e1);
                e1 = std::move(e2);
            }
            if (std::isfinite(last) && std::fabs(best - last) < 1e-13) { total = best; break; }
            last = best;
            total = best;
        }
    }
    const double p = 0.5 + total / M_PI;
    if (!(p >= 1e-7 && p <= 1.0 + 1e-9)) return mixture_tail_spa(lam_in, q);
    return std::min(1.0, p);
}

struct GroupResult {
    double burden_p = std::numeric_limits<double>::quiet_NaN();
    double skat_p = std::numeric_limits<double>::quiet_NaN();
    double skato_p = std::numeric_limits<double>::quiet_NaN();
    double skato_rho = std::numeric_limits<double>::quiet_NaN();
};

// The mixture weights of F'F, from the smaller of F'F and F F', which share their nonzero eigenvalues:
// a group of more records than samples never forms a records-square matrix.
std::vector<double> factor_weights(const Eigen::MatrixXd& F) {
    return F.cols() <= F.rows() ? mixture_weights(F.transpose() * F) : mixture_weights(F * F.transpose());
}

// s: the weighted rows' scores (n x m design Zs, s ~ N(0, Zs'Zs) under the null).
GroupResult group_tests(const Eigen::VectorXd& s, const Eigen::MatrixXd& Zs) {
    GroupResult out;
    // A failed null model shows as non-finite scores: no p, never a small one.
    if (!s.allFinite() || !Zs.allFinite()) return out;
    const Eigen::Index m = Zs.cols();
    const Eigen::VectorXd z1 = Zs.rowwise().sum();   // Zs 1
    const double sum_s = s.sum();
    const double var_b = z1.squaredNorm();
    out.burden_p = var_b > 0.0 ? std::erfc(std::sqrt(sum_s * sum_s / var_b / 2.0)) : 1.0;
    const double q_skat = s.squaredNorm();
    out.skat_p = mixture_tail(factor_weights(Zs), q_skat);
    // SKAT-O: Q_rho = (1 - rho) Q_skat + rho Q_burden over SKAT's grid, rho = 1 taken as 0.999.
    const std::vector<double> rho = {0.0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.5, 0.999};
    const std::size_t nr = rho.size();
    std::vector<double> q_rho(nr), p_rho(nr), mu_rho(nr), var_rho(nr), df_rho(nr);
    for (std::size_t r = 0; r < nr; ++r) {
        q_rho[r] = (1.0 - rho[r]) * q_skat + rho[r] * sum_s * sum_s;
        // Zs R_rho^(1/2), R_rho = (1 - rho) I + rho 11': its Gram matrix carries Q_rho's weights.
        const double a = std::sqrt(1.0 - rho[r]);
        const double b = (std::sqrt(1.0 - rho[r] + rho[r] * static_cast<double>(m)) - a) / static_cast<double>(m);
        const Eigen::MatrixXd F = a * Zs + b * z1 * Eigen::RowVectorXd::Ones(m);
        const std::vector<double> lam = factor_weights(F);
        p_rho[r] = mixture_tail(lam, q_rho[r]);
        double c1 = 0.0, c2 = 0.0, c4 = 0.0;
        for (double l : lam) { c1 += l; c2 += l * l; c4 += l * l * l * l; }
        mu_rho[r] = c1; var_rho[r] = 2.0 * c2; df_rho[r] = c4 > 0.0 ? c2 * c2 / c4 : 1.0;
    }
    std::size_t best = 0;
    for (std::size_t r = 1; r < nr; ++r) if (p_rho[r] < p_rho[best]) best = r;
    const double pmin = p_rho[best];
    if (!std::isfinite(pmin)) return out;
    out.skato_rho = rho[best] >= 0.999 ? 1.0 : rho[best];
    // The burden direction split off: Zs = z_mean cof' + the rest (SKAT_Optimal_Param).
    const Eigen::VectorXd zm = Zs.rowwise().mean();
    const double zz = zm.squaredNorm();
    const Eigen::VectorXd cof = zz > 0.0 ? Eigen::VectorXd(Zs.transpose() * zm / zz) : Eigen::VectorXd::Zero(m);
    const Eigen::MatrixXd item2 = Zs - zm * cof.transpose();
    const std::vector<double> lam = factor_weights(item2);
    // Nothing beyond the burden direction (one record, or copies of one): every Q_rho is the same test.
    if (lam.empty()) { out.skato_p = pmin; return out; }
    // Each Q_rho's value at the smallest p, from its kurtosis-matched chi-square (SKAT's Liu "Mod").
    std::vector<double> q_at(nr);
    for (std::size_t r = 0; r < nr; ++r) {
        const double qo = chisq_upper_quantile(pmin, df_rho[r]);
        q_at[r] = (qo - df_rho[r]) / std::sqrt(2.0 * df_rho[r]) * std::sqrt(var_rho[r]) + mu_rho[r];
    }
    double mu_q = 0.0, var_q = 0.0;
    for (double l : lam) { mu_q += l; var_q += 2.0 * l * l; }
    // 4 sum((item1'item1) .* (item2'item2)), item1 = z_mean cof': 4 |z_mean|^2 |item2 cof|^2
    const double var_remain = 4.0 * zz * (item2 * cof).squaredNorm();
    var_q += var_remain;
    const double lam_sum = mu_q;
    std::vector<double> tau(nr);
    for (std::size_t r = 0; r < nr; ++r)
        tau[r] = (static_cast<double>(m * m) * rho[r] + cof.squaredNorm() * (1.0 - rho[r])) * zz;
    // p = P(chi2_1 > X) + integral over x in [0, X] of P(the rest exceeds its bound | x) dchisq(x, 1), the
    // complement form of SKAT's 1 - integral, so a small p keeps its digits. SKAT stops at X = 40; here X
    // also reaches where chi2_1 leaves 1e-4 of the smallest p, so the first term cannot floor a strong
    // group's p. x = t^2 removes the density's singularity at 0; composite 10-point Gauss-Legendre on t.
    auto tail_given = [&](double x) {
        double mn = std::numeric_limits<double>::infinity();
        for (std::size_t r = 0; r < nr; ++r) mn = std::min(mn, (q_at[r] - tau[r] * x) / (1.0 - rho[r]));
        if (mn > lam_sum * 1e4) return 0.0;
        const double st = (mn - mu_q) * std::sqrt(var_q - var_remain) / std::sqrt(var_q) + mu_q;
        return std::min(1.0, mixture_tail(lam, st));
    };
    const double xmax = std::max(40.0, chisq_upper_quantile(pmin * 1e-4, 1.0));
    // Adaptive Gauss-Kronrod on t, to an absolute error of 1e-3 of the smallest p.
    auto integrand = [&](double t) { return std::sqrt(2.0 / M_PI) * std::exp(-0.5 * t * t) * tail_given(t * t); };
    double integral = 0.0;
    const double tm = std::sqrt(xmax);
    // The tail inside is good to about 1e-11 absolute: no tolerance below that, or the splitting never ends.
    const double tol = std::max(1e-3 * pmin, 1e-10) / 16.0;
    for (int k = 0; k < 16; ++k) integral += adaptive_gk(integrand, tm * k / 16.0, tm * (k + 1) / 16.0, tol, 12);
    double p = std::erfc(std::sqrt(xmax / 2.0)) + integral;
    // SKAT's checks: never above min p times the grid size; a failed integral falls back to min p x 3.
    if (pmin * static_cast<double>(nr) < p) p = pmin * static_cast<double>(nr);
    if (!(p > 0.0)) p = pmin * 3.0;
    out.skato_p = std::min(1.0, p);
    return out;
}

// Case/control groups with few carriers: the scores s = G'(y - mu) depend only on the carriers' outcomes,
// so these are enumerated exactly at their fitted chances (as SKAT's efficient resampling, done in full),
// and each test's p is the chance of a statistic at least as extreme: burden (sum s)^2, SKAT |s|^2, and for
// SKAT-O any Q_rho at or past the value whose asymptotic tail equals the observed smallest p.
// G: n x m weighted rows (the live columns), Zs: their covariance factor, mu: the null fit's chances.
GroupResult group_exact_binary(const Eigen::MatrixXd& G, const Eigen::MatrixXd& Zs, const Eigen::VectorXd& yv,
                               const Eigen::VectorXd& mu, const GroupResult& asym, bool intercept_only) {
    GroupResult out = asym;
    std::vector<Eigen::Index> car;
    for (Eigen::Index i = 0; i < G.rows(); ++i) if (G.row(i).cwiseAbs().maxCoeff() > 0.0) car.push_back(i);
    const Eigen::Index m = G.cols();
    auto scores = [&](const std::vector<double>& yc) {
        Eigen::VectorXd sc = Eigen::VectorXd::Zero(m);
        for (std::size_t k = 0; k < car.size(); ++k) sc += G.row(car[k]).transpose() * (yc[k] - mu(car[k]));
        return sc;
    };
    std::vector<double> yobs(car.size());
    for (std::size_t k = 0; k < car.size(); ++k) yobs[k] = yv(car[k]);
    const Eigen::VectorXd so = scores(yobs);
    const double bo = so.sum() * so.sum(), qo = so.squaredNorm();
    // SKAT-O: each Q_rho's value whose asymptotic tail is the observed smallest p (bisection on the tail).
    const std::vector<double> rho = {0.0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.5, 0.999};
    const Eigen::VectorXd z1 = Zs.rowwise().sum();
    std::vector<std::vector<double>> lam(rho.size());
    double pmin = 1.0;
    for (std::size_t r = 0; r < rho.size(); ++r) {
        const double a = std::sqrt(1.0 - rho[r]);
        const double b = (std::sqrt(1.0 - rho[r] + rho[r] * static_cast<double>(m)) - a) / static_cast<double>(m);
        lam[r] = factor_weights(a * Zs + b * z1 * Eigen::RowVectorXd::Ones(m));
        pmin = std::min(pmin, mixture_tail(lam[r], (1.0 - rho[r]) * qo + rho[r] * so.sum() * so.sum()));
    }
    std::vector<double> cut(rho.size());
    for (std::size_t r = 0; r < rho.size(); ++r) {
        double lo = 0.0, hi = 1.0;
        while (mixture_tail(lam[r], hi) > pmin && hi < 1e300) hi *= 2.0;
        for (int k = 0; k < 80; ++k) { const double mid = 0.5 * (lo + hi); if (mixture_tail(lam[r], mid) > pmin) lo = mid; else hi = mid; }
        cut[r] = lo;   // the lower end: ties with the observed value count as extreme
    }
    double pb = 0.0, ps = 0.0, po = 0.0;
    const std::size_t nc = car.size();
    std::vector<double> yc(nc);
    // Without covariates, given the number of cases K: a set of s carriers being the cases has chance
    // C(n - m, K - s) / C(n, K), as for one record.
    const double nn = static_cast<double>(G.rows()), K = yv.sum(), mc = static_cast<double>(nc);
    auto lchoose = [](double a, double b) {
        return (b < 0.0 || b > a) ? -std::numeric_limits<double>::infinity()
                                  : std::lgamma(a + 1.0) - std::lgamma(b + 1.0) - std::lgamma(a - b + 1.0);
    };
    const double lall = lchoose(nn, K);
    for (std::uint64_t cfg = 0; cfg < (std::uint64_t{1} << nc); ++cfg) {
        double pr = 1.0, cases = 0.0;
        for (std::size_t k = 0; k < nc; ++k) {
            yc[k] = (cfg >> k) & 1 ? 1.0 : 0.0;
            cases += yc[k];
            pr *= yc[k] > 0.0 ? mu(car[k]) : 1.0 - mu(car[k]);
        }
        if (intercept_only) pr = std::exp(lchoose(nn - mc, K - cases) - lall);
        const Eigen::VectorXd sc = scores(yc);
        const double bq = sc.sum() * sc.sum(), qq = sc.squaredNorm();
        const double tol = 1e-9 * std::max(qo, bo) + 1e-300;
        if (bq >= bo - tol) pb += pr;
        if (qq >= qo - tol) ps += pr;
        for (std::size_t r = 0; r < rho.size(); ++r)
            if ((1.0 - rho[r]) * qq + rho[r] * bq >= cut[r] * (1.0 - 1e-9)) { po += pr; break; }
    }
    out.burden_p = std::min(1.0, pb);
    out.skat_p = std::min(1.0, ps);
    out.skato_p = std::min(1.0, po);
    return out;
}

struct Options {
    // One or more matrices, each --genotypes with the --samples given in the same position.
    std::vector<std::string> genotypes;
    std::vector<std::string> samples;
    std::vector<std::string> feature_annot;
    std::string node_genes;
    std::string phenotype;
    std::string out_prefix = "associate";
    std::string kinship;         // precomputed (genome-wide) GRM file (--kinship)
    int pca = 0;                 // top-N kinship PCs as GLM covariates (--pca)
    double min_maf = 0.01;
    std::string model = "auto";  // auto|linear|logistic|lmm
    std::string unit = "auto";   // auto|variant|feature -- the multiple-testing unit (see below)
    double ld_r2 = 0.8;          // variant tier: genotype r^2 above which a variant is an LD shadow of a lead
    long min_ac = 3;             // variant tier: samples off the modal dosage below which low_af is set
    std::size_t exact_carriers = 16;   // logistic: enumerate the score exactly up to this many carriers
    double cojo_p = -1.0;        // variant tier: forward-stepwise (COJO) entry p; <0 -> the family-wise threshold
    std::size_t permutations = 0;   // shuffles for the min-p family-wise threshold; 0 = none
    std::uint64_t seed = 1;
    std::string haplotypes;      // cosigt table: the matrices are per haplotype, combined per sample here
    std::vector<std::string> encodings{"sum"};   // sum|min|max of a sample's haplotype values
    std::size_t finemap_l = 10;  // variant unit: SuSiE effects; 0 = no fine-mapping
    double cs_coverage = 0.95;   // probability a credible set holds
    std::string group;           // rare-variant groups: gene|bubble|layer|all; empty = none
    bool lmm_fast = false;       // LMM: the null variance ratio for every row (EMMAX), not each row's own
    bool rint = false;           // rank-based inverse-normal transform of a measured phenotype
    double rare_maf = 0.01;      // a record is rare at minor_freq <= this
    bool quiet = false;
};

// The node ids in a sidecar's `nodes` cell: ids separated by ';' (k-mers, graph) or ',' (variant rows),
// and an oriented edge key "from±>to±" names its two nodes.
std::vector<std::string> node_ids_in(const std::string& s) {
    std::vector<std::string> out;
    if (s.empty() || s == ".") return out;
    const auto add_node = [&](std::string t) {
        t = trim(t);
        if (!t.empty()) out.push_back(std::move(t));
    };
    std::size_t start = 0;
    for (std::size_t p = 0; p <= s.size(); ++p) {
        if (p < s.size() && s[p] != ';' && s[p] != ',') continue;
        const std::string tok = s.substr(start, p - start);
        start = p + 1;
        const std::size_t gt = tok.find('>');
        if (gt == std::string::npos) { add_node(tok); continue; }
        for (std::string end : {tok.substr(0, gt), tok.substr(gt + 1)}) {
            if (!end.empty() && (end.back() == '+' || end.back() == '-')) end.pop_back();
            add_node(end);
        }
    }
    return out;
}

void print_help() {
    std::cout
        << "Usage:\n"
        << "  panvar associate --genotypes <bimbam.gz> --samples <samples.txt.gz> --phenotype <tsv>\n"
        << "                   [-o <prefix>] [options]\n\n"
        << "Test each row of a describe matrix for association with a phenotype: linear or logistic\n"
        << "regression, or a linear mixed model with a supplied kinship matrix. A dosage is tested as a\n"
        << "count, so a copy number is a dosage. Logistic p-values come from a score test; log_or and se\n"
        << "are a separate effect estimate.\n\n"
        << "Several matrices (--genotypes and --samples repeated, in pairs) are tested as one: rows of\n"
        << "every matrix, on the samples all of them share, matched by name. Their rows then compete in\n"
        << "the same corrections and conditional analysis. <prefix>.layers.tsv puts each layer's\n"
        << "strongest row side by side.\n\n"
        << "Options:\n"
        << "      --genotypes <path>            BIMBAM matrix, from describe or another tool; repeat for\n"
        << "                                    several (required)\n"
        << "      --samples <path>              Its column order, one name per line, e.g. describe's\n"
        << "                                    samples.txt.gz; one per --genotypes, in the same order\n"
        << "                                    (required)\n"
        << "      --haplotypes <path>           Sample table (sample, then one column per haplotype, as describe\n"
        << "                                    --samples): the matrices are describe's haplotype level and\n"
        << "                                    each sample's value is formed from its haplotypes\n"
        << "      --encoding <list>             With --haplotypes: sum, min and/or max of a sample's haplotype\n"
        << "                                    values, comma-separated; min and max add rows <id>:min and\n"
        << "                                    <id>:max (default sum)\n"
        << "      --phenotype <path>            Table with a header: sample, phenotype, then covariates;\n"
        << "                                    NA (or -9 for the phenotype) is missing; a 0/1 or PLINK 1/2\n"
        << "                                    phenotype is binary (required)\n"
        << "  -o, --out-prefix <prefix>         Output prefix (default associate)\n"
        << "      --feature-annot <path>        describe's feature_annot sidecar, or a table like it: id, layer,\n"
        << "                                    encoding, sites, nodes, then by name svtype, gene, AF, AN,\n"
        << "                                    chrom, pos, end (any present); repeat for several\n"
        << "      --node-genes <path>           call's node_genes.tsv (call --gtf): genes from each row's nodes\n"
        << "      --model <m>                   auto, linear, logistic or lmm (default auto: logistic for a\n"
        << "                                    0/1 phenotype, else linear)\n"
        << "      --kinship <path>              Genome-wide kinship matrix, n x n in the (first) --samples\n"
        << "                                    order (--haplotypes: its samples), no header, for lmm and --pca\n"
        << "      --pca <N>                     Add the top N kinship eigenvectors as covariates\n"
        << "      --lmm-fast                    LMM: the null variance ratio for every row (EMMAX) instead of\n"
        << "                                    each row's own fit (as GEMMA)\n"
        << "      --rint                        Rank-based inverse-normal transform of a measured phenotype\n"
        << "      --min-maf <X>                 Drop rows whose non-modal genotype frequency is below X\n"
        << "                                    (default 0.01)\n"
        << "      --unit <u>                    Multiple-testing unit: auto, variant or feature (default auto:\n"
        << "                                    variant when most rows are records, not k-mers, nodes,\n"
        << "                                    edges or unannotated rows)\n"
        << "      --permutations <N>            N null replicates (shuffled residuals; drawn phenotypes for\n"
        << "                                    case/control) for a family-wise min-p threshold and\n"
        << "                                    p_perm = (1 + shuffles at or below p) / (N + 1), p under the\n"
        << "                                    LMM at the null variance ratio; N >= 20 for a threshold\n"
        << "                                    (default 0 = none)\n"
        << "      --seed <N>                    Seed of the shuffles (default 1)\n"
        << "      --ld-r2 <X>                   r^2 above which a record joins a stronger one's clump (variant\n"
        << "                                    unit), and above which a conditioning missingness indicator\n"
        << "                                    copying the tested row is left out (default 0.8)\n"
        << "      --min-ac <N>                  Variant unit: flag low_af below N carriers of the minority\n"
        << "                                    genotype (default 3)\n"
        << "      --cojo-p <X>                  Variant unit: entry p of the stepwise conditional selection\n"
        << "                                    (default: the permutation threshold with 20 or more\n"
        << "                                    --permutations, else 0.05/tests)\n"
        << "      --finemap-l <N>               Variant unit: fine-map with up to N effects (SuSiE): pip,\n"
        << "                                    credible sets in <prefix>.cs.tsv; 0 = none (default 10)\n"
        << "      --cs-coverage <X>             Probability a credible set holds (default 0.95)\n"
        << "      --group <g>                   Variant unit: test rare records together per gene, bubble,\n"
        << "                                    layer or all of them (burden, SKAT, SKAT-O) into\n"
        << "                                    <prefix>.groups.tsv\n"
        << "      --rare-maf <X>                With --group: a record is rare at minor_freq <= X (default 0.01)\n"
        << "      --exact-carriers <N>          Logistic: enumerate the score's tail exactly up to N carriers,\n"
        << "                                    0 to 20 (default 16; 0 = saddlepoint only); with covariates,\n"
        << "                                    only while carriers are at most 5% of the samples, else the\n"
        << "                                    saddlepoint at any z\n"
        << "  -q, --quiet                       No progress logs\n"
        << "  -h, --help                        Show this help\n";
}

struct FeatureAnnot {
    std::string layer = ".";
    std::string bubbles = ".";
    std::string nodes = ".";
    // Columns only some sidecars carry: a record's svtype, gene, AF and AN; and a row's place on the
    // reference, chrom, pos and end.
    std::string svtype = ".";
    std::string gene = ".";
    std::string af = ".";
    std::string an = ".";
    std::string chrom = ".", pos = ".", end = ".";
};

// k-mers, nodes and edges: many correlated features per event. Rows of any other layer are records.
bool is_feature_layer(const std::string& layer) { return layer == "kmer" || layer == "graph"; }

// One BIMBAM matrix and its column order.
struct Matrix {
    std::string genotypes, samples;
    std::vector<std::string> names;   // its samples
    std::size_t width = 0;            // dosage columns per row: the --samples lines
    std::vector<std::size_t> col;     // the sample of each shared sample, in the shared order
    // --haplotypes: per sample, per haplotype, the columns of its contigs.
    std::vector<std::vector<std::vector<std::size_t>>> slots;
};

// A sample's value from its haplotypes' values: sum, smaller or larger.
double encode(const std::vector<double>& h, const std::string& how) {
    if (how == "min") return *std::min_element(h.begin(), h.end());
    if (how == "max") return *std::max_element(h.begin(), h.end());
    return std::accumulate(h.begin(), h.end(), 0.0);
}

// The id of a row written with an encoding other than the sum.
std::string encoded_id(const std::string& id, const std::string& how) {
    return how == "sum" ? id : id + ":" + how;
}

// The matrices read as one: the rows of each in turn, every row's dosages put in the order of the
// samples all matrices share, so every pass reads a row the same way. With --haplotypes each line
// becomes one row per encoding.
class GenotypeRows {
public:
    explicit GenotypeRows(const std::vector<Matrix>& m, const std::vector<std::string>& encodings = {"sum"})
        : m_(m), enc_(encodings) {}
    // id, A, B and one dosage per shared sample; false at the end of the last matrix.
    bool next(std::vector<std::string>& f) {
        if (!queue_.empty()) { f = std::move(queue_.back()); queue_.pop_back(); return true; }
        for (;;) {
            if (!reader_) {
                if (k_ >= m_.size()) return false;
                reader_ = std::make_unique<GzLineReader>(m_[k_].genotypes);
                if (!reader_->ok()) throw std::runtime_error("cannot open --genotypes: " + m_[k_].genotypes);
                line_no_ = 0;
            }
            std::string line;
            if (!reader_->getline(line)) { reader_.reset(); ++k_; continue; }
            ++line_no_;
            if (trim(line).empty()) continue;
            std::vector<std::string> raw = split(line, ',');
            check_width(raw);
            const Matrix& mx = m_[k_];
            if (!mx.slots.empty()) { combine(raw); f = std::move(queue_.back()); queue_.pop_back(); return true; }
            if (identity()) { f = std::move(raw); return true; }
            f.assign(raw.begin(), raw.begin() + 3);
            for (std::size_t c : mx.col) f.push_back(std::move(raw[3 + c]));
            return true;
        }
    }
    // "<file> line <n>", for messages about the row just read.
    std::string where() const {
        return "--genotypes " + m_[std::min(k_, m_.size() - 1)].genotypes + " line " + std::to_string(line_no_);
    }

private:
    // One haplotype-level line as a row per encoding, queued in reverse so they come out in order. A
    // haplotype is known where exactly one of its contigs is (as describe's sample level reads it), and a
    // sample where each of its haplotypes is.
    void combine(const std::vector<std::string>& raw) {
        const Matrix& mx = m_[k_];
        std::vector<double> col(mx.width);
        for (std::size_t c = 0; c < mx.width; ++c) {
            try {
                col[c] = parse_num(raw[3 + c]);
            } catch (const std::runtime_error& e) {
                throw std::runtime_error(where() + ": " + e.what());
            }
        }
        const std::string id = trim(raw[0]);
        std::vector<std::vector<std::string>> out(enc_.size());
        for (std::size_t e = 0; e < enc_.size(); ++e)
            out[e] = {encoded_id(id, enc_[e]), raw[1], raw[2]};
        std::vector<double> h;
        char buf[32];
        for (std::size_t s : mx.col) {
            h.clear();
            for (const std::vector<std::size_t>& slot : mx.slots[s]) {
                std::size_t seen = 0;
                double v = kNaN;
                for (std::size_t c : slot) if (std::isfinite(col[c])) { ++seen; v = col[c]; }
                if (seen != 1) { h.clear(); break; }
                h.push_back(v);
            }
            for (std::size_t e = 0; e < enc_.size(); ++e) {
                if (h.empty()) { out[e].emplace_back("NA"); continue; }
                std::snprintf(buf, sizeof buf, "%.17g", encode(h, enc_[e]));
                out[e].emplace_back(buf);
            }
        }
        for (std::size_t e = enc_.size(); e-- > 0;) queue_.push_back(std::move(out[e]));
    }
    bool identity() const {
        const Matrix& mx = m_[k_];
        if (mx.col.size() != mx.names.size()) return false;
        for (std::size_t c = 0; c < mx.col.size(); ++c) if (mx.col[c] != c) return false;
        return true;
    }
    // Dosage column c belongs to the c-th --samples name by position, so every row must have that
    // width; the usual cause of a mismatch is a samples.txt from the other describe folder.
    void check_width(const std::vector<std::string>& f) const {
        const Matrix& mx = m_[k_];
        const std::size_t want = 3 + mx.width;
        if (f.size() == want) return;
        std::string msg = where() + " has " + std::to_string(f.size()) + " comma-separated field(s), expected " +
                          std::to_string(want) + " (id, A, B, then one dosage per --samples name). ";
        if (f.size() < 3) {
            msg += "BIMBAM rows must be comma-separated.";
        } else {
            msg += "The row carries " + std::to_string(f.size() - 3) + " dosage(s) but --samples " +
                   mx.samples + " lists " + std::to_string(mx.width) +
                   " name(s); check that both come from the same describe folder (haplotype/ vs sample/).";
        }
        throw std::runtime_error(msg);
    }
    const std::vector<Matrix>& m_;
    std::vector<std::string> enc_;
    std::vector<std::vector<std::string>> queue_;
    std::size_t k_ = 0, line_no_ = 0;
    std::unique_ptr<GzLineReader> reader_;
};

}  // namespace

int run_associate_command(const std::vector<std::string>& args) {
    if (args.empty()) { print_help(); return 0; }
    for (const std::string& a : args)
        if (a == "-h" || a == "--help") { print_help(); return 0; }
    Options opt;
    auto need = [&](std::size_t& i) -> std::string {
        if (i + 1 >= args.size()) throw std::runtime_error("missing value for " + args[i]);
        return args[++i];
    };
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& a = args[i];
        if (a == "--genotypes") opt.genotypes.push_back(need(i));
        else if (a == "--samples") opt.samples.push_back(need(i));
        else if (a == "--feature-annot") opt.feature_annot.push_back(need(i));
        else if (a == "--node-genes") opt.node_genes = need(i);
        else if (a == "--phenotype") opt.phenotype = need(i);
        else if (a == "-o" || a == "--out-prefix") opt.out_prefix = need(i);
        else if (a == "--min-maf") opt.min_maf = cli::parse_unit_fraction_arg(a, need(i));
        else if (a == "--unit") opt.unit = need(i);
        else if (a == "--ld-r2") opt.ld_r2 = cli::parse_unit_fraction_arg(a, need(i));
        else if (a == "--min-ac") opt.min_ac = static_cast<long>(cli::parse_size_arg(a, need(i)));
        else if (a == "--exact-carriers") {
            const std::size_t v = cli::parse_size_arg(a, need(i));
            // 2^N (sum, probability) pairs per feature: N = 20 is 32 MB.
            if (v > 20) throw std::runtime_error("--exact-carriers must be in [0, 20] (2^N terms)");
            opt.exact_carriers = v;
        }
        else if (a == "--cojo-p") opt.cojo_p = cli::parse_unit_fraction_arg(a, need(i));
        else if (a == "--model") opt.model = need(i);
        else if (a == "--kinship") opt.kinship = need(i);
        else if (a == "--permutations") opt.permutations = cli::parse_size_arg(a, need(i));
        else if (a == "--seed") opt.seed = static_cast<std::uint64_t>(cli::parse_size_arg(a, need(i)));
        else if (a == "--pca") opt.pca = static_cast<int>(std::min<std::size_t>(cli::parse_size_arg(a, need(i)), 100000));
        else if (a == "--haplotypes") opt.haplotypes = need(i);
        else if (a == "--lmm-fast") opt.lmm_fast = true;
        else if (a == "--rint") opt.rint = true;
        else if (a == "--group") {
            opt.group = need(i);
            if (opt.group != "gene" && opt.group != "bubble" && opt.group != "layer" && opt.group != "all")
                throw std::runtime_error("--group must be gene, bubble, layer or all");
        }
        else if (a == "--rare-maf") opt.rare_maf = cli::parse_unit_fraction_arg(a, need(i));
        else if (a == "--finemap-l") {
            opt.finemap_l = cli::parse_size_arg(a, need(i));
            if (opt.finemap_l > 50) throw std::runtime_error("--finemap-l must be in [0, 50]");
        }
        else if (a == "--cs-coverage") {
            opt.cs_coverage = cli::parse_unit_fraction_arg(a, need(i));
            if (!(opt.cs_coverage > 0.0 && opt.cs_coverage < 1.0))
                throw std::runtime_error("--cs-coverage must be between 0 and 1");
        }
        else if (a == "--encoding") {
            opt.encodings.clear();
            for (const std::string& e : split(need(i), ',')) {
                const std::string t = trim(e);
                if (t != "sum" && t != "min" && t != "max")
                    throw std::runtime_error("--encoding takes sum, min and max, comma-separated (got '" + t + "')");
                if (std::find(opt.encodings.begin(), opt.encodings.end(), t) != opt.encodings.end())
                    throw std::runtime_error("--encoding lists '" + t + "' twice");
                opt.encodings.push_back(t);
            }
        }
        else if (a == "-q" || a == "--quiet") opt.quiet = true;
        else throw std::runtime_error("Unknown option for associate: " + a);
    }
    // A group pools each record once, as its sum.
    if (!opt.group.empty() && std::find(opt.encodings.begin(), opt.encodings.end(), "sum") == opt.encodings.end())
        throw std::runtime_error("--group pools each record's summed dosage: --encoding must include sum");
    // min and max need each haplotype's value, so the matrices must be per haplotype.
    if (opt.haplotypes.empty() && (opt.encodings.size() != 1 || opt.encodings[0] != "sum"))
        throw std::runtime_error("--encoding min or max needs --haplotypes and describe's haplotype-level "
                                 "matrices: a sample-level matrix holds only the sum");
    if (opt.genotypes.empty() || opt.samples.empty() || opt.phenotype.empty())
        throw std::runtime_error("associate requires --genotypes, --samples, --phenotype");
    if (opt.genotypes.size() != opt.samples.size())
        throw std::runtime_error("associate takes one --samples per --genotypes, in the same order (" +
                                 std::to_string(opt.genotypes.size()) + " --genotypes, " +
                                 std::to_string(opt.samples.size()) + " --samples)");
    // --genotypes is read again for the conditional pass (and for --unit auto's id pass).
    for (const std::string& g : opt.genotypes) cli::require_regular_input("associate", "--genotypes", g);
    if (!(opt.min_maf >= 0.0 && opt.min_maf <= 1.0))
        throw std::runtime_error("--min-maf must be between 0 and 1");
    if (!(opt.ld_r2 >= 0.0 && opt.ld_r2 <= 1.0))
        throw std::runtime_error("--ld-r2 must be between 0 and 1");
    if (opt.min_ac < 0) throw std::runtime_error("--min-ac must be >= 0");
    if (opt.pca < 0) throw std::runtime_error("--pca must be >= 0");
    // <0 is the "unset" sentinel meaning 0.05/Meff; anything supplied must be a probability.
    if (opt.cojo_p >= 0.0 && !(opt.cojo_p > 0.0 && opt.cojo_p <= 1.0))
        throw std::runtime_error("--cojo-p must be in (0, 1]");
    if (opt.model != "auto" && opt.model != "linear" && opt.model != "logistic" && opt.model != "lmm")
        throw std::runtime_error("--model must be auto|linear|logistic|lmm");
    if (opt.unit != "auto" && opt.unit != "variant" && opt.unit != "feature")
        throw std::runtime_error("--unit must be auto|variant|feature");
    const bool need_kinship = (opt.model == "lmm") || (opt.pca > 0);
    if (need_kinship && opt.kinship.empty())
        throw std::runtime_error("--model lmm / --pca require an external --kinship <file> "
                                 "(panvar is local and does not build a GRM; or use PC covariates)");

    cli::RunLog log("associate", opt.quiet);
    log.info("reading: " + opt.phenotype);

    // ---- sample (genotype column) order ----
    // The samples tested are those every matrix has, in the first matrix's order.
    std::vector<Matrix> matrices(opt.genotypes.size());
    for (std::size_t k = 0; k < matrices.size(); ++k) {
        Matrix& mx = matrices[k];
        mx.genotypes = opt.genotypes[k];
        mx.samples = opt.samples[k];
        GzLineReader r(mx.samples);
        if (!r.ok()) throw std::runtime_error("cannot open --samples: " + mx.samples);
        std::string line;
        while (r.getline(line)) { line = trim(line); if (!line.empty()) mx.names.push_back(line); }
        if (mx.names.empty()) throw std::runtime_error("no samples in " + mx.samples);
        mx.width = mx.names.size();
        // A repeated column would give one individual's phenotype to two genotype columns.
        std::unordered_set<std::string> seen;
        for (const std::string& sname : mx.names)
            if (!seen.insert(sname).second)
                throw std::runtime_error("--samples " + mx.samples + " has a duplicate sample id: " + sname);
    }
    // --haplotypes: each matrix's columns are haplotypes (contigs), and its samples the table's, each
    // read from its haplotypes' columns.
    if (!opt.haplotypes.empty()) {
        const CosigtTable table = read_cosigt_table(opt.haplotypes, "associate", nullptr, "", opt.quiet);
        if (table.sample_order.empty()) throw std::runtime_error("no samples in --haplotypes " + opt.haplotypes);
        for (Matrix& mx : matrices) {
            std::unordered_map<std::string, std::size_t> column;
            for (std::size_t c = 0; c < mx.names.size(); ++c) column[mx.names[c]] = c;
            mx.slots.clear();
            for (const std::string& sm : table.sample_order) {
                std::vector<std::vector<std::size_t>> cols;
                for (const std::vector<std::string>& slot : table.slots.at(sm)) {
                    cols.emplace_back();
                    for (const std::string& contig : slot) {
                        const auto it = column.find(contig);
                        if (it == column.end())
                            throw std::runtime_error("--haplotypes names '" + contig + "' (sample " + sm +
                                                     "), which --samples " + mx.samples + " does not list; "
                                                     "give describe's haplotype-level matrices");
                        cols.back().push_back(it->second);
                    }
                }
                mx.slots.push_back(std::move(cols));
            }
            mx.names = table.sample_order;
        }
    }
    std::vector<std::string> geno_samples;
    {
        std::vector<std::unordered_map<std::string, std::size_t>> index(matrices.size());
        for (std::size_t k = 0; k < matrices.size(); ++k)
            for (std::size_t c = 0; c < matrices[k].names.size(); ++c) index[k][matrices[k].names[c]] = c;
        for (const std::string& sname : matrices[0].names) {
            bool everywhere = true;
            for (const auto& ix : index) if (!ix.count(sname)) { everywhere = false; break; }
            if (!everywhere) continue;
            geno_samples.push_back(sname);
            for (std::size_t k = 0; k < matrices.size(); ++k) matrices[k].col.push_back(index[k][sname]);
        }
        if (geno_samples.empty())
            throw std::runtime_error("the --samples files share no sample name, so no sample has every matrix");
    }
    if (!opt.haplotypes.empty()) {
        std::string enc;
        for (const std::string& e : opt.encodings) enc += (enc.empty() ? "" : ",") + e;
        log.info("haplotype matrices: " + std::to_string(matrices[0].names.size()) + " samples from --haplotypes; "
                 "each row tested as the " + enc + " of a sample's haplotypes");
    }
    {
        std::string sizes;
        for (const Matrix& mx : matrices) sizes += (sizes.empty() ? "" : ", ") + std::to_string(mx.names.size());
        log.info("input: matrices: " + std::to_string(matrices.size()) + "; samples in every matrix: " +
                 std::to_string(geno_samples.size()) + (matrices.size() > 1 ? " (of " + sizes + ")" : std::string()));
    }

    // ---- phenotype + covariates: sample -> (phenotype, covariates[]) ----
    std::vector<std::string> covar_names;
    std::size_t n_pheno_minus9 = 0;
    std::unordered_map<std::string, double> pheno;
    std::unordered_map<std::string, std::vector<double>> covars;
    {
        GzLineReader r(opt.phenotype);
        if (!r.ok()) throw std::runtime_error("cannot open --phenotype: " + opt.phenotype);
        std::string line;
        bool header = true;
        std::size_t ncov = 0;
        std::size_t pline = 0;
        while (r.getline(line)) {
            ++pline;
            if (trim(line).empty()) continue;
            std::vector<std::string> f = split(line, '\t');
            if (f.size() < 2) f = split(line, ',');
            // A trailing separator is no column; an empty cell within the header's width is a missing value.
            while (f.size() > (header ? 1 : ncov + 2) && trim(f.back()).empty()) f.pop_back();
            if (header) {
                header = false;
                if (f.size() < 2)
                    throw std::runtime_error("--phenotype " + opt.phenotype + ": the header needs a sample and "
                                             "a phenotype column");
                for (std::size_t j = 2; j < f.size(); ++j) covar_names.push_back(trim(f[j]));
                ncov = covar_names.size();
                continue;
            }
            const std::string s = trim(f[0]);
            if (s.empty()) continue;
            if (f.size() != ncov + 2)
                throw std::runtime_error("--phenotype " + opt.phenotype + " line " + std::to_string(pline) +
                                         " has " + std::to_string(f.size()) + " column(s), the header " +
                                         std::to_string(ncov + 2) + "; write NA for a missing value");
            // A repeated row would make the phenotype depend on file order.
            if (pheno.find(s) != pheno.end())
                throw std::runtime_error("--phenotype has a duplicate sample id: " + s);
            std::vector<double> cv(ncov, kNaN);
            try {
                double v = parse_num(f[1]);
                if (v == -9.0) { v = kNaN; ++n_pheno_minus9; }   // PLINK's missing phenotype
                pheno[s] = v;
                for (std::size_t j = 0; j < ncov; ++j) cv[j] = parse_num(f[j + 2]);
            } catch (const std::runtime_error& e) {
                throw std::runtime_error("--phenotype " + opt.phenotype + " line " + std::to_string(pline) + ": " + e.what());
            }
            covars[s] = std::move(cv);
        }
    }

    // ---- usable samples: present in genotype order AND with non-NA phenotype + all covariates ----
    // col_to_used[c] = index into the compact phenotype/covariate arrays, or -1.
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
    {
        std::size_t not_in_matrix = 0;
        const std::unordered_set<std::string> in_matrix(geno_samples.begin(), geno_samples.end());
        for (const auto& [sname, v] : pheno) if (!in_matrix.count(sname)) ++not_in_matrix;
        if (not_in_matrix > 0 || n_pheno_minus9 > 0)
            log.info("phenotype table: samples not in --samples: " + std::to_string(not_in_matrix) +
                     "; phenotypes of -9 read as missing: " + std::to_string(n_pheno_minus9));
    }
    if (n_used < ncov + 3) throw std::runtime_error("too few usable samples after NA filtering (" +
                                                    std::to_string(n_used) + ")");

    // ---- phenotype type / model ----
    bool binary = true;
    for (double v : y) if (v != 0.0 && v != 1.0) { binary = false; break; }
    // PLINK codes control/case as 1/2: read as 0/1, the same test for the linear model and a binary
    // trait for the others.
    if (!binary) {
        bool has1 = false, has2 = false, plink = true;
        for (double v : y) {
            if (v == 1.0) has1 = true;
            else if (v == 2.0) has2 = true;
            else { plink = false; break; }
        }
        if (plink && has1 && has2) {
            for (double& v : y) v -= 1.0;
            binary = true;
            log.info("phenotype: coded 1/2, read as PLINK control/case (1 = control, 2 = case)");
        }
    }
    // --rint: a measured phenotype replaced by the normal quantile of its rank (ties share their mean
    // rank), so a few extreme values cannot drive the fits; the covariates are kept as they are.
    if (opt.rint) {
        if (binary) throw std::runtime_error("--rint transforms a measured phenotype; this one is case/control");
        std::vector<std::size_t> idx(y.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return y[a] < y[b]; });
        std::vector<double> rk(y.size());
        for (std::size_t a = 0; a < idx.size();) {
            std::size_t b = a;
            while (b + 1 < idx.size() && y[idx[b + 1]] == y[idx[a]]) ++b;
            for (std::size_t k = a; k <= b; ++k) rk[idx[k]] = 0.5 * static_cast<double>(a + b) + 1.0;
            a = b + 1;
        }
        const double nn = static_cast<double>(y.size());
        for (std::size_t u = 0; u < y.size(); ++u) {
            const double pr = (rk[u] - 0.5) / nn;
            const double z = std::sqrt(chisq_upper_quantile(2.0 * std::min(pr, 1.0 - pr), 1.0));
            y[u] = pr < 0.5 ? -z : (pr > 0.5 ? z : 0.0);
        }
        log.info("phenotype: rank-based inverse-normal transform (--rint)");
    }
    std::string model = opt.model;
    if (model == "auto") model = binary ? "logistic" : "linear";
    if (model == "logistic" && !binary)
        throw std::runtime_error("--model logistic but phenotype is not binary (0/1)");
    if (model == "lmm" && binary)
        throw std::runtime_error("--model lmm needs a quantitative phenotype (binary->logistic-LMM not yet implemented)");

    // Kinship for --model lmm / --pca: external GRM over used samples (rows/cols in --samples order).
    // A kinship from this region's own genotypes would absorb the signal it is meant to adjust for, so it
    // must come from the genome. --pca appends PCs as fixed covariates; the LMM uses K as a random effect.
    Eigen::MatrixXd K;
    Eigen::MatrixXd kinship_pcs;   // top --pca eigenvectors of K, largest eigenvalue first
    if (need_kinship) {
        std::vector<std::size_t> used_cols;
        for (std::size_t c = 0; c < geno_samples.size(); ++c) if (col_to_used[c] >= 0) used_cols.push_back(c);
        K = Eigen::MatrixXd::Zero(n_used, n_used);
        // dense n_geno x n_geno matrix (whitespace/comma), subset to used rows+cols
        std::vector<std::vector<double>> M;
        GzLineReader r(opt.kinship);
        if (!r.ok()) throw std::runtime_error("cannot open --kinship: " + opt.kinship);
        std::string line;
        std::size_t kline = 0;
        while (r.getline(line)) {
            ++kline;
            if (trim(line).empty()) continue;
            std::vector<std::string> f = split(line, '\t');
            if (f.size() < 2) f = split(line, ',');
            if (f.size() < 2) f = split(line, ' ');
            std::vector<double> row; row.reserve(f.size());
            try {
                for (const std::string& t : f) { std::string u = trim(t); if (!u.empty()) row.push_back(parse_num(u)); }
            } catch (const std::runtime_error& e) {
                throw std::runtime_error("--kinship " + opt.kinship + " line " + std::to_string(kline) + ": " + e.what());
            }
            if (!row.empty()) M.push_back(std::move(row));
        }
        // Exactly one row and one column per --samples name, in that order: the matrix has no sample ids
        // and is joined by position.
        // With several matrices, the rows follow the first --samples; with --haplotypes, its samples.
        const std::vector<std::string>& kin_names = matrices[0].names;
        const std::string kin_source = opt.haplotypes.empty() ? "name of the first --samples" : "--haplotypes sample";
        if (M.size() != kin_names.size())
            throw std::runtime_error("--kinship matrix has " + std::to_string(M.size()) + " rows, expected " +
                                     std::to_string(kin_names.size()) + " (one per " + kin_source + ", in that order)");
        for (std::size_t c = 0; c < M.size(); ++c)
            if (M[c].size() != kin_names.size())
                throw std::runtime_error("--kinship row " + std::to_string(c + 1) + " has " +
                                         std::to_string(M[c].size()) + " entries, expected " +
                                         std::to_string(kin_names.size()) + " (one per " + kin_source + ")");
        for (std::size_t a = 0; a < n_used; ++a)
            for (std::size_t b = 0; b < n_used; ++b)
                K(a, b) = M[matrices[0].col[used_cols[a]]][matrices[0].col[used_cols[b]]];

        // Finite and symmetric, as the eigendecomposition below assumes.
        double asym = 0.0, scale = 0.0;
        for (std::size_t a = 0; a < n_used; ++a) {
            for (std::size_t b = 0; b < n_used; ++b) {
                if (!std::isfinite(K(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b))))
                    throw std::runtime_error("--kinship has a non-finite entry at row " +
                                             std::to_string(a + 1) + ", column " + std::to_string(b + 1));
                scale = std::max(scale, std::fabs(K(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b))));
                if (b > a)
                    asym = std::max(asym, std::fabs(K(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) -
                                                    K(static_cast<Eigen::Index>(b), static_cast<Eigen::Index>(a))));
            }
        }
        if (scale > 0.0 && asym > 1e-6 * scale)
            throw std::runtime_error("--kinship is not symmetric (largest |K(i,j)-K(j,i)| = " +
                                     std::to_string(asym) + " against a maximum |K| of " +
                                     std::to_string(scale) + ")");
        // Positive semi-definiteness: a negative eigenvalue means it is not a covariance, and the LMM
        // variance ratio it feeds is then meaningless. Symmetrise the rounding first so the check tests
        // the matrix rather than the last bit of the input's decimals.
        for (std::size_t a = 0; a < n_used; ++a)
            for (std::size_t b = a + 1; b < n_used; ++b) {
                const double m = 0.5 * (K(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) +
                                        K(static_cast<Eigen::Index>(b), static_cast<Eigen::Index>(a)));
                K(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) = m;
                K(static_cast<Eigen::Index>(b), static_cast<Eigen::Index>(a)) = m;
            }
        {
            // With --pca the eigenvectors come from this same full decomposition.
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
                K, opt.pca > 0 ? Eigen::ComputeEigenvectors : Eigen::EigenvaluesOnly);
            if (es.info() != Eigen::Success)
                throw std::runtime_error("--kinship eigendecomposition failed; the matrix is not usable");
            const double lo = es.eigenvalues().minCoeff(), hi = es.eigenvalues().maxCoeff();
            // The tolerance scales with max(1, |hi|), so a negative definite matrix is caught too.
            const double tol = 1e-6 * std::max(1.0, std::fabs(hi));
            for (std::size_t a = 0; a < n_used; ++a)
                if (K(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(a)) < -tol)
                    throw std::runtime_error("--kinship has a negative diagonal entry at " +
                                             std::to_string(a + 1) + "; it is not a covariance matrix");
            if (lo < -tol)
                throw std::runtime_error("--kinship is not positive semi-definite (smallest eigenvalue " +
                                         std::to_string(lo) + " against a largest of " + std::to_string(hi) +
                                         "); it is not a valid GRM");
            if (opt.pca > 0) {
                const Eigen::Index k = std::min<Eigen::Index>(opt.pca, K.rows());
                // The top k eigenvectors are defined only if the k-th eigenvalue differs from the next: a
                // tie leaves which vectors of that eigenspace are taken to the sample order.
                if (k < K.rows()) {
                    const double kth = es.eigenvalues()(K.rows() - k), next = es.eigenvalues()(K.rows() - k - 1);
                    if (std::fabs(kth - next) <= 1e-8 * std::max(1.0, std::fabs(hi)))
                        throw std::runtime_error("--pca " + std::to_string(opt.pca) + ": the kinship's eigenvalue " +
                                                 std::to_string(k) + " equals eigenvalue " + std::to_string(k + 1) +
                                                 " (" + fmt_value(kth) + "), so the top " + std::to_string(k) +
                                                 " PCs are not unique; choose another number of PCs");
                }
                kinship_pcs.resize(K.rows(), k);
                for (Eigen::Index j = 0; j < k; ++j) {   // eigenvalues ascend: take them from the end
                    Eigen::VectorXd v = es.eigenvectors().col(K.rows() - 1 - j);
                    Eigen::Index arg = 0;                 // a fixed sign, for a reproducible covariate
                    v.cwiseAbs().maxCoeff(&arg);
                    if (v(arg) < 0.0) v = -v;
                    kinship_pcs.col(j) = v;
                }
            }
        }
    }

    // PCA: append the top-N kinship eigenvectors as fixed covariates (cheap structure control).
    if (opt.pca > 0) {
        const Eigen::MatrixXd& V = kinship_pcs;
        for (Eigen::Index j = 0; j < V.cols(); ++j) covar_names.push_back("kPC" + std::to_string(j + 1));
        for (std::size_t u = 0; u < n_used; ++u)
            for (Eigen::Index j = 0; j < V.cols(); ++j) Z[u].push_back(V(static_cast<Eigen::Index>(u), j));
    }
    const std::size_t ncov_eff = covar_names.size();  // base covariates (+ PCs)
    // Each covariate is centred and scaled to unit SD. Only the genotype's coefficient is reported, so no
    // output changes, and the fits' singularity tests, convergence and step limits do not depend on the
    // covariates' units or offsets.
    for (std::size_t j = 0; j < ncov_eff; ++j) {
        double m = 0.0;
        for (std::size_t u = 0; u < n_used; ++u) m += Z[u][j];
        m /= static_cast<double>(n_used);
        double ss = 0.0;
        for (std::size_t u = 0; u < n_used; ++u) { Z[u][j] -= m; ss += Z[u][j] * Z[u][j]; }
        const double sd = std::sqrt(ss / static_cast<double>(n_used));
        if (sd > 0.0)
            for (std::size_t u = 0; u < n_used; ++u) Z[u][j] /= sd;
    }

    // The phenotype and covariates, before any feature: a constant phenotype, dependent covariates, or
    // covariates that reproduce or separate the phenotype are refused.
    {
        const double y0 = y.front();
        bool constant = true;
        for (double v : y) if (v != y0) { constant = false; break; }
        if (constant)
            throw std::runtime_error("the phenotype is " + fmt_value(y0) + " for all " + std::to_string(n_used) +
                                     " usable samples; there is nothing to associate" +
                                     (binary ? " (a binary phenotype needs both cases and controls)" : ""));
        const std::size_t q = 1 + ncov_eff;
        std::vector<double> X0(n_used * q), XtX(q * q, 0.0), Xty(q, 0.0);
        for (std::size_t u = 0; u < n_used; ++u) {
            X0[u * q] = 1.0;
            for (std::size_t j = 0; j < ncov_eff; ++j) X0[u * q + 1 + j] = Z[u][j];
        }
        for (std::size_t u = 0; u < n_used; ++u)
            for (std::size_t a = 0; a < q; ++a) {
                Xty[a] += X0[u * q + a] * y[u];
                for (std::size_t b = 0; b < q; ++b) XtX[a * q + b] += X0[u * q + a] * X0[u * q + b];
            }
        std::vector<double> b0, inv0;
        if (!solve_and_invert(XtX, Xty, q, b0, inv0))
            throw std::runtime_error("the covariates are linearly dependent on each other or on the intercept "
                                     "(a constant covariate is the intercept again); drop the redundant ones");
        if (model != "logistic") {
            double ybar = 0.0;
            for (double v : y) ybar += v;
            ybar /= static_cast<double>(n_used);
            double tss = 0.0, rss = 0.0;
            for (std::size_t u = 0; u < n_used; ++u) {
                double fit = 0.0;
                for (std::size_t a = 0; a < q; ++a) fit += X0[u * q + a] * b0[a];
                rss += (y[u] - fit) * (y[u] - fit);
                tss += (y[u] - ybar) * (y[u] - ybar);
            }
            if (!(rss > 1e-10 * tss))
                throw std::runtime_error("the covariates reproduce the phenotype exactly (R^2 = 1); is the "
                                         "phenotype also given as a covariate?");
        } else {
            std::vector<double> bl, invl;
            if (!logistic_irls(X0, y, n_used, q, bl, invl, nullptr))
                throw std::runtime_error("the covariates separate the cases from the controls, so the null "
                                         "model cannot be fitted and no score test exists");
        }
    }

    // LMM null: eigendecompose K once, rotate phenotype + covariates, estimate the variance ratio.
    LmmNull lmm;
    double lmm_h2 = kNaN;
    if (model == "lmm") {
        log.info("LMM: eigendecomposing the kinship of " + std::to_string(n_used) + " samples");
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
        lmm.delta = estimate_delta(lmm, &lmm.at_boundary);
        lmm_h2 = lmm.d.mean() / (lmm.d.mean() + lmm.delta);
        log.info("LMM: heritability explained by the kinship: " + fmt_value(lmm_h2) +
                 "; variance ratio delta: " + fmt_value(lmm.delta));
        if (lmm.at_boundary)
            std::cerr << cli::log_prefix("associate") << "warning: the LMM's variance ratio sits at the edge of its search range "
                         "(heritability about " << fmt_value(lmm_h2) << "): the kinship explains "
                      << (lmm_h2 < 0.5 ? "almost none" : "almost all") << " of the phenotype, and the fit "
                         "is effectively the " << (lmm_h2 < 0.5 ? "linear model" : "kinship alone") << "\n";
    }

    // ---- optional feature annotation ----
    std::unordered_map<std::string, FeatureAnnot> annot;
    for (const std::string& annot_path : opt.feature_annot) {
        GzLineReader r(annot_path);
        if (!r.ok()) throw std::runtime_error("cannot open --feature-annot: " + annot_path);
        std::string line; bool header = true;
        // The first five columns are fixed; the others are found by name (case aside), as sidecars carry
        // different ones.
        std::unordered_map<std::string, std::size_t> col;
        const auto get = [&](const std::vector<std::string>& f, const char* name) {
            const auto it = col.find(name);
            return it != col.end() && it->second < f.size() && !f[it->second].empty() ? f[it->second] : std::string(".");
        };
        while (r.getline(line)) {
            if (line.empty()) continue;
            if (header) {
                header = false;
                const std::vector<std::string> h = split(line, '\t');
                for (std::size_t i = 5; i < h.size(); ++i) {
                    std::string k = h[i];
                    for (char& c : k) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
                    if (!k.empty() && k.back() == '\r') k.pop_back();
                    col.emplace(k, i);
                }
                continue;
            }
            std::vector<std::string> f = split(line, '\t');
            if (!f.empty() && !f.back().empty() && f.back().back() == '\r') f.back().pop_back();
            if (f.size() < 5) continue;
            FeatureAnnot a; a.layer = f[1]; a.bubbles = f[3]; a.nodes = f[4];
            a.svtype = get(f, "svtype"); a.gene = get(f, "gene"); a.af = get(f, "af"); a.an = get(f, "an");
            a.chrom = get(f, "chrom"); a.pos = get(f, "pos"); a.end = get(f, "end");
            // The same id in two sidecars, or twice in one, would leave a row's layer to file order.
            if (annot.count(f[0]) != 0)
                throw std::runtime_error("--feature-annot lists the id '" + f[0] + "' more than once");
            // An encoded row shares its record's annotation.
            for (const std::string& e : opt.encodings)
                if (e != "sum") annot.emplace(encoded_id(f[0], e), a);
            annot.emplace(f[0], std::move(a));
        }
    }

    // ---- multiple-testing unit: variant (one test per SV call) vs feature (k-mer/node/edge) ----
    // Variant mode tests one row per record, with an effective test count (Meff); feature mode tests
    // k-mers, nodes and edges, with --permutations for a family-wise threshold.
    bool variant_mode;
    if (opt.unit == "variant") variant_mode = true;
    else if (opt.unit == "feature") variant_mode = false;
    else {  // auto: variant when most genotype rows are records (not k-mers, nodes or edges)
        // Decided from the rows of the matrices, so sidecar rows of another matrix do not switch it; a row
        // without annotation counts as a feature. The pass is skipped when no sidecar row is a record.
        variant_mode = false;
        bool any_record = false;
        for (const auto& kv : annot) if (!is_feature_layer(kv.second.layer)) { any_record = true; break; }
        if (any_record) {
            GenotypeRows ids(matrices, opt.encodings);
            std::vector<std::string> f;
            std::size_t n_rec = 0, n_feat = 0, n_none = 0;
            while (ids.next(f)) {
                const auto it = annot.find(trim(f[0]));
                if (it == annot.end()) ++n_none;
                else if (is_feature_layer(it->second.layer)) ++n_feat;
                else ++n_rec;
            }
            variant_mode = n_rec > 0 && n_rec >= n_feat + n_none;
            if (!variant_mode && n_rec > 0)
                std::cerr << cli::log_prefix("associate") << "warning: of the genotype rows, " << n_rec
                          << " are records, " << n_feat << " k-mers, nodes or edges and " << n_none
                          << " unannotated; testing per feature (--unit variant forces the variant tier)\n";
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
        std::string af = ".", an = ".", svtype = ".";
        std::string chrom = ".", pos = ".", end = ".";
        int low_af = -1;          // 1/0 in variant mode (AF/AN known), -1 = unknown -> "." in output
        int clump = -1;           // LD-clump id (variant mode), -1 -> "."
        int is_lead = -1;         // 1/0 lead-of-clump (variant mode), -1 -> "."
        std::size_t n = 0;
        double minor_freq = kNaN, beta = kNaN, se = kNaN, z = kNaN, p = kNaN;
        double lp = kNaN;         // ln p, unfloored where the tail allows: rows are ranked on it
        double p_geno = kNaN, p_nonadd = kNaN;   // dosage classes as a general effect (linear model)
        double p_null_delta = kNaN;   // LMM: the p at the null variance share, the shuffles' own test
        // `ok`: the maximum-likelihood fit converged. `separation`: it did not (logistic), so log_or/se
        // come from Firth or are absent, while the score test's p stays valid.
        std::string effect_status = "ok";
        // Which tail produced `p`: "t" (Student-t, quantitative), "score" (Rao score, normal tail),
        // "score_exact_carriers" (the score's tail enumerated over its few carriers), "score_spa" (score
        // with the saddlepoint), "score_exact" (the score is at the edge of its support, where the
        // probability is a single Bernoulli configuration and is written down exactly),
        // "score_normal_fallback" (past the cutoff, but no correction could be formed), "lmm".
        std::string p_method = "t";
        // Binary traits: the samples off the modal dosage, among cases and among controls; 1 / 19 is far
        // weaker evidence than 10 / 10.
        long mac_case = -1, mac_ctrl = -1;
    };
    std::vector<Row> rows;
    // The one ranking of rows used everywhere (output order, clumps, leads, layers): the smaller
    // unfloored ln p, then the larger |stat| (tails floored at 1e-300), then the earlier row.
    auto stronger = [&rows](std::size_t a, std::size_t b) {
        const Row& ra = rows[a];
        const Row& rb = rows[b];
        const bool fa = std::isfinite(ra.p), fb = std::isfinite(rb.p);
        if (fa != fb) return fa;
        if (!fa) return a < b;
        const double la = std::isfinite(ra.lp) ? ra.lp : std::log(ra.p);
        const double lb = std::isfinite(rb.lp) ? rb.lp : std::log(rb.p);
        if (la != lb) return la < lb;
        if (std::fabs(ra.z) != std::fabs(rb.z)) return std::fabs(ra.z) > std::fabs(rb.z);
        return a < b;
    };
    struct LayerBest { bool found = false; std::size_t row = 0; Eigen::VectorXd dose; std::vector<char> obs; };
    std::map<std::string, LayerBest> layer_best;
    // Variant tier only: full-length mean-imputed dosage per retained row, for LD r^2 between variants.
    // Kept only in variant mode (few variants), so the k-mer substrate never holds a giant matrix.
    std::vector<Eigen::VectorXd> var_dose;
    // Which entries of var_dose were observed rather than mean-imputed: a conditional test uses the
    // target's observed samples, as its marginal test does.
    std::vector<std::vector<char>> var_obs;
    // --permutations: each tested row's dosage per used sample (NaN where missing).
    std::vector<std::vector<double>> perm_geno;
    std::size_t n_geno_rows = 0, n_dropped_maf = 0, n_dropped_fit = 0, n_normal_fallback = 0;
    std::size_t n_annotated = 0;   // genotype rows the --feature-annot sidecar knows
    const std::size_t p_dim = 2 + ncov_eff;  // intercept + genotype + covariates(+PCs)
    const bool is_lmm = (model == "lmm");

    // --group: the rare records, whatever --min-maf, each with its dosage's distance from the most common
    // dosage (the minor alleles, for a copy number the copies gained or lost) per used sample, NaN missing.
    struct RareRow { std::string id, gene, bubbles, layer; double minor_freq = 0.0; Eigen::VectorXd d; };
    std::vector<RareRow> rare;
    if (!opt.group.empty() && !variant_mode)
        throw std::runtime_error("--group tests records: it needs the variant unit (--unit variant)");

    // BIMBAM mean genotype: id, A, B, dose1, dose2, ...   (comma-separated, possibly space-padded)
    GenotypeRows gr(matrices, opt.encodings);
    std::vector<std::string> f;
    std::unordered_set<std::string> ids_read;
    while (gr.next(f)) {
        ++n_geno_rows;
        const auto dose_of = [&](std::size_t c) {
            try {
                return parse_num(f[3 + c]);
            } catch (const std::runtime_error& e) {
                throw std::runtime_error(gr.where() + ": " + e.what());
            }
        };
        const std::string id = trim(f[0]);
        // Every row read, tested or not: a filtered duplicate would still be found by the id lookups.
        if (!ids_read.insert(id).second)
            throw std::runtime_error("--genotypes has the feature id '" + id +
                                     "' on more than one row; ids must be unique across matrices");
        if (!annot.empty() && annot.count(id) != 0) ++n_annotated;

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
                const double dose = dose_of(c);
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
                const double dose = dose_of(c);
                if (!std::isfinite(dose)) continue;  // genuinely missing (NA = non-traversing)
                const std::size_t u = static_cast<std::size_t>(col_to_used[c]);
                g.push_back(dose);
                yy.push_back(y[u]);
                zz.push_back(Z[u]);
            }
        }
        const std::size_t n = is_lmm ? n_used : g.size();
        if (n < p_dim + 1) { ++n_dropped_fit; continue; }

        // Minor frequency: the share of observed samples whose dosage is not the most common one. It
        // reads the same for presence (0/1) and copy-number dosages, and for a rescaled matrix.
        std::map<double, std::size_t> cat;   // ordered: a tie for the mode goes to the smaller dosage
        for (double v : g) ++cat[v];
        std::size_t mode = 0;
        double modal = 0.0;
        for (const auto& kv : cat) if (kv.second > mode) { mode = kv.second; modal = kv.first; }
        const std::size_t minor_count = g.size() - mode;
        const double minor_freq = static_cast<double>(minor_count) / static_cast<double>(g.size());
        // An encoded row repeats its record, so only the sum enters a group.
        const bool encoded = !opt.haplotypes.empty() &&
            (id.size() > 4 && (id.compare(id.size() - 4, 4, ":min") == 0 || id.compare(id.size() - 4, 4, ":max") == 0));
        if (!opt.group.empty() && minor_freq > 0.0 && minor_freq <= opt.rare_maf && !encoded) {
            RareRow rr;
            rr.id = id; rr.minor_freq = minor_freq;
            std::string nodes = ".", gene = ".";
            if (auto it = annot.find(id); it != annot.end()) {
                rr.bubbles = it->second.bubbles; rr.layer = it->second.layer; nodes = it->second.nodes; gene = it->second.gene;
            } else { rr.bubbles = "."; rr.layer = "."; }
            rr.gene = (gene != "." && !gene.empty()) ? gene : genes_for(nodes);
            rr.d = Eigen::VectorXd::Constant(static_cast<Eigen::Index>(n_used), kNaN);
            for (std::size_t c = 0; c < geno_samples.size(); ++c) {
                if (col_to_used[c] < 0) continue;
                const double dose = dose_of(c);
                if (std::isfinite(dose)) rr.d(static_cast<Eigen::Index>(col_to_used[c])) = std::fabs(dose - modal);
            }
            rare.push_back(std::move(rr));
        }
        if (minor_freq < opt.min_maf) { ++n_dropped_maf; continue; }

        FitResult fr;
        std::string row_effect_status = "ok";
        std::string row_p_method = is_lmm ? "lmm" : "t";
        if (is_lmm) {
            fr = lmm_test(lmm, glmm, !opt.lmm_fast);
        } else {
            // design matrix: intercept, genotype, covariates(+PCs)
            std::vector<double> X(n * p_dim);
            for (std::size_t i = 0; i < n; ++i) {
                X[i * p_dim + 0] = 1.0;
                X[i * p_dim + 1] = g[i];
                for (std::size_t j = 0; j < ncov_eff; ++j) X[i * p_dim + 2 + j] = zz[i][j];
            }
            fr = (model == "logistic") ? fit_logistic(X, yy, n, p_dim) : fit_linear(X, yy, n, p_dim);
            if (model == "logistic") {
                bool separated = false;
                // The ML fit's beta/se are the effect size; the p is the score test's.
                const FitResult sc = score_logistic(X, yy, n, p_dim, opt.exact_carriers);
                // Without a score test there is no logistic p: the feature is dropped.
                if (!sc.ok) { fr.ok = false; }
                if (sc.ok) {
                    if (!fr.ok) {
                        // ML diverged: the effect size from Firth's penalised likelihood.
                        std::vector<double> fb, finv;
                        if (firth_logistic(X, yy, n, p_dim, fb, finv) && finv.size() == p_dim * p_dim &&
                            finv[1 * p_dim + 1] > 0.0) {
                            fr.beta = fb[1];
                            fr.se = std::sqrt(finv[1 * p_dim + 1]);
                        } else {
                            fr.beta = kNaN; fr.se = kNaN;
                        }
                        separated = true;
                    }
                    fr.z = sc.z; fr.p = sc.p; fr.lp = sc.lp; fr.ok = true;
                }
                if (separated) row_effect_status = "separation";
                row_p_method = sc.carriers_exact ? "score_exact_carriers"
                             : sc.exact_tail ? "score_exact"
                             : sc.spa ? "score_spa"
                             : sc.normal_fallback ? "score_normal_fallback" : "score";
                if (sc.normal_fallback) ++n_normal_fallback;
            }
        }
        if (!fr.ok) { ++n_dropped_fit; continue; }

        Row row;
        if (model == "linear") genotypic_test(g, yy, zz, &row.p_geno, &row.p_nonadd);
        if (is_lmm) row.p_null_delta = opt.lmm_fast ? fr.p : lmm_test(lmm, glmm, false).p;
        row.effect_status = row_effect_status;
        row.p_method = row_p_method;
        row.id = id;
        std::string ann_gene;
        if (auto it = annot.find(id); it != annot.end()) {
            row.layer = it->second.layer; row.bubbles = it->second.bubbles; row.nodes = it->second.nodes;
            row.svtype = it->second.svtype; row.af = it->second.af; row.an = it->second.an;
            row.chrom = it->second.chrom; row.pos = it->second.pos; row.end = it->second.end;
            ann_gene = it->second.gene;
        } else { row.layer = "."; row.bubbles = "."; row.nodes = "."; }
        // gene: the variant sidecar's GENES if present, else the node->gene join.
        row.gene = (ann_gene != "." && !ann_gene.empty()) ? ann_gene : genes_for(row.nodes);
        // low_af: fewer than --min-ac samples off the modal dosage, read the same for presence and
        // copy-number dosages (a copy-number record's carrier AF can be near 1 with ample variance).
        if (variant_mode)
            row.low_af = (static_cast<long>(minor_count) < opt.min_ac) ? 1 : 0;
        if (binary && g.size() == yy.size() && !g.empty()) {
            long ca = 0, co = 0;
            for (std::size_t t = 0; t < g.size(); ++t)
                if (g[t] != modal) { if (yy[t] > 0.5) ++ca; else ++co; }
            row.mac_case = ca; row.mac_ctrl = co;
        }
        // The LMM fits every usable sample, its missing dosages mean-imputed; n counts the observed ones.
        row.n = is_lmm ? g.size() : n; row.minor_freq = minor_freq;
        row.beta = fr.beta; row.se = fr.se; row.z = fr.z; row.p = fr.p; row.lp = fr.lp;
        rows.push_back(std::move(row));
        {   // Each layer's strongest row so far, with its dosage (mean-imputed) for the layers table.
            LayerBest& lb = layer_best[rows.back().layer];
            if (std::isfinite(rows.back().p) && (!lb.found || stronger(rows.size() - 1, lb.row))) {
                lb.found = true;
                lb.row = rows.size() - 1;
                lb.dose = Eigen::VectorXd::Constant(static_cast<Eigen::Index>(n_used), kNaN);
                double sum = 0.0; std::size_t nf = 0;
                for (std::size_t c = 0; c < geno_samples.size(); ++c) {
                    if (col_to_used[c] < 0) continue;
                    const double dose = dose_of(c);
                    if (!std::isfinite(dose)) continue;
                    lb.dose(static_cast<Eigen::Index>(col_to_used[c])) = dose; sum += dose; ++nf;
                }
                const double mean = nf ? sum / static_cast<double>(nf) : 0.0;
                lb.obs.assign(n_used, 1);
                for (Eigen::Index u = 0; u < lb.dose.size(); ++u)
                    if (!std::isfinite(lb.dose(u))) { lb.dose(u) = mean; lb.obs[static_cast<std::size_t>(u)] = 0; }
            }
        }
        if (opt.permutations > 0) {
            std::vector<double> full(n_used, kNaN);
            for (std::size_t c = 0; c < geno_samples.size(); ++c)
                if (col_to_used[c] >= 0) full[static_cast<std::size_t>(col_to_used[c])] = dose_of(c);
            perm_geno.push_back(std::move(full));
        }

        // Retain a full-length mean-imputed dosage for variant-tier LD-clumping (cheap: few variants).
        if (variant_mode) {
            Eigen::VectorXd vd(n_used);
            std::vector<char> obs(n_used, 1);
            if (is_lmm) {
                vd = glmm;  // already full-length, mean-imputed (the LMM rotation needs a fixed set)
                for (std::size_t c = 0; c < geno_samples.size(); ++c)
                    if (col_to_used[c] >= 0 && !std::isfinite(dose_of(c))) obs[static_cast<std::size_t>(col_to_used[c])] = 0;
            } else {
                for (std::size_t u = 0; u < n_used; ++u) vd(static_cast<Eigen::Index>(u)) =
                    std::numeric_limits<double>::quiet_NaN();
                double sum = 0.0; std::size_t nf = 0;
                for (std::size_t c = 0; c < geno_samples.size(); ++c) {
                    if (col_to_used[c] < 0) continue;
                    const double dose = dose_of(c);
                    if (!std::isfinite(dose)) continue;
                    vd(static_cast<Eigen::Index>(col_to_used[c])) = dose; sum += dose; ++nf;
                }
                const double mean = nf ? sum / static_cast<double>(nf) : 0.0;
                for (std::size_t u = 0; u < n_used; ++u)
                    if (!std::isfinite(vd(static_cast<Eigen::Index>(u)))) {
                        vd(static_cast<Eigen::Index>(u)) = mean;
                        obs[u] = 0;
                    }
            }
            var_dose.push_back(std::move(vd));
            var_obs.push_back(std::move(obs));
        }
    }
    if (n_geno_rows == 0)
        throw std::runtime_error("--genotypes contains no genotype rows");
    if (!annot.empty() && n_annotated == 0)
        throw std::runtime_error("--feature-annot annotates none of the " + std::to_string(n_geno_rows) +
                                 " genotype rows; use the feature_annot file from the same describe folder");
    if (!annot.empty() && n_annotated < n_geno_rows)
        std::cerr << cli::log_prefix("associate") << "warning: --feature-annot annotates " << n_annotated << " of the "
                  << n_geno_rows << " genotype rows; the others are reported with layer '.'\n";

    // ---- permutation min-p (--permutations) ----
    // A family-wise threshold that holds however correlated the rows are: the run is repeated on phenotypes
    // with no association but the same covariates, every row tested again, and the smallest p of each kept
    // (Westfall & Young's single-step min-p). A row's p_perm is (1 + k) / (N + 1), k the shuffles whose
    // smallest p is at or below its own, the observed data counting as one more draw. The genotypes stay with their samples and covariates, so the rows keep the
    // correlation they have given the covariates, which is the one their tests have:
    //   - linear: the residuals of the phenotype on the covariates are shuffled (Freedman & Lane);
    //   - logistic: phenotypes are drawn from the covariates-only fit and the fit redone (parametric
    //     bootstrap);
    //   - LMM: the linear scheme in the kinship-rotated coordinates where the samples are independent, with
    //     the null fit's variance ratio (as MVNpermute).
    // Shuffling the genotypes against phenotype and covariates instead loses the link between genotypes and
    // covariates, and with ancestry PCs at a differentiated locus the threshold would come out too permissive.
    std::vector<double> perm_min;
    std::vector<double> p_perm(rows.size(), kNaN);
    if (opt.permutations > 0 && !rows.empty()) {
        const std::size_t q = 1 + ncov_eff, nn = n_used;
        const bool logistic = (model == "logistic");
        Eigen::MatrixXd X0(nn, q);
        Eigen::VectorXd yv(nn);
        for (std::size_t u = 0; u < nn; ++u) {
            X0(static_cast<Eigen::Index>(u), 0) = 1.0;
            for (std::size_t j = 0; j < ncov_eff; ++j) X0(static_cast<Eigen::Index>(u), 1 + j) = Z[u][j];
            yv(static_cast<Eigen::Index>(u)) = y[u];
        }
        // LMM: everything rotated, every row mean-imputed (as its marginal test).
        Eigen::ArrayXd rw;
        if (is_lmm) {
            rw = (lmm.d.array() + lmm.delta).rsqrt();
            X0 = rw.matrix().asDiagonal() * (lmm.U.transpose() * X0);
            yv = (rw * (lmm.U.transpose() * yv).array()).matrix();
            for (auto& gr : perm_geno) {
                double sum = 0.0; std::size_t nf = 0;
                for (double v : gr) if (std::isfinite(v)) { sum += v; ++nf; }
                const double mean = nf ? sum / static_cast<double>(nf) : 0.0;
                Eigen::VectorXd gv(static_cast<Eigen::Index>(nn));
                for (std::size_t u = 0; u < nn; ++u) gv(static_cast<Eigen::Index>(u)) = std::isfinite(gr[u]) ? gr[u] : mean;
                gv = (rw * (lmm.U.transpose() * gv).array()).matrix();
                for (std::size_t u = 0; u < nn; ++u) gr[u] = gv(static_cast<Eigen::Index>(u));
            }
        }
        std::vector<double> xv(nn * q);
        for (std::size_t u = 0; u < nn; ++u)
            for (std::size_t a = 0; a < q; ++a) xv[u * q + a] = X0(static_cast<Eigen::Index>(u), a);
        // The null fit on every used sample.
        std::vector<double> mu0;
        Eigen::VectorXd fitted, resid0;
        if (logistic) {
            std::vector<double> b, inv;
            if (!logistic_irls(xv, y, nn, q, b, inv, &mu0))
                throw std::runtime_error("--permutations: the null logistic model could not be fitted");
        } else {
            fitted = X0 * (X0.transpose() * X0).ldlt().solve(X0.transpose() * yv);
            resid0 = yv - fitted;
        }
        const double df = static_cast<double>(nn) - static_cast<double>(q) - 1.0;

        // Complete rows are tested together as a matrix; rows with missing dosages one by one.
        std::vector<std::size_t> comp, part;
        for (std::size_t i = 0; i < rows.size(); ++i) {
            bool all = true;
            for (double v : perm_geno[i]) if (!std::isfinite(v)) { all = false; break; }
            (all ? comp : part).push_back(i);
        }
        Eigen::MatrixXd Gc(comp.size(), nn);
        for (std::size_t r = 0; r < comp.size(); ++r)
            for (std::size_t u = 0; u < nn; ++u)
                Gc(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(u)) = perm_geno[comp[r]][u];
        // Rows the marginal test counts exactly (few carriers) are counted exactly here too.
        std::vector<char> few(rows.size(), 0);
        if (logistic)
            for (std::size_t i = 0; i < rows.size(); ++i) {
                std::map<double, std::size_t> cnt;
                std::size_t obs = 0;
                for (double v : perm_geno[i]) if (std::isfinite(v)) { ++cnt[v]; ++obs; }
                std::size_t best = 0;
                for (const auto& kv : cnt) best = std::max(best, kv.second);
                few[i] = opt.exact_carriers > 0 && obs - best <= opt.exact_carriers;
            }
        // Linear and LMM: the genotypes adjusted for the covariates once.
        Eigen::MatrixXd R;
        if (!logistic && !comp.empty()) {
            const Eigen::MatrixXd P = (X0.transpose() * X0).ldlt().solve(X0.transpose());   // q x n
            R = Gc - (Gc * P.transpose()) * X0.transpose();
        }

        // One row's p by the marginal test's own function, on a phenotype yk (linear and logistic rows).
        const auto full_test = [&](std::size_t i, const std::vector<double>& ystar) {
            std::vector<double> X, yk;
            std::size_t nk = 0;
            for (std::size_t u = 0; u < nn; ++u) {
                const double gv = perm_geno[i][u];
                if (!std::isfinite(gv)) continue;
                X.push_back(1.0);
                X.push_back(gv);
                for (std::size_t j = 0; j < ncov_eff; ++j) X.push_back(Z[u][j]);
                yk.push_back(ystar[u]);
                ++nk;
            }
            if (nk < p_dim + 1) return kNaN;
            const FitResult fr = logistic ? score_logistic(X, yk, nk, p_dim, opt.exact_carriers)
                                          : fit_linear(X, yk, nk, p_dim);
            return fr.ok ? fr.p : kNaN;
        };

        std::mt19937_64 rng(opt.seed);
        const auto uniform = [&]() { return static_cast<double>(rng() >> 11) * 0x1.0p-53; };
        std::vector<std::size_t> perm(nn);
        std::vector<double> ystar(nn);
        perm_min.reserve(opt.permutations);
        std::size_t failed = 0;
        for (std::size_t k = 0; k < opt.permutations; ++k) {
            double best = 1.0;
            if (logistic) {
                // A phenotype from the covariates-only fit, and that fit redone on it.
                for (std::size_t u = 0; u < nn; ++u) ystar[u] = uniform() < mu0[u] ? 1.0 : 0.0;
                std::vector<double> b, inv, mu;
                if (!logistic_irls(xv, ystar, nn, q, b, inv, &mu)) { ++failed; --k; if (failed > 10 * opt.permutations) break; continue; }
                if (!comp.empty()) {
                    Eigen::VectorXd w(nn), e(nn);
                    for (std::size_t u = 0; u < nn; ++u) {
                        w(static_cast<Eigen::Index>(u)) = std::max(mu[u] * (1.0 - mu[u]), 1e-12);
                        e(static_cast<Eigen::Index>(u)) = ystar[u] - mu[u];
                    }
                    const Eigen::MatrixXd XtW = (X0.array().colwise() * w.array()).matrix().transpose();
                    const Eigen::MatrixXd Pw = (XtW * X0).ldlt().solve(XtW);
                    const Eigen::MatrixXd Rw = Gc - (Gc * Pw.transpose()) * X0.transpose();
                    const Eigen::VectorXd U = Rw * e;
                    for (std::size_t r = 0; r < comp.size(); ++r) {
                        const Eigen::Index ri = static_cast<Eigen::Index>(r);
                        const double V = (Rw.row(ri).array().square() * w.transpose().array()).sum();
                        if (!(V > 0.0)) continue;
                        const double z = U(ri) / std::sqrt(V);
                        // Where the marginal test corrects the tail, so does this.
                        const double pr = (few[comp[r]] || std::fabs(z) > 2.0) ? full_test(comp[r], ystar) : wald_p(z);
                        if (std::isfinite(pr)) best = std::min(best, pr);
                    }
                }
            } else {
                // Freedman-Lane: the null fit plus its shuffled residuals.
                std::iota(perm.begin(), perm.end(), std::size_t{0});
                for (std::size_t u = nn; u > 1; --u) {   // Fisher-Yates, the same on every platform
                    const std::uint64_t bound = u, lim = std::numeric_limits<std::uint64_t>::max() -
                                                         std::numeric_limits<std::uint64_t>::max() % bound;
                    std::uint64_t x;
                    do { x = rng(); } while (x >= lim);
                    std::swap(perm[u - 1], perm[static_cast<std::size_t>(x % bound)]);
                }
                Eigen::VectorXd e(nn);
                for (std::size_t u = 0; u < nn; ++u) e(static_cast<Eigen::Index>(u)) = resid0(static_cast<Eigen::Index>(perm[u]));
                // The shuffled residuals are not orthogonal to the covariates: their own residual is the test's.
                const Eigen::VectorXd es = e - X0 * (X0.transpose() * X0).ldlt().solve(X0.transpose() * e);
                const double yy = es.squaredNorm();
                for (std::size_t u = 0; u < nn; ++u) ystar[u] = fitted(static_cast<Eigen::Index>(u)) + e(static_cast<Eigen::Index>(u));
                if (!comp.empty()) {
                    const Eigen::VectorXd U = R * es;
                    for (std::size_t r = 0; r < comp.size(); ++r) {
                        const Eigen::Index ri = static_cast<Eigen::Index>(r);
                        const double sgg = R.row(ri).squaredNorm();
                        const double cen = (Gc.row(ri).array() - Gc.row(ri).mean()).square().sum();
                        if (!(sgg > 1e-10 * cen)) continue;   // the covariates reproduce it
                        const double rss = yy - U(ri) * U(ri) / sgg;
                        if (!(rss > 0.0)) continue;
                        const double t = (U(ri) / sgg) / std::sqrt(rss / df / sgg);
                        const double pr = student_p(t, df);
                        if (std::isfinite(pr)) best = std::min(best, pr);
                    }
                }
            }
            for (std::size_t i : part) {
                const double pr = full_test(i, ystar);
                if (std::isfinite(pr)) best = std::min(best, pr);
            }
            perm_min.push_back(best);
        }
        if (perm_min.size() < opt.permutations)
            throw std::runtime_error("--permutations: the null logistic model could not be refitted on the drawn phenotypes");
        std::vector<double> sorted_min = perm_min;
        std::sort(sorted_min.begin(), sorted_min.end());
        // The LMM's shuffles test at the null variance share (refitting it per row in every shuffle would
        // cost the run many times over), so each row is compared on that same test.
        for (std::size_t i = 0; i < rows.size(); ++i) {
            const double pr = is_lmm ? rows[i].p_null_delta : rows[i].p;
            if (!std::isfinite(pr)) continue;
            const std::size_t at_or_below = static_cast<std::size_t>(
                std::upper_bound(sorted_min.begin(), sorted_min.end(), pr) - sorted_min.begin());
            p_perm[i] = static_cast<double>(1 + at_or_below) / static_cast<double>(opt.permutations + 1);
        }
        perm_geno.clear();
    }

    // ---- multiple testing over the ACTUALLY TESTED features (not genome-wide) ----
    const std::size_t n_tests = rows.size();
    std::vector<double> pv(n_tests);
    for (std::size_t i = 0; i < n_tests; ++i) pv[i] = rows[i].p;
    const std::vector<double> qv = bh_qvalues(pv);

    // Effective number of independent tests (Meff), variant mode only: Li-Ji when feasible, else the LD
    // clumps, which also give the lead/shadow labels.
    std::size_t meff = n_tests;
    std::size_t meff_clump = 0, meff_eigen = 0;   // clumping heuristic vs the phenotype-blind estimate
    if (variant_mode && n_tests > 0) {
        std::vector<std::size_t> byp(n_tests);
        std::iota(byp.begin(), byp.end(), 0);
        std::sort(byp.begin(), byp.end(), stronger);
        int next_clump = 0;
        for (std::size_t oi = 0; oi < n_tests; ++oi) {
            const std::size_t i = byp[oi];
            if (rows[i].clump != -1 || !std::isfinite(rows[i].p)) continue;  // shadow or unfittable
            // A low_af record has an unstable r^2: it does not lead a clump, but can be another's shadow.
            if (rows[i].low_af == 1) continue;
            rows[i].clump = next_clump; rows[i].is_lead = 1;
            for (std::size_t oj = oi + 1; oj < n_tests; ++oj) {
                const std::size_t j = byp[oj];
                if (rows[j].clump != -1 || !std::isfinite(rows[j].p)) continue;
                if (r2_vec(var_dose[i], var_dose[j]) > opt.ld_r2) { rows[j].clump = next_clump; rows[j].is_lead = 0; }
            }
            ++next_clump;
        }
        // Every tested record counts toward Meff: those left in no clump get one of their own.
        for (std::size_t oi = 0; oi < n_tests; ++oi) {
            const std::size_t i = byp[oi];
            if (rows[i].clump != -1 || !std::isfinite(rows[i].p)) continue;
            rows[i].clump = next_clump++;
            rows[i].is_lead = 1;
        }
        meff_clump = static_cast<std::size_t>(std::max(1, next_clump));
        // Phenotype-blind effective tests (Li & Ji 2005): from the eigenvalues of the records' dosage
        // correlation matrix, Meff = sum_i [ I(lambda_i >= 1) + frac(lambda_i) ]. Clumping depends on the
        // order of the p-values; this does not look at the phenotype. The eigendecomposition is of a
        // min(tests, samples) square matrix; past the cap below the run keeps the clumping estimate.
        constexpr std::size_t kMeffEigenMax = 20000;
        const std::size_t eig_dim = std::min(n_tests, n_used);
        if (eig_dim > kMeffEigenMax) {
            std::cerr << cli::log_prefix("associate") << "warning: a " << eig_dim << "-square eigendecomposition is past the "
                      << kMeffEigenMax << " limit for the Li-Ji Meff; using "
                      << (meff_clump > 0 ? "the LD-clumping Meff" : "raw Bonferroni") << " (meff_method)\n";
        } else if (n_tests >= 2) {
            // Rows standardised so that S S' is the correlation matrix. With fewer samples than tests,
            // S'S (samples x samples) has the same nonzero eigenvalues and is far smaller.
            const std::size_t ns = n_used;
            Eigen::MatrixXd S(n_tests, ns);
            bool usable = true;
            for (std::size_t i = 0; i < n_tests && usable; ++i) {
                const Eigen::VectorXd& v = var_dose[i];
                const double mean = v.mean();
                const double ss = (v.array() - mean).square().sum();
                if (!(ss > 0.0)) { usable = false; break; }
                S.row(static_cast<Eigen::Index>(i)) = ((v.array() - mean) / std::sqrt(ss)).matrix().transpose();
            }
            if (usable) {
                const Eigen::MatrixXd C = n_tests <= ns ? Eigen::MatrixXd(S * S.transpose())
                                                        : Eigen::MatrixXd(S.transpose() * S);
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C, Eigen::EigenvaluesOnly);
                if (es.info() == Eigen::Success) {
                    // Li-Ji sums an integer and a fractional part, so an eigenvalue a few ulps off a
                    // whole number (k identical records give k - 1e-16 and zeros of +-1e-16) reads one
                    // too high: numeric zeros and near-integers are snapped first.
                    double m = 0.0;
                    const double scale = std::max(1.0, es.eigenvalues().cwiseAbs().maxCoeff());
                    for (Eigen::Index k = 0; k < es.eigenvalues().size(); ++k) {
                        double lam = es.eigenvalues()(k);
                        if (!(lam > 1e-8 * scale)) lam = 0.0;
                        const double near = std::round(lam);
                        if (std::fabs(lam - near) < 1e-6 * scale) lam = near;
                        m += (lam >= 1.0 ? 1.0 : 0.0) + (lam - std::floor(lam));
                    }
                    meff_eigen = std::min<std::size_t>(
                        n_tests, std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(m - 1e-9))));
                }
            }
        }
        meff = meff_eigen > 0 ? meff_eigen : meff_clump;
    } else {
        // k-mers, nodes and edges are not independent in the way an effective-test count assumes, and no
        // count tried held its error rate; the family-wise threshold here is --permutations.
        meff = 0;
    }

    // Variant-tier forward-stepwise conditional/joint (COJO), since r^2-clumping cannot establish
    // independence: a weak tag of a strong locus stays significant below --ld-r2. Add the smallest
    // conditional p until none clears --cojo-p, then report each conditioned on the rest, so shadows
    // collapse and true signals survive. Conditioning dosages enter as covariates with the genotype at
    // column 1, so FitResult is its conditional Wald; under the LMM the columns are rotated first.
    // The p below which p_perm < 0.05: a row counts at most j shuffles whose minimum is at or below it.
    double perm_threshold = kNaN;
    if (!perm_min.empty()) {
        std::vector<double> sm = perm_min;
        std::sort(sm.begin(), sm.end());
        // p_perm = (1 + k) / (B + 1) < 0.05 holds for k <= ceil((B + 1) / 20) - 2; integer arithmetic, as a
        // floating-point product rounds across the boundary when B + 1 is a multiple of 20.
        const long j = static_cast<long>((perm_min.size() + 1 + 19) / 20) - 2;
        perm_threshold = j < 0 ? 0.0 : sm[static_cast<std::size_t>(std::min<long>(j, static_cast<long>(sm.size()) - 1))];
    }
    // The p a record needs to enter the conditional selection: --cojo-p, else the family-wise threshold,
    // the permutation one when shuffles were run, else Bonferroni. Meff is not used: on dense, correlated
    // records it undercounts the tests, and a null run then selects a signal too often. Under 20 shuffles
    // no p_perm can fall below 0.05, so Bonferroni stands in.
    const bool perm_entry = !perm_min.empty() && perm_threshold > 0.0;
    const double cojo_entry = opt.cojo_p > 0.0 ? opt.cojo_p
                            : perm_entry ? perm_threshold
                            : 0.05 / static_cast<double>(std::max<std::size_t>(rows.size(), 1));
    const std::string cojo_entry_source = opt.cojo_p > 0.0 ? "--cojo-p" : perm_entry ? "permutation" : "bonferroni";
    if (variant_mode && opt.cojo_p <= 0.0 && !perm_min.empty() && !perm_entry)
        std::cerr << cli::log_prefix("associate") << "warning: " << perm_min.size() << " permutations cannot give "
                  << "p_perm below 0.05 (20 or more needed); the conditional selection uses Bonferroni\n";
    std::vector<double> p_cond(rows.size(), kNaN);
    // Samples the CONDITIONAL model used: the target's own complete cases, the marginal test's `n`.
    std::vector<std::size_t> n_cond(rows.size(), 0);
    // cond_role: variant tier -> "signal" (COJO-selected) / "shadow" (conditioned on the signals) /
    // "collinear" (the signals and covariates reproduce it on its samples, so no conditional fit exists);
    // feature tier -> "lead" / "collinear" (same-event redundancy, not scored) / "conditioned";
    // "." when not computed (e.g. no signal was selected).
    std::vector<std::string> cond_role(rows.size(), ".");
    std::size_t cojo_n_signals = 0;
    // LMM: the conditional tests run where the samples are independent, the kinship rotation scaled by the
    // null fit's variance ratio, on every sample (the marginal LMM test's set). Rotated columns are cached.
    Eigen::ArrayXd lmm_w;
    Eigen::MatrixXd lmm_C;    // rotated intercept and covariates
    Eigen::VectorXd lmm_y;
    std::unordered_map<std::size_t, Eigen::VectorXd> lmm_rot;
    auto rotate = [&](const Eigen::VectorXd& v) -> Eigen::VectorXd {
        return (lmm_w * (lmm.U.transpose() * v).array()).matrix();
    };
    if (is_lmm) {
        const auto nu = static_cast<Eigen::Index>(n_used);
        lmm_w = (lmm.d.array() + lmm.delta).rsqrt();
        Eigen::MatrixXd C(nu, 1 + static_cast<Eigen::Index>(ncov_eff));
        Eigen::VectorXd yv(nu);
        for (Eigen::Index u = 0; u < nu; ++u) {
            C(u, 0) = 1.0;
            for (std::size_t j = 0; j < ncov_eff; ++j) C(u, 1 + static_cast<Eigen::Index>(j)) = Z[static_cast<std::size_t>(u)][j];
            yv(u) = y[static_cast<std::size_t>(u)];
        }
        lmm_C = lmm_w.matrix().asDiagonal() * (lmm.U.transpose() * C);
        lmm_y = rotate(yv);
    }
    auto rotated = [&](std::size_t i) -> const Eigen::VectorXd& {
        auto it = lmm_rot.find(i);
        if (it == lmm_rot.end()) it = lmm_rot.emplace(i, rotate(var_dose[i])).first;
        return it->second;
    };
    if (variant_mode && !rows.empty()) {
        // conditional p of variant i given a set of conditioning records (empty -> its marginal test).
        auto cond_p = [&](std::size_t i, const std::vector<std::size_t>& cond, std::size_t* out_n,
                          double* out_lp = nullptr, double* out_z = nullptr) -> double {
            // The target's complete cases, as in its marginal test (every sample under the LMM); the
            // conditioning records enter mean-imputed. call leaves a record missing exactly where a larger
            // deletion covers it, so requiring them observed too would lose the deletion's carriers; a
            // conditioning record missing in some of these samples brings an indicator of where, so that
            // the missingness, itself informative, is accounted for too.
            std::vector<std::size_t> keep;
            keep.reserve(n_used);
            for (std::size_t u = 0; u < n_used; ++u) if (is_lmm || var_obs[i][u]) keep.push_back(u);
            std::size_t nk = keep.size();
            // Conditioning columns: each record, and its missingness where that varies; a column constant
            // on these samples is the intercept again and is left out.
            std::vector<std::pair<std::size_t, bool>> cols;   // (record, indicator?)
            auto varies = [&](auto value) {
                const double v0 = value(keep.empty() ? 0 : keep[0]);
                for (std::size_t u : keep) if (value(u) != v0) return true;
                return false;
            };
            // An indicator that nearly copies the target (a record nested in the target's own deletion, missing
            // exactly where the target is carried) would be the target's effect itself: it is left out.
            auto r2_with_target = [&](std::size_t c) {
                double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, nn = 0;
                for (std::size_t u : keep) {
                    const double xv = var_obs[c][u] ? 0.0 : 1.0, yv = var_dose[i](static_cast<Eigen::Index>(u));
                    sx += xv; sy += yv; sxx += xv * xv; syy += yv * yv; sxy += xv * yv; ++nn;
                }
                const double cxy = sxy - sx * sy / nn, cxx = sxx - sx * sx / nn, cyy = syy - sy * sy / nn;
                return (cxx > 0.0 && cyy > 0.0) ? cxy * cxy / (cxx * cyy) : 0.0;
            };
            for (std::size_t c : cond) {
                // The dosage enters where it varies among the samples it measures: constant there, its mean-filled
                // values would only copy its missingness indicator.
                bool dv = false;
                {
                    bool seen = false; double v0 = 0.0;
                    for (std::size_t u : keep) {
                        if (!var_obs[c][u]) continue;
                        const double v = var_dose[c](static_cast<Eigen::Index>(u));
                        if (!seen) { v0 = v; seen = true; } else if (v != v0) { dv = true; break; }
                    }
                }
                if (dv) cols.emplace_back(c, false);
                if (varies([&](std::size_t u) { return var_obs[c][u] ? 0.0 : 1.0; }) && r2_with_target(c) <= opt.ld_r2)
                    cols.emplace_back(c, true);
            }
            std::size_t pc = p_dim + cols.size();            // intercept, g_i, covariates, conditioning
            std::vector<double> X(nk * pc), yk(nk);
            if (is_lmm) {
                // Rotated design: intercept, target, covariates, then each conditioning column rotated.
                std::vector<Eigen::VectorXd> extra;
                for (const auto& [c, ind] : cols) {
                    if (!ind) { extra.push_back(rotated(c)); continue; }
                    Eigen::VectorXd m(static_cast<Eigen::Index>(n_used));
                    for (std::size_t u = 0; u < n_used; ++u) m(static_cast<Eigen::Index>(u)) = var_obs[c][u] ? 0.0 : 1.0;
                    extra.push_back(rotate(m));
                }
                const Eigen::VectorXd& gi = rotated(i);
                for (std::size_t k = 0; k < nk; ++k) {
                    const auto u = static_cast<Eigen::Index>(k);
                    std::size_t col = 0;
                    yk[k] = lmm_y(u);
                    X[k * pc + col++] = lmm_C(u, 0);
                    X[k * pc + col++] = gi(u);
                    for (std::size_t j = 0; j < ncov_eff; ++j) X[k * pc + col++] = lmm_C(u, 1 + static_cast<Eigen::Index>(j));
                    for (const Eigen::VectorXd& e : extra) X[k * pc + col++] = e(u);
                }
            } else {
                for (std::size_t k = 0; k < nk; ++k) {
                    const std::size_t u = keep[k];
                    std::size_t col = 0;
                    yk[k] = y[u];
                    X[k * pc + col++] = 1.0;
                    X[k * pc + col++] = var_dose[i](static_cast<Eigen::Index>(u));   // target -> column 1
                    for (std::size_t j = 0; j < ncov_eff; ++j) X[k * pc + col++] = Z[u][j];
                    for (const auto& [c, ind] : cols)
                        X[k * pc + col++] = ind ? (var_obs[c][u] ? 0.0 : 1.0) : var_dose[c](static_cast<Eigen::Index>(u));
                }
            }
            // Case/control: a conditioning record whose carriers are all cases (or all controls) leaves the
            // test on the other samples (drop_separated).
            std::size_t pd = p_dim;   // the columns before the conditioning ones
            if (model == "logistic" && !cols.empty() && !drop_separated(X, yk, nk, pc, pd)) return kNaN;
            // A conditioning column that repeats the covariates or the columns before it on these samples (two
            // records missing in the same samples, an indicator equal to another record's dosage) adds nothing
            // to the conditioning and is left out, so that only the target itself can make the design singular.
            // The columns kept enter centred (spanning_columns).
            if (pc > pd) {
                const std::size_t ne = pc - pd;
                const auto nke = static_cast<Eigen::Index>(nk);
                Eigen::MatrixXd B(nke, static_cast<Eigen::Index>(pd - 1)), E(nke, static_cast<Eigen::Index>(ne));
                for (std::size_t k = 0; k < nk; ++k) {
                    const auto r = static_cast<Eigen::Index>(k);
                    B(r, 0) = X[k * pc];
                    for (std::size_t j = 2; j < pd; ++j) B(r, static_cast<Eigen::Index>(j - 1)) = X[k * pc + j];
                    for (std::size_t j = 0; j < ne; ++j) E(r, static_cast<Eigen::Index>(j)) = X[k * pc + pd + j];
                }
                const std::vector<char> kept = spanning_columns(B, E);
                const std::size_t pk = pd + static_cast<std::size_t>(std::count(kept.begin(), kept.end(), 1));
                std::vector<double> Xk(nk * pk);
                for (std::size_t k = 0; k < nk; ++k) {
                    std::size_t col = 0;
                    for (std::size_t j = 0; j < pc; ++j) {
                        if (j < pd) Xk[k * pk + col++] = X[k * pc + j];
                        else if (kept[j - pd]) Xk[k * pk + col++] = E(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(j - pd));
                    }
                }
                X.swap(Xk);
                pc = pk;
            }
            if (out_n != nullptr) *out_n = nk;
            if (nk <= pc) return kNaN;   // each fit sets its own minimum beyond this
            const FitResult fc = (model == "logistic") ? score_logistic(X, yk, nk, pc, opt.exact_carriers)
                                                       : fit_linear(X, yk, nk, pc);
            if (out_z != nullptr) *out_z = fc.ok ? std::fabs(fc.z) : kNaN;
            if (out_lp != nullptr) *out_lp = fc.ok ? (std::isfinite(fc.lp) ? fc.lp : std::log(fc.p)) : kNaN;
            return fc.ok ? fc.p : kNaN;
        };
        // Forward selection of jointly independent signals, ranked by the unfloored ln p, then |stat|; after
        // each entry, as in GCTA-COJO, a selected record no longer below the bar given the others is dropped
        // for good (two partial tags of one event: the second can make the first redundant). Rows that test the
        // same given the selected ones (a sample's min and max once their sum is in, max = sum - min) differ
        // only by rounding: within a relative 1e-9 they tie, and the one stronger on its own enters.
        const auto marginal_lp = [&](std::size_t i) {
            return std::isfinite(rows[i].lp) ? rows[i].lp : std::log(rows[i].p);
        };
        const auto same = [](double a, double b) { return std::fabs(a - b) <= 1e-9 * std::max(1.0, std::fabs(b)); };
        std::vector<std::size_t> selected;
        std::vector<char> in_sel(rows.size(), 0), dropped(rows.size(), 0);
        for (std::size_t step = 0; step < 2 * rows.size(); ++step) {
            std::size_t best = rows.size();
            double best_p = kNaN, best_lp = std::numeric_limits<double>::infinity(), best_z = 0.0;
            for (std::size_t i = 0; i < rows.size(); ++i) {
                if (in_sel[i] || dropped[i] || !std::isfinite(rows[i].p)) continue;
                double z_i = 0.0, lp_i = kNaN;
                const double pc_i = cond_p(i, selected, nullptr, &lp_i, &z_i);
                if (!std::isfinite(pc_i)) continue;
                const bool tie = best < rows.size() && same(lp_i, best_lp) && same(z_i, best_z);
                if (tie ? marginal_lp(i) < marginal_lp(best) : (lp_i < best_lp || (lp_i == best_lp && z_i > best_z))) {
                    best_p = pc_i; best_lp = lp_i; best_z = z_i; best = i;
                }
            }
            if (best == rows.size() || !(best_p < cojo_entry)) break;   // nothing new clears the bar
            selected.push_back(best); in_sel[best] = 1;
            while (selected.size() > 1) {
                std::size_t worst = rows.size();
                double worst_p = -1.0;
                for (std::size_t s : selected) {
                    std::vector<std::size_t> others;
                    for (std::size_t o : selected) if (o != s) others.push_back(o);
                    const double ps = cond_p(s, others, nullptr);
                    const double pv = std::isfinite(ps) ? ps : 1.0;   // no fit beside the others: redundant
                    if (pv > worst_p) { worst_p = pv; worst = s; }
                }
                if (!(worst_p >= cojo_entry)) break;
                selected.erase(std::find(selected.begin(), selected.end(), worst));
                in_sel[worst] = 0; dropped[worst] = 1;
            }
        }
        cojo_n_signals = selected.size();
        // report each variant conditioned on the selected set minus itself
        for (std::size_t i = 0; i < rows.size(); ++i) {
            if (!std::isfinite(rows[i].p) || selected.empty()) continue;   // no signal: nothing to condition on
            std::vector<std::size_t> cond;
            for (std::size_t c : selected) if (c != i) cond.push_back(c);
            if (!cond.empty()) p_cond[i] = cond_p(i, cond, &n_cond[i]);   // empty cond (sole signal) -> NA
            cond_role[i] = in_sel[i] ? "signal" : std::isfinite(p_cond[i]) ? "shadow" : "collinear";
        }
    } else if (!variant_mode && !is_lmm && !rows.empty()) {
        // Feature tier: each feature conditioned on the strongest one. Features with r^2 > 0.95 to it
        // measure the same event and are marked `collinear` rather than tested; a mere tag of the lead's
        // event collapses. Two streaming passes, no feature x sample matrix kept.
        auto parse_full = [&](const std::vector<std::string>& f, Eigen::VectorXd& vd,
                              std::vector<char>* obs = nullptr) -> bool {
            vd = Eigen::VectorXd::Constant(static_cast<Eigen::Index>(n_used), kNaN);
            double sum = 0.0; std::size_t nf = 0;
            for (std::size_t c = 0; c < geno_samples.size(); ++c) {
                if (col_to_used[c] < 0) continue;
                const double dose = parse_num(f[3 + c]);   // validated by the first pass
                if (!std::isfinite(dose)) continue;
                vd(static_cast<Eigen::Index>(col_to_used[c])) = dose; sum += dose; ++nf;
            }
            if (nf == 0) return false;
            const double mean = sum / static_cast<double>(nf);
            if (obs != nullptr) obs->assign(n_used, 1);
            for (std::size_t u = 0; u < n_used; ++u)
                if (!std::isfinite(vd(static_cast<Eigen::Index>(u)))) {
                    vd(static_cast<Eigen::Index>(u)) = mean;
                    if (obs != nullptr) (*obs)[u] = 0;
                }
            return true;
        };
        std::size_t lead_idx = 0;
        for (std::size_t i = 1; i < rows.size(); ++i)
            if (stronger(i, lead_idx)) lead_idx = i;
        const std::string lead_id = rows[lead_idx].id;
        std::unordered_map<std::string, std::size_t> id_to_row;
        for (std::size_t i = 0; i < rows.size(); ++i) id_to_row[rows[i].id] = i;
        Eigen::VectorXd gl;  // pass 1: capture the lead feature's dosage, and where it is missing
        std::vector<char> gl_obs;

        { GenotypeRows g2(matrices, opt.encodings); std::vector<std::string> f;
          while (g2.next(f))
              if (trim(f[0]) == lead_id) { parse_full(f, gl, &gl_obs); break; } }
        if (gl.size() == static_cast<Eigen::Index>(n_used)) {
            GenotypeRows g3(matrices, opt.encodings); std::vector<std::string> f;
            while (g3.next(f)) {
                const auto it = id_to_row.find(trim(f[0]));
                if (it == id_to_row.end()) continue;
                const std::size_t i = it->second;
                if (i == lead_idx) { cond_role[i] = "lead"; continue; }
                Eigen::VectorXd vd;
                std::vector<char> obs;
                if (!parse_full(f, vd, &obs)) continue;
                if (r2_vec(vd, gl) > 0.95) { cond_role[i] = "collinear"; continue; }  // same-event redundancy
                std::vector<std::size_t> keep;
                keep.reserve(n_used);
                // The target's own samples, as in its marginal test; the lead enters mean-imputed.
                for (std::size_t u = 0; u < n_used; ++u)
                    if (obs[u]) keep.push_back(u);
                std::size_t nk = keep.size();
                // The lead's missingness, where it varies on these samples, enters as an indicator too.
                bool ind = false;
                for (std::size_t k = 1; k < nk && !ind; ++k) if (gl_obs[keep[k]] != gl_obs[keep[0]]) ind = true;
                if (ind) {   // left out when it nearly copies the target, as in the variant unit
                    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
                    for (std::size_t u : keep) {
                        const double xv = gl_obs[u] ? 0.0 : 1.0, yv = vd(static_cast<Eigen::Index>(u));
                        sx += xv; sy += yv; sxx += xv * xv; syy += yv * yv; sxy += xv * yv;
                    }
                    const double nn = static_cast<double>(nk);
                    const double cxy = sxy - sx * sy / nn, cxx = sxx - sx * sx / nn, cyy = syy - sy * sy / nn;
                    if (cxx > 0.0 && cyy > 0.0 && cxy * cxy / (cxx * cyy) > opt.ld_r2) ind = false;
                }
                // A lead constant on these samples (a deletion its target is measured only outside of)
                // carries nothing here: conditioning on it changes nothing, so it is left out rather than
                // making the design singular.
                // Varying is judged where the lead is observed: its mean-filled values would only copy the indicator.
                bool lead_var = false;
                {
                    bool seen = false; double v0 = 0.0;
                    for (std::size_t u : keep) {
                        if (!gl_obs[u]) continue;
                        const double v = gl(static_cast<Eigen::Index>(u));
                        if (!seen) { v0 = v; seen = true; } else if (v != v0) { lead_var = true; break; }
                    }
                }
                // intercept, g_i, covariates[, g_lead][, missing]
                std::size_t pc = p_dim + (lead_var ? 1 : 0) + (ind ? 1 : 0);
                if (nk <= pc) continue;
                std::vector<double> X(nk * pc), yk(nk);
                for (std::size_t k = 0; k < nk; ++k) {
                    const std::size_t u = keep[k];
                    yk[k] = y[u];
                    X[k * pc + 0] = 1.0;
                    X[k * pc + 1] = vd(static_cast<Eigen::Index>(u));
                    for (std::size_t j = 0; j < ncov_eff; ++j) X[k * pc + 2 + j] = Z[u][j];
                    std::size_t c = 2 + ncov_eff;
                    if (lead_var) X[k * pc + c++] = gl(static_cast<Eigen::Index>(u));
                    if (ind) X[k * pc + c] = gl_obs[u] ? 0.0 : 1.0;
                }
                // A lead whose carriers are all cases (or all controls): the test on the other samples.
                std::size_t pd = p_dim;
                if (model == "logistic" && pc > p_dim && !drop_separated(X, yk, nk, pc, pd)) continue;
                const FitResult fc = (model == "logistic") ? score_logistic(X, yk, nk, pc, opt.exact_carriers)
                                                           : fit_linear(X, yk, nk, pc);
                if (fc.ok) { p_cond[i] = fc.p; n_cond[i] = nk; cond_role[i] = "conditioned"; }
            }
        }
    }

    // ---- fine-mapping (variant unit) ----
    // The records' mean-imputed dosages and the phenotype, the covariates projected out of both; for the
    // LMM in the kinship-rotated coordinates where the samples are independent, and for case/control the
    // linear model on the 0/1 phenotype.
    std::vector<double> pip(rows.size(), kNaN);
    std::vector<int> cs_of(rows.size(), -1);
    SusieFit fm;
    std::vector<std::size_t> fm_rows;   // the rows fitted, in the fit's column order
    const bool finemap = variant_mode && opt.finemap_l > 0 && n_tests > 0;
    if (finemap) {
        const auto nu = static_cast<Eigen::Index>(n_used);
        Eigen::MatrixXd C(nu, 1 + static_cast<Eigen::Index>(ncov_eff));
        Eigen::VectorXd yf(nu);
        for (Eigen::Index u = 0; u < nu; ++u) {
            C(u, 0) = 1.0;
            for (std::size_t j = 0; j < ncov_eff; ++j) C(u, 1 + static_cast<Eigen::Index>(j)) = Z[static_cast<std::size_t>(u)][j];
            yf(u) = y[static_cast<std::size_t>(u)];
        }
        Eigen::MatrixXd G(nu, static_cast<Eigen::Index>(rows.size()));
        for (std::size_t i = 0; i < rows.size(); ++i) G.col(static_cast<Eigen::Index>(i)) = var_dose[i];
        if (is_lmm) {
            const Eigen::ArrayXd w = (lmm.d.array() + lmm.delta).rsqrt();
            C = w.matrix().asDiagonal() * (lmm.U.transpose() * C);
            yf = w.matrix().asDiagonal() * (lmm.U.transpose() * yf);
            G = w.matrix().asDiagonal() * (lmm.U.transpose() * G);
        }
        const Eigen::HouseholderQR<Eigen::MatrixXd> qr(C);
        const Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(nu, C.cols());
        yf -= Q * (Q.transpose() * yf);
        G -= Q * (Q.transpose() * G);
        Eigen::MatrixXd X(nu, G.cols());
        Eigen::Index k = 0;
        for (Eigen::Index i = 0; i < G.cols(); ++i) {
            const double ss = G.col(i).squaredNorm();
            if (!(ss > 1e-12 * static_cast<double>(nu))) continue;   // no variation left by the covariates
            X.col(k++) = G.col(i) / std::sqrt(ss / static_cast<double>(nu - 1));
            fm_rows.push_back(static_cast<std::size_t>(i));
        }
        X.conservativeResize(nu, k);
        fm = susie(X, yf, opt.finemap_l, opt.cs_coverage);
        for (std::size_t c = 0; c < fm_rows.size(); ++c) pip[fm_rows[c]] = fm.pip[c];
        // Rows map back to the matrix; sets are numbered by their most probable row's strength.
        for (CredibleSet& cs : fm.sets)
            for (std::size_t& r : cs.rows) r = fm_rows[r];
        std::sort(fm.sets.begin(), fm.sets.end(), [&](const CredibleSet& a, const CredibleSet& b) {
            return stronger(a.rows.front(), b.rows.front());
        });
        for (std::size_t c = 0; c < fm.sets.size(); ++c)
            for (std::size_t r : fm.sets[c].rows) if (cs_of[r] < 0) cs_of[r] = static_cast<int>(c + 1);
        if (!fm.converged)
            std::cerr << cli::log_prefix("associate") << "warning: fine-mapping did not converge in "
                      << fm.iterations << " iterations; pip and credible sets are approximate\n";
    }

    // ---- rare-variant groups ----
    // Each group's rare records, weighted by the Beta(1, 25) density of their minor_freq (SKAT's default),
    // missing dosages filled with the record's mean. Scores against the null model's residuals: least
    // squares (for the LMM in the rotated coordinates), or the logistic null fit for case/control.
    struct GroupOut {
        std::string name; std::vector<std::size_t> rows; std::size_t carriers = 0; GroupResult r;
        double skato_p_given = std::numeric_limits<double>::quiet_NaN();
        double acat_p = std::numeric_limits<double>::quiet_NaN();
    };
    std::vector<GroupOut> groups;
    std::size_t n_group_signals = 0;
    if (!opt.group.empty() && !rare.empty()) {
        std::map<std::string, std::vector<std::size_t>> by_group;
        for (std::size_t i = 0; i < rare.size(); ++i) {
            if (opt.group == "all") { by_group["all"].push_back(i); continue; }
            const std::string& key = opt.group == "gene" ? rare[i].gene : opt.group == "layer" ? rare[i].layer : rare[i].bubbles;
            std::set<std::string> names;
            for (const std::string& part : split(key, ';'))
                for (const std::string& t : split(part, ',')) { const std::string u = trim(t); if (!u.empty() && u != ".") names.insert(u); }
            for (const std::string& g : names) by_group[g].push_back(i);
        }
        const auto nu = static_cast<Eigen::Index>(n_used);
        Eigen::MatrixXd C(nu, 1 + static_cast<Eigen::Index>(ncov_eff));
        Eigen::VectorXd yv(nu);
        for (Eigen::Index u = 0; u < nu; ++u) {
            C(u, 0) = 1.0;
            for (std::size_t j = 0; j < ncov_eff; ++j) C(u, 1 + static_cast<Eigen::Index>(j)) = Z[static_cast<std::size_t>(u)][j];
            yv(u) = y[static_cast<std::size_t>(u)];
        }
        Eigen::ArrayXd lw;   // LMM: the rotation's scale
        if (is_lmm) {
            lw = (lmm.d.array() + lmm.delta).rsqrt();
            C = lw.matrix().asDiagonal() * (lmm.U.transpose() * C);
            yv = lw.matrix().asDiagonal() * (lmm.U.transpose() * yv);
        }
        // The null model on covariates Cn: residuals and the scores' covariance factor Zs for a design G.
        struct NullFit { Eigen::MatrixXd C, Q, CWCi; Eigen::VectorXd resid, sw, mu; double sigma = 1.0; bool ok = true; };
        const bool logit = model == "logistic";
        auto null_fit = [&](const Eigen::MatrixXd& Cn) {
            NullFit nf;
            nf.C = Cn;
            // A design without full column rank (a signal the covariates reproduce) or a logistic fit that
            // does not converge (a signal separating cases from controls) has no null model: no test.
            const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Cn);
            if (qr.rank() < Cn.cols()) { nf.ok = false; return nf; }
            if (logit) {
                Eigen::VectorXd beta = Eigen::VectorXd::Zero(Cn.cols()), mu(nu);
                bool converged = false;
                for (int it = 0; it < 100 && !converged; ++it) {
                    mu = (1.0 / (1.0 + (-(Cn * beta).array()).exp())).matrix();
                    const Eigen::VectorXd w = (mu.array() * (1.0 - mu.array())).matrix();
                    const Eigen::LDLT<Eigen::MatrixXd> h(Cn.transpose() * w.asDiagonal() * Cn);
                    if (h.info() != Eigen::Success || !(h.vectorD().minCoeff() > 1e-10 * h.vectorD().maxCoeff())) break;
                    const Eigen::VectorXd step = h.solve(Cn.transpose() * (yv - mu));
                    if (!step.allFinite()) break;
                    beta += step;
                    converged = step.cwiseAbs().maxCoeff() < 1e-8 * (1.0 + beta.cwiseAbs().maxCoeff());
                }
                mu = (1.0 / (1.0 + (-(Cn * beta).array()).exp())).matrix();
                if (!converged || !(mu.minCoeff() > 1e-10 && mu.maxCoeff() < 1.0 - 1e-10)) { nf.ok = false; return nf; }
                nf.resid = yv - mu;
                nf.mu = mu;
                nf.sw = (mu.array() * (1.0 - mu.array())).sqrt().matrix();
                nf.CWCi = (Cn.transpose() * nf.sw.array().square().matrix().asDiagonal() * Cn).inverse();
            } else {
                nf.Q = qr.householderQ() * Eigen::MatrixXd::Identity(nu, Cn.cols());
                nf.resid = yv - nf.Q * (nf.Q.transpose() * yv);
                nf.sigma = std::sqrt(nf.resid.squaredNorm() / static_cast<double>(nu - Cn.cols()));
            }
            return nf;
        };
        auto test_with = [&](const NullFit& nf, const Eigen::MatrixXd& Gz) {
            Eigen::MatrixXd Zs, scaled;
            if (logit) {
                scaled = nf.sw.asDiagonal() * Gz;
                const Eigen::MatrixXd SC = nf.sw.asDiagonal() * nf.C;
                Zs = scaled - SC * (nf.CWCi * (SC.transpose() * scaled));
            } else {
                scaled = nf.sigma * Gz;
                Zs = nf.sigma * (Gz - nf.Q * (nf.Q.transpose() * Gz));
            }
            // A record the null model already contains (a selected signal in its own group) projects to
            // rounding residue: it is dropped, judged against its own size before projection; a group left
            // with nothing is untestable (NA), never significant.
            std::vector<Eigen::Index> live;
            for (Eigen::Index j = 0; j < Zs.cols(); ++j)
                if (Zs.col(j).norm() > 1e-8 * scaled.col(j).norm()) live.push_back(j);
            if (live.empty()) return GroupResult{};
            Eigen::MatrixXd Zl(Zs.rows(), static_cast<Eigen::Index>(live.size()));
            Eigen::MatrixXd Gl(Gz.rows(), static_cast<Eigen::Index>(live.size()));
            for (std::size_t k = 0; k < live.size(); ++k) {
                Zl.col(static_cast<Eigen::Index>(k)) = Zs.col(live[k]);
                Gl.col(static_cast<Eigen::Index>(k)) = Gz.col(live[k]);
            }
            const Eigen::VectorXd S = Gl.transpose() * nf.resid;
            Eigen::MatrixXd Zc = Zl;   // the scores' covariance factor, calibrated for case/control below
            if (logit) {
                // Case/control: the scores' normal tails are wrong for rare records, badly so when cases are
                // few. As SAIGE-GENE+ does, each record's score variance past |z| = 2 is set so that its normal
                // tail gives its saddlepoint p, then the whole group is scaled so that the burden's does too;
                // burden, SKAT and SKAT-O then all use the calibrated variances.
                const std::vector<double> mu(nf.mu.data(), nf.mu.data() + nf.mu.size());
                const std::vector<double> yy(yv.data(), yv.data() + yv.size());
                auto spa_var = [&](double sc, const Eigen::VectorXd& zcol, double v) {
                    // zcol = sw * adjusted genotype; the score is gt'(y - mu), the same as G'(y - mu). As in
                    // SAIGE, only past |z| = 2: near the centre the normal tail is right and the saddlepoint's
                    // is unstable, and a score that is rounding residue must not set the scale.
                    if (!(v > 0.0) || !(sc * sc > 4.0 * v)) return v;
                    std::vector<double> gt(static_cast<std::size_t>(zcol.size()));
                    double gmax = 0.0;
                    for (Eigen::Index i = 0; i < zcol.size(); ++i) {
                        gt[static_cast<std::size_t>(i)] = zcol(i) / nf.sw(i);
                        gmax = std::max(gmax, std::fabs(gt[static_cast<std::size_t>(i)]));
                    }
                    // Scaled to a largest |gt| of 1, the score with it, as for one record: the same tail, with
                    // the saddlepoint's checks independent of the dosage's units.
                    if (!(gmax > 0.0)) return v;
                    for (double& x : gt) x /= gmax;
                    const double ps = spa_two_sided(sc / gmax, gt, mu, &yy, nullptr);
                    if (!(ps > 0.0) || !(ps < 1.0)) return v;
                    const double q = chisq1_upper_quantile(ps);
                    return q > 0.0 ? sc * sc / q : v;
                };
                for (Eigen::Index j = 0; j < Zc.cols(); ++j) {
                    const double v = Zc.col(j).squaredNorm(), va = spa_var(S(j), Zl.col(j), v);
                    if (va != v) Zc.col(j) *= std::sqrt(va / v);
                }
                const Eigen::VectorXd zb = Zc.rowwise().sum();
                const double vb = zb.squaredNorm(), sb = S.sum();
                // The burden's own adjusted genotype is the sum of the records' unscaled ones.
                const Eigen::VectorXd zb0 = Zl.rowwise().sum();
                // A burden the null model already contains (the records summing to a selected signal) projects
                // to rounding residue: its score and variance are both noise, so it sets no scale.
                double col_norms = 0.0;
                for (Eigen::Index j = 0; j < Zl.cols(); ++j) col_norms += Zl.col(j).norm();
                const bool burden_live = zb0.norm() > 1e-6 * col_norms;
                const double vb_target = burden_live ? spa_var(sb, zb0, zb0.squaredNorm()) : zb0.squaredNorm();
                if (vb > 0.0 && vb_target > 0.0 && vb_target != zb0.squaredNorm()) Zc *= std::sqrt(vb_target / vb);
            }
            const GroupResult asym = group_tests(S, Zc);
            if (logit && opt.exact_carriers > 0) {
                std::size_t nc = 0;
                for (Eigen::Index i = 0; i < Gl.rows(); ++i) if (Gl.row(i).cwiseAbs().maxCoeff() > 0.0) ++nc;
                const bool io = nf.C.cols() == 1;
                // As for one record: with covariates, only while the carriers are a small share.
                if (nc <= opt.exact_carriers && (io || static_cast<double>(nc) <= 0.05 * static_cast<double>(Gl.rows())))
                    return group_exact_binary(Gl, Zl, yv, nf.mu, asym, io);
            }
            return asym;
        };
        const NullFit base = null_fit(C);
        // Given the conditional selection's signals, entered as covariates as in the conditional tests: do
        // the rare records add anything beyond them?
        std::vector<std::size_t> signals;
        for (std::size_t i = 0; i < rows.size(); ++i) if (cond_role[i] == "signal") signals.push_back(i);
        NullFit given;
        if (!signals.empty()) {
            // Each signal, and where it is missing when that varies (as in the conditional tests); rotated
            // with the rest under the LMM.
            std::vector<Eigen::VectorXd> extra;
            for (std::size_t sg : signals) {
                extra.push_back(var_dose[sg]);
                Eigen::VectorXd m(nu);
                for (Eigen::Index u = 0; u < nu; ++u) m(u) = var_obs[sg][static_cast<std::size_t>(u)] ? 0.0 : 1.0;
                if (m.maxCoeff() > m.minCoeff()) extra.push_back(m);
            }
            Eigen::MatrixXd Es(nu, static_cast<Eigen::Index>(extra.size()));
            for (std::size_t k = 0; k < extra.size(); ++k)
                Es.col(static_cast<Eigen::Index>(k)) =
                    is_lmm ? Eigen::VectorXd(lw.matrix().asDiagonal() * (lmm.U.transpose() * extra[k])) : extra[k];
            // Columns repeating the covariates or one another are left out, as in the conditional tests.
            const std::vector<char> kept = spanning_columns(C, Es);
            Eigen::MatrixXd Cs(nu, C.cols() + static_cast<Eigen::Index>(std::count(kept.begin(), kept.end(), 1)));
            Cs.leftCols(C.cols()) = C;
            Eigen::Index col = C.cols();
            for (std::size_t k = 0; k < extra.size(); ++k)
                if (kept[k]) Cs.col(col++) = Es.col(static_cast<Eigen::Index>(k));
            given = null_fit(Cs);
        }
        n_group_signals = signals.size();
        if (!base.ok)
            std::cerr << cli::log_prefix("associate") << "warning: the group tests' null model cannot be fitted "
                         "(covariates collinear or separating); group p-values are NA\n";
        if (!signals.empty() && !given.ok)
            std::cerr << cli::log_prefix("associate") << "warning: the selected signals with the covariates admit no "
                         "null model (they separate cases from controls); skato_p_conditional is NA\n";
        for (const auto& [name, members] : by_group) {
            const auto m = static_cast<Eigen::Index>(members.size());
            Eigen::MatrixXd Gz(nu, m);
            std::vector<char> carrier(n_used, 0);
            for (Eigen::Index k = 0; k < m; ++k) {
                const RareRow& rr = rare[members[static_cast<std::size_t>(k)]];
                double sum = 0.0; std::size_t nf = 0;
                for (Eigen::Index u = 0; u < nu; ++u)
                    if (std::isfinite(rr.d(u))) { sum += rr.d(u); ++nf; if (rr.d(u) > 0.0) carrier[static_cast<std::size_t>(u)] = 1; }
                const double mean = nf ? sum / static_cast<double>(nf) : 0.0;
                const double w = 25.0 * std::pow(1.0 - rr.minor_freq, 24.0);   // dbeta(f, 1, 25)
                for (Eigen::Index u = 0; u < nu; ++u) Gz(u, k) = w * (std::isfinite(rr.d(u)) ? rr.d(u) : mean);
            }
            if (is_lmm) Gz = lw.matrix().asDiagonal() * (lmm.U.transpose() * Gz);
            GroupOut go;
            go.name = name;
            go.rows = members;
            go.carriers = static_cast<std::size_t>(std::count(carrier.begin(), carrier.end(), 1));
            if (base.ok) {
                go.r = test_with(base, Gz);
                // ACAT-O: burden, SKAT and SKAT-O with the Beta(1, 25) weights and with equal ones, combined
                // by the Cauchy rule (Liu et al. 2019), valid whatever their correlation.
                Eigen::MatrixXd G1 = Gz;
                for (Eigen::Index k = 0; k < m; ++k) {
                    const double w = 25.0 * std::pow(1.0 - rare[members[static_cast<std::size_t>(k)]].minor_freq, 24.0);
                    G1.col(k) /= w;
                }
                const GroupResult r1 = test_with(base, G1);
                double t = 0.0; int nt = 0;
                for (double pv : {go.r.burden_p, go.r.skat_p, go.r.skato_p, r1.burden_p, r1.skat_p, r1.skato_p}) {
                    if (!std::isfinite(pv)) continue;
                    // A p of 1 (common for exact tests) would dominate the sum: capped at 0.999, as is usual.
                    const double pc = std::min(std::max(pv, 1e-300), 0.999);
                    t += pc < 1e-15 ? 1.0 / (pc * M_PI) : std::tan((0.5 - pc) * M_PI);
                    ++nt;
                }
                if (nt > 0) {
                    t /= nt;
                    go.acat_p = t > 1e15 ? 1.0 / (t * M_PI) : 0.5 - std::atan(t) / M_PI;
                }
            }
            if (!signals.empty() && given.ok) go.skato_p_given = test_with(given, Gz).skato_p;
            groups.push_back(std::move(go));
        }
    }

    const double bonf_threshold = n_tests > 0 ? 0.05 / static_cast<double>(n_tests) : kNaN;
    const double bonf_threshold_meff = meff > 0 ? 0.05 / static_cast<double>(meff) : kNaN;
    std::size_t n_sig_bonf = 0, n_sig_fdr = 0, n_sig_meff = 0, n_sig_perm = 0;
    for (std::size_t i = 0; i < n_tests; ++i) {
        if (std::isfinite(rows[i].p) && rows[i].p < bonf_threshold) ++n_sig_bonf;
        if (std::isfinite(rows[i].p) && rows[i].p < bonf_threshold_meff) ++n_sig_meff;
        if (std::isfinite(qv[i]) && qv[i] < 0.05) ++n_sig_fdr;
        if (std::isfinite(p_perm[i]) && p_perm[i] < 0.05) ++n_sig_perm;
    }
    const double lambda_gc = genomic_lambda(pv);
    // lambda_GC reads the median chi-square assuming most tests are null and independent. At one locus a
    // signal and everything correlated with it can be most of the tests, and k-mers, nodes and edges of a
    // few events are not independent: it is labelled an inflation estimate only in variant mode, with
    // enough records, most of them null.
    const bool lambda_meaningful = variant_mode &&
        n_tests >= 100 && static_cast<double>(n_sig_fdr) < 0.25 * static_cast<double>(n_tests);

    // ---- write assoc.tsv (sorted by p) ----
    // Staged, committed together, and never over an input.
    const std::string assoc_final = opt.out_prefix + ".assoc.tsv";
    const std::string summary_final = opt.out_prefix + ".summary.tsv";
    const std::string layers_final = opt.out_prefix + ".layers.tsv";
    const std::string cs_final = opt.out_prefix + ".cs.tsv";
    const std::string groups_final = opt.out_prefix + ".groups.tsv";
    std::vector<std::string> inputs = {opt.node_genes, opt.phenotype, opt.kinship, opt.haplotypes};
    for (const auto* list : {&opt.genotypes, &opt.samples, &opt.feature_annot})
        inputs.insert(inputs.end(), list->begin(), list->end());
    std::vector<std::string> outputs = {assoc_final, summary_final, layers_final, cs_final};
    if (!opt.group.empty()) outputs.push_back(groups_final);
    cli::reject_output_collisions("associate", outputs, inputs);
    cli::StagedOutputs staged("associate");
    const std::string assoc_staged = staged.stage(assoc_final);
    const std::string summary_staged = staged.stage(summary_final);
    const std::string layers_staged = staged.stage(layers_final);
    const std::string cs_staged = staged.stage(cs_final);
    const std::string groups_staged = opt.group.empty() ? std::string() : staged.stage(groups_final);

    std::vector<std::size_t> order(n_tests);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), stronger);
    auto fmt = [](double v) {
        if (!std::isfinite(v)) return std::string("NA");
        std::ostringstream o; o << v; return o.str();
    };
    const std::string effect = (model == "logistic") ? "log_or" : "beta";
    auto flag = [](int v) { return v < 0 ? std::string(".") : std::to_string(v); };
    {
        std::ofstream out(assoc_staged);
        if (!out) throw std::runtime_error("cannot write " + assoc_final);
        out << "feature_id\tlayer\tbubbles\tnodes\tchrom\tpos\tend\tn\tn_conditional\tminor_freq\t" << effect
            << "\tse\tstat\tp\tp_method\teffect_status\tmac_case\tmac_ctrl\tp_bonf\tp_bonf_meff\tq_bh\tp_perm\taf\tan\tlow_af\tclump\tis_lead\tgene"
               "\tp_conditional\tcond_role\tpip\tcs\tp_genotypic\tp_nonadditive\n";
        for (std::size_t k = 0; k < n_tests; ++k) {
            const std::size_t i = order[k];
            const Row& r = rows[i];
            const double p_bonf = std::isfinite(r.p) ? std::min(1.0, r.p * static_cast<double>(n_tests)) : kNaN;
            const double p_bonf_meff = std::isfinite(r.p) && meff > 0 ? std::min(1.0, r.p * static_cast<double>(meff)) : kNaN;
            out << r.id << '\t' << r.layer << '\t' << r.bubbles << '\t' << r.nodes << '\t'
                << r.chrom << '\t' << r.pos << '\t' << r.end << '\t'
                // A conditional n only beside a conditional p.
                << r.n << '\t' << (n_cond[i] && std::isfinite(p_cond[i]) ? std::to_string(n_cond[i]) : std::string("."))
                << '\t' << fmt(r.minor_freq) << '\t' << fmt(r.beta) << '\t' << fmt(r.se)
                << '\t' << fmt(r.z) << '\t' << fmt(r.p) << '\t' << r.p_method << '\t' << r.effect_status
                << '\t' << (r.mac_case < 0 ? std::string(".") : std::to_string(r.mac_case))
                << '\t' << (r.mac_ctrl < 0 ? std::string(".") : std::to_string(r.mac_ctrl))
                << '\t' << fmt(p_bonf) << '\t' << (meff > 0 ? fmt(p_bonf_meff) : std::string("."))
                << '\t' << fmt(qv[i]) << '\t' << (perm_min.empty() ? std::string(".") : fmt(p_perm[i]))
                << '\t' << r.af << '\t' << r.an << '\t' << flag(r.low_af)
                << '\t' << flag(r.clump) << '\t' << flag(r.is_lead) << '\t' << r.gene
                << '\t' << fmt(p_cond[i]) << '\t' << cond_role[i]
                << '\t' << (finemap ? fmt(pip[i]) : std::string(".")) << '\t' << flag(cs_of[i])
                << '\t' << fmt(r.p_geno) << '\t' << fmt(r.p_nonadd) << '\n';
        }
        // Closed and checked: a destructor would not report a failed write.
        out.close();
        if (!out) throw std::runtime_error("failed writing " + assoc_final);
    }
    // ---- summary (also the per-plot threshold lines) ----
    {
        std::ofstream out(summary_staged);
        if (!out) throw std::runtime_error("cannot write " + summary_final);
        out << "key\tvalue\n"
            << "model\t" << model << '\n'
            << "phenotype_type\t" << (binary ? "binary" : "quantitative") << '\n'
            << "covariates\t" << ncov_eff << '\n'
            << "pca_covariates\t" << kinship_pcs.cols() << '\n'
            << "samples_used\t" << n_used << '\n'
            << "genotype_rows\t" << n_geno_rows << '\n'
            << "features_tested\t" << n_tests << '\n'
            << "unit\t" << (variant_mode ? "variant" : "feature") << '\n'
            << "meff\t" << (meff > 0 ? std::to_string(meff) : std::string("NA")) << '\n'
            << "meff_method\t" << (meff_eigen > 0 ? "eigenvalue (Li-Ji, phenotype-blind)"
                                                   : (variant_mode ? "ld_clumping (heuristic)" : "none (feature unit)"))
               << '\n'
            << "meff_eigen\t" << (meff_eigen > 0 ? std::to_string(meff_eigen) : std::string("NA")) << '\n'
            << "meff_ld_clumping\t" << (meff_clump > 0 ? std::to_string(meff_clump) : std::string("NA")) << '\n'
            << "dropped_min_maf\t" << n_dropped_maf << '\n'
            << "dropped_fit\t" << n_dropped_fit << '\n'
            << "score_normal_fallback\t" << n_normal_fallback << '\n'
            << "bonferroni_threshold\t" << fmt(bonf_threshold) << '\n'
            << "bonferroni_threshold_meff\t" << fmt(bonf_threshold_meff) << '\n'
            << "nominal_threshold\t0.05\n"
            << "significant_bonferroni\t" << n_sig_bonf << '\n'
            << "significant_bonferroni_meff\t" << (meff > 0 ? std::to_string(n_sig_meff) : std::string("NA")) << '\n'
            << "significant_fdr05\t" << n_sig_fdr << '\n'
            << "permutations\t" << perm_min.size() << '\n'
            << "permutation_threshold\t" << fmt(perm_threshold) << '\n'
            << "significant_permutation\t" << (perm_min.empty() ? std::string("NA") : std::to_string(n_sig_perm)) << '\n'
            << "lambda_gc\t" << fmt(lambda_gc) << '\n';
        if (variant_mode) {
            out << "cojo_entry_threshold\t" << fmt(cojo_entry) << '\n'
                << "cojo_entry_source\t" << cojo_entry_source << '\n'
                << "cojo_independent_signals\t" << cojo_n_signals << '\n';
        }
        if (!opt.group.empty())
            out << "group_unit\t" << opt.group << '\n'
                << "rare_maf\t" << fmt(opt.rare_maf) << '\n'
                << "rare_rows\t" << rare.size() << '\n'
                << "groups_tested\t" << groups.size() << '\n'
                << "groups_conditioned_on\t" << n_group_signals << '\n';
        if (variant_mode)
            out << "finemap_effects\t" << opt.finemap_l << '\n'
                << "credible_sets\t" << (finemap ? std::to_string(fm.sets.size()) : std::string("NA")) << '\n'
                << "finemap_converged\t" << (finemap ? (fm.converged ? "1" : "0") : "NA") << '\n';
        if (model == "lmm")
            out << "lmm_delta\t" << fmt(lmm.delta) << '\n'
                << "lmm_h2\t" << fmt(lmm_h2) << '\n'
                << "lmm_delta_at_boundary\t" << (lmm.at_boundary ? 1 : 0) << '\n';
        out.close();
        if (!out) throw std::runtime_error("failed writing " + summary_final);
    }

    // ---- layers: each layer's strongest row, side by side ----
    {
        std::ofstream out(layers_staged);
        if (!out) throw std::runtime_error("cannot write " + layers_final);
        std::vector<const std::pair<const std::string, LayerBest>*> by_p;
        for (const auto& kv : layer_best) if (kv.second.found) by_p.push_back(&kv);
        std::sort(by_p.begin(), by_p.end(), [&](auto* a, auto* b) { return stronger(a->second.row, b->second.row); });
        const Eigen::VectorXd top = by_p.empty() ? Eigen::VectorXd() : by_p.front()->second.dose;
        const std::vector<char> top_obs = by_p.empty() ? std::vector<char>() : by_p.front()->second.obs;
        // The two rows against each other, on the samples both measure (every sample, mean-imputed, under
        // the LMM): the p of one given the other, either way. A source whose best row is a near-copy of the
        // top row cannot be told apart from it in either direction.
        auto given = [&](const Eigen::VectorXd& a, const std::vector<char>& ao,
                         const Eigen::VectorXd& b, const std::vector<char>& bo) -> double {
            std::vector<std::size_t> keep;
            for (std::size_t u = 0; u < n_used; ++u) if (is_lmm || (ao[u] && bo[u])) keep.push_back(u);
            std::size_t nk = keep.size(), pc = p_dim + 1;
            if (nk <= pc) return kNaN;
            std::vector<double> X(nk * pc), yk(nk);
            Eigen::VectorXd ra, rb;
            if (is_lmm) { ra = rotate(a); rb = rotate(b); }
            for (std::size_t k = 0; k < nk; ++k) {
                const std::size_t u = keep[k];
                const auto ui = static_cast<Eigen::Index>(u);
                std::size_t col = 0;
                if (is_lmm) {
                    yk[k] = lmm_y(ui);
                    X[k * pc + col++] = lmm_C(ui, 0);
                    X[k * pc + col++] = ra(ui);
                    for (std::size_t j = 0; j < ncov_eff; ++j) X[k * pc + col++] = lmm_C(ui, 1 + static_cast<Eigen::Index>(j));
                    X[k * pc + col++] = rb(ui);
                } else {
                    yk[k] = y[u];
                    X[k * pc + col++] = 1.0;
                    X[k * pc + col++] = a(ui);
                    for (std::size_t j = 0; j < ncov_eff; ++j) X[k * pc + col++] = Z[u][j];
                    X[k * pc + col++] = b(ui);
                }
            }
            std::size_t pd = p_dim;
            if (model == "logistic" && !drop_separated(X, yk, nk, pc, pd)) return kNaN;   // as in the conditional tests
            const FitResult f = (model == "logistic") ? score_logistic(X, yk, nk, pc, opt.exact_carriers)
                                                      : fit_linear(X, yk, nk, pc);
            return f.ok ? f.p : kNaN;
        };
        // Linear model, variant unit: how much more of the phenotype (beyond the covariates) this layer's best
        // row explains than the top layer's (negative when the top layer's explains more). Each best row was picked on the phenotype, so its
        // in-sample share carries a winner's bias that grows with the layer's number of rows. The difference
        // is therefore taken out of sample: each of 1,000 bootstrap resamples picks each layer's best row
        // on the resampled samples, and scores the two on the samples the resample left out. delta_r2 is the
        // median, the interval the 2.5% and 97.5% points; an interval spanning 0 means the two layers cannot
        // be ranked. Dosages are mean-imputed, as in fine-mapping.
        const std::string top_layer = by_p.empty() ? std::string() : by_p.front()->first;
        Eigen::MatrixXd Rres;          // residuals on the covariates: column 0 the phenotype, then rows
        std::vector<Eigen::Index> col_of(rows.size(), -1);
        const bool oob = model == "linear" && variant_mode && !by_p.empty() && n_used > p_dim + 2;
        if (oob) {
            std::vector<std::size_t> cand;
            for (std::size_t i = 0; i < n_tests; ++i) if (std::isfinite(rows[i].p)) cand.push_back(i);
            const auto nu = static_cast<Eigen::Index>(n_used);
            Eigen::MatrixXd C(nu, 1 + static_cast<Eigen::Index>(ncov_eff));
            Rres.resize(nu, 1 + static_cast<Eigen::Index>(cand.size()));
            for (Eigen::Index u = 0; u < nu; ++u) {
                C(u, 0) = 1.0;
                for (std::size_t j = 0; j < ncov_eff; ++j) C(u, 1 + static_cast<Eigen::Index>(j)) = Z[static_cast<std::size_t>(u)][j];
                Rres(u, 0) = y[static_cast<std::size_t>(u)];
            }
            for (std::size_t k = 0; k < cand.size(); ++k) {
                Rres.col(1 + static_cast<Eigen::Index>(k)) = var_dose[cand[k]];
                col_of[cand[k]] = 1 + static_cast<Eigen::Index>(k);
            }
            const Eigen::HouseholderQR<Eigen::MatrixXd> qr(C);
            const Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(nu, C.cols());
            Rres -= Q * (Q.transpose() * Rres);
        }
        auto delta_r2 = [&](const std::string& layer, double ci[2]) -> double {
            ci[0] = ci[1] = kNaN;
            if (!oob) return kNaN;
            std::vector<Eigen::Index> ca, ct;   // the two layers' candidate columns
            for (std::size_t i = 0; i < n_tests; ++i) {
                if (col_of[i] < 0) continue;
                if (rows[i].layer == layer) ca.push_back(col_of[i]);
                else if (rows[i].layer == top_layer) ct.push_back(col_of[i]);
            }
            if (ca.empty() || ct.empty()) return kNaN;
            const Eigen::Index nu = Rres.rows();
            Eigen::MatrixXd A(nu, static_cast<Eigen::Index>(ca.size() + ct.size()));
            for (std::size_t k = 0; k < ca.size(); ++k) A.col(static_cast<Eigen::Index>(k)) = Rres.col(ca[k]);
            for (std::size_t k = 0; k < ct.size(); ++k) A.col(static_cast<Eigen::Index>(ca.size() + k)) = Rres.col(ct[k]);
            const Eigen::VectorXd ry = Rres.col(0);
            const Eigen::MatrixXd A2 = A.array().square().matrix();
            const auto na = static_cast<Eigen::Index>(ca.size());
            std::mt19937_64 rng(opt.seed);
            std::vector<double> ds;
            Eigen::VectorXd w(nu);
            for (int b = 0; b < 1000; ++b) {
                w.setZero();
                for (Eigen::Index k = 0; k < nu; ++k) w(static_cast<Eigen::Index>(rng() % static_cast<std::uint64_t>(nu))) += 1.0;
                const Eigen::VectorXd wy = w.cwiseProduct(ry);
                const double syy = wy.dot(ry);
                if (!(syy > 0.0)) continue;
                const Eigen::VectorXd sxy = A.transpose() * wy, sxx = A2.transpose() * w;
                auto best = [&](Eigen::Index from, Eigen::Index to) {
                    Eigen::Index arg = -1; double top_r2 = -1.0;
                    for (Eigen::Index c = from; c < to; ++c)
                        if (sxx(c) > 0.0) { const double r2 = sxy(c) * sxy(c) / (sxx(c) * syy); if (r2 > top_r2) { top_r2 = r2; arg = c; } }
                    return arg;
                };
                const Eigen::Index ba = best(0, na), bt = best(na, A.cols());
                if (ba < 0 || bt < 0) continue;
                // Scored on the samples this resample left out (about a third).
                double oyy = 0, oaa = 0, ott = 0, oya = 0, oyt = 0; std::size_t no = 0;
                for (Eigen::Index u = 0; u < nu; ++u) {
                    if (w(u) != 0.0) continue;
                    const double yv = ry(u), av = A(u, ba), tv = A(u, bt);
                    oyy += yv * yv; oaa += av * av; ott += tv * tv; oya += yv * av; oyt += yv * tv; ++no;
                }
                if (no < 3 || !(oyy > 0.0) || !(oaa > 0.0) || !(ott > 0.0)) continue;
                ds.push_back(oya * oya / (oyy * oaa) - oyt * oyt / (oyy * ott));
            }
            if (ds.size() < 100) return kNaN;
            std::sort(ds.begin(), ds.end());
            const auto at = [&](double q) { return ds[static_cast<std::size_t>(q * static_cast<double>(ds.size() - 1))]; };
            ci[0] = at(0.025); ci[1] = at(0.975);
            return at(0.5);
        };
        out << "layer\trows_tested\tsignificant_bonferroni\tsignificant_bonferroni_meff\tsignificant_permutation"
               "\tsignificant_fdr05\tbest_feature\tbest_" << effect << "\tbest_p\tbest_r2\tr2_with_top"
               "\tbest_p_conditional\tbest_cond_role\tcs_prob\tp_best_given_top\tp_top_given_best"
               "\tdelta_r2\tdelta_r2_low\tdelta_r2_high\n";
        for (const auto* kv : by_p) {
            const std::string& layer = kv->first;
            std::size_t tested = 0, sb = 0, sm = 0, sp = 0, sf = 0;
            // The layer's probability within the credible sets: how many located signals it holds.
            double cs_prob = 0.0;
            for (const CredibleSet& cs : fm.sets)
                for (std::size_t k = 0; k < cs.rows.size(); ++k)
                    if (rows[cs.rows[k]].layer == layer) cs_prob += cs.prob[k];
            for (std::size_t i = 0; i < n_tests; ++i) {
                if (rows[i].layer != layer) continue;
                ++tested;
                if (std::isfinite(rows[i].p) && rows[i].p < bonf_threshold) ++sb;
                if (std::isfinite(rows[i].p) && rows[i].p < bonf_threshold_meff) ++sm;
                if (std::isfinite(p_perm[i]) && p_perm[i] < 0.05) ++sp;
                if (std::isfinite(qv[i]) && qv[i] < 0.05) ++sf;
            }
            const std::size_t i = kv->second.row;
            const Row& r = rows[i];
            // The share of the phenotype's variance left by the covariates that the row explains, from
            // its t statistic. None for the logistic score test, nor for the LMM, whose t lives in the
            // kinship-rotated coordinates, not on the phenotype's scale.
            double r2 = kNaN;
            if (model == "linear") {
                const double df = static_cast<double>(r.n) - static_cast<double>(p_dim);
                if (df > 0.0 && std::isfinite(r.z)) r2 = r.z * r.z / (r.z * r.z + df);
            }
            out << layer << '\t' << tested << '\t' << sb << '\t' << (meff > 0 ? std::to_string(sm) : "NA")
                << '\t' << (perm_min.empty() ? std::string("NA") : std::to_string(sp)) << '\t' << sf
                << '\t' << r.id << '\t' << fmt(r.beta) << '\t' << fmt(r.p) << '\t' << fmt(r2)
                << '\t' << fmt(r2_vec(kv->second.dose, top)) << '\t' << fmt(p_cond[i]) << '\t' << cond_role[i]
                << '\t' << (finemap ? fmt(cs_prob) : std::string("NA"));
            if (kv == by_p.front()) {
                out << "\tNA\tNA\tNA\tNA\tNA\n";
            } else {
                double ci[2];
                const double d = delta_r2(layer, ci);
                out << '\t' << fmt(given(kv->second.dose, kv->second.obs, top, top_obs))
                    << '\t' << fmt(given(top, top_obs, kv->second.dose, kv->second.obs))
                    << '\t' << fmt(d) << '\t' << fmt(ci[0]) << '\t' << fmt(ci[1]) << '\n';
            }
        }
        out.close();
        if (!out) throw std::runtime_error("failed writing " + layers_final);
    }

    // ---- credible sets: one line each, with how its probability splits across layers ----
    {
        std::ofstream out(cs_staged);
        if (!out) throw std::runtime_error("cannot write " + cs_final);
        out << "cs\trows\tcoverage\tpurity\ttop_feature\ttop_prob\tlayers\tfeatures\n";
        for (std::size_t c = 0; c < fm.sets.size(); ++c) {
            const CredibleSet& cs = fm.sets[c];
            std::map<std::string, double> by_layer;
            std::string ids;
            for (std::size_t k = 0; k < cs.rows.size(); ++k) {
                by_layer[rows[cs.rows[k]].layer] += cs.prob[k];
                ids += (k ? "," : "") + rows[cs.rows[k]].id;
            }
            std::vector<std::pair<std::string, double>> lv(by_layer.begin(), by_layer.end());
            std::stable_sort(lv.begin(), lv.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
            std::string layers;
            for (const auto& [layer, pr] : lv) layers += (layers.empty() ? "" : ";") + layer + ":" + fmt(pr);
            out << (c + 1) << '\t' << cs.rows.size() << '\t' << fmt(cs.coverage) << '\t' << fmt(cs.purity)
                << '\t' << rows[cs.rows.front()].id << '\t' << fmt(cs.prob.front()) << '\t' << layers
                << '\t' << ids << '\n';
        }
        out.close();
        if (!out) throw std::runtime_error("failed writing " + cs_final);
    }

    // ---- rare-variant groups ----
    if (!opt.group.empty()) {
        std::vector<std::size_t> ord(groups.size());
        std::iota(ord.begin(), ord.end(), 0);
        std::stable_sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) {
            const double pa = groups[a].r.skato_p, pb = groups[b].r.skato_p;
            if (std::isfinite(pa) != std::isfinite(pb)) return std::isfinite(pa);
            return std::isfinite(pa) && pa < pb;
        });
        std::ofstream out(groups_staged);
        if (!out) throw std::runtime_error("cannot write " + groups_final);
        out << "group\trows\tcarriers\tburden_p\tskat_p\tskato_p\tskato_rho\tp_bonf\tskato_p_conditional\tacat_p\tfeatures\n";
        for (std::size_t k : ord) {
            const GroupOut& g = groups[k];
            std::string ids;
            for (std::size_t i : g.rows) ids += (ids.empty() ? "" : ",") + rare[i].id;
            out << g.name << '\t' << g.rows.size() << '\t' << g.carriers << '\t' << fmt(g.r.burden_p) << '\t'
                << fmt(g.r.skat_p) << '\t' << fmt(g.r.skato_p) << '\t' << fmt(g.r.skato_rho) << '\t'
                << fmt(std::isfinite(g.r.skato_p) ? std::min(1.0, g.r.skato_p * static_cast<double>(groups.size())) : kNaN) << '\t'
                << fmt(g.skato_p_given) << '\t' << fmt(g.acat_p) << '\t' << ids << '\n';
        }
        out.close();
        if (!out) throw std::runtime_error("failed writing " + groups_final);
    }

    // The files installed as one transaction.
    staged.commit();

    log.info("model: " + model + " (" + (binary ? "binary" : "quantitative") + "); samples: " +
             std::to_string(n_used) + "; covariates: " + std::to_string(ncov_eff) +
             (kinship_pcs.cols() > 0 ? " (" + std::to_string(kinship_pcs.cols()) + " of them kinship PCs)" : ""));
    log.info("unit: " + std::string(variant_mode ? "variant" : "feature"));
    log.info("rows tested: " + std::to_string(n_tests) + (meff > 0 ? "; Meff: " + std::to_string(meff) : std::string()));
    log.info("rows not tested: below --min-maf: " + std::to_string(n_dropped_maf) + "; not fitted: " +
             std::to_string(n_dropped_fit));
    if (binary) {
        // With few carriers the score's distribution is lumpy, and its tails rest on the saddlepoint or
        // on enumeration; --group pools such records.
        std::size_t n_rare = 0;
        for (const Row& r : rows)
            if (r.mac_case >= 0 && r.mac_case + r.mac_ctrl < 10) ++n_rare;
        if (n_rare > 0)
            log.info("rows with fewer than 10 samples off the modal dosage: " + std::to_string(n_rare) +
                     " (p from the exact or saddlepoint tail: exploratory; --group tests rare records together)");
    }
    if (finemap)
        log.info("fine-mapping: credible sets: " + std::to_string(fm.sets.size()) + "; effects fitted: " +
                 std::to_string(opt.finemap_l));
    if (!opt.group.empty()) {
        std::set<std::size_t> grouped;   // a rare record without a gene or site is in no group
        for (const GroupOut& go : groups) grouped.insert(go.rows.begin(), go.rows.end());
        log.info("groups by " + opt.group + ": " + std::to_string(groups.size()) + "; rare records in them: " +
                 std::to_string(grouped.size()));
    }
    std::string sig = "significant: Bonferroni: " + std::to_string(n_sig_bonf);
    if (meff > 0) sig += "; Bonferroni on Meff: " + std::to_string(n_sig_meff);
    if (!perm_min.empty()) sig += "; permutation: " + std::to_string(n_sig_perm);
    sig += "; FDR<0.05: " + std::to_string(n_sig_fdr);
    log.info(sig);
    if (variant_mode) log.info("signals (conditional selection): " + std::to_string(cojo_n_signals));
    const auto fmt3 = [](double v) {
        if (!std::isfinite(v)) return std::string("NA");
        std::ostringstream o; o << std::setprecision(3) << v; return o.str();
    };
    log.info("lambda_GC: " + fmt3(lambda_gc) + (lambda_meaningful ? std::string()
             : std::string(" (not an inflation estimate here: the tests are few or correlated)")));
    {
        std::vector<std::string> written = {assoc_final, summary_final, layers_final, cs_final};
        if (!opt.group.empty()) written.push_back(groups_final);
        log.wrote(written);
    }
    log.done();
    return 0;
}

}  // namespace panvar
