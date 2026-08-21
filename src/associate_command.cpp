#include "panvar/associate_command.hpp"
#include "panvar/cli_utils.hpp"

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

// Regularized incomplete beta I_x(a,b), by the continued fraction of Numerical Recipes 6.4. Needed for
// the Student-t tail below; the normal tail is not the right reference distribution for a linear model
// whose residual variance was estimated from the same data.
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

// Two-sided Student-t tail on `df` degrees of freedom. Falls back to the normal when df is large enough
// that the difference is below the printing precision anyway.
double student_p(double t, double df) {
    if (!std::isfinite(t) || !(df > 0.0)) return kNaN;
    if (df > 1e6) return wald_p(t);
    const double p = betai(0.5 * df, 0.5, df / (df + t * t));
    return p < 1e-300 ? 1e-300 : (p > 1.0 ? 1.0 : p);
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
    bool spa = false;      // p came from the saddlepoint rather than the normal tail
    bool exact_tail = false;  // p is the exact boundary probability, not any approximation
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
    // Student-t on n - p degrees of freedom, not the normal: sigma^2 is estimated from the same
    // residuals, so the statistic is t-distributed. Immaterial at GWAS n, but panvar is a locus tool and
    // a few dozen samples is a realistic cohort, where the normal tail is anti-conservative.
    r.p = student_p(r.z, static_cast<double>(n - p));
    r.ok = std::isfinite(r.p);
    return r;
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
        double delta = 0.0;
        for (std::size_t a = 0; a < p; ++a) { delta += std::fabs(nb[a] - beta[a]); beta[a] = nb[a]; }
        if (delta < 1e-8) { converged = true; break; }
    }
    // A run that hit the iteration cap has not found the maximum -- under separation it is still
    // diverging. Reporting the last iterate as an estimate would be reporting where we ran out of
    // patience, so the caller counts it as an unfittable feature instead.
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
// Under separation the ordinary maximum likelihood estimate diverges -- there is no finite coefficient
// that maximises the likelihood -- so an MLE fit returns nothing and the effect size is simply absent.
// Firth adds Jeffreys' invariant prior, |I(beta)|^{1/2}, as a penalty. The penalised score is
//
//   U*_j = sum_i x_ij [ y_i - mu_i + h_i (0.5 - mu_i) ],    h = diag( W^{1/2} X (X'WX)^-1 X' W^{1/2} )
//
// which stays finite because the h_i(0.5 - mu_i) term pulls the fitted probabilities back off 0 and 1.
// The estimate is also first-order unbiased, which ordinary ML is not at small counts.
//
// This is used ONLY for the reported effect size when ML fails. The p-value stays the score test's,
// which needs no alternative fit at all.
bool firth_logistic(const std::vector<double>& X, const std::vector<double>& y, std::size_t n,
                    std::size_t p, std::vector<double>& beta, std::vector<double>& inv) {
    beta.assign(p, 0.0);
    inv.clear();
    std::vector<double> mu(n), w(n);
    for (int iter = 0; iter < 200; ++iter) {
        std::vector<double> XtWX(p * p, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double eta = 0.0;
            for (std::size_t a = 0; a < p; ++a) eta += X[i * p + a] * beta[a];
            mu[i] = 1.0 / (1.0 + std::exp(-eta));
            w[i] = std::max(mu[i] * (1.0 - mu[i]), 1e-10);
            for (std::size_t a = 0; a < p; ++a)
                for (std::size_t b = a; b < p; ++b) XtWX[a * p + b] += X[i * p + a] * w[i] * X[i * p + b];
        }
        for (std::size_t a = 0; a < p; ++a) for (std::size_t b = 0; b < a; ++b) XtWX[a * p + b] = XtWX[b * p + a];
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
        double step = 0.0;
        for (std::size_t a = 0; a < p; ++a) {
            double d = 0.0;
            for (std::size_t b = 0; b < p; ++b) d += inv[a * p + b] * Ustar[b];
            // Half-step when the move is large: the penalised likelihood is well behaved but a raw
            // Newton step can still overshoot from a cold start.
            if (d > 5.0) d = 5.0;
            if (d < -5.0) d = -5.0;
            beta[a] += d;
            step += std::fabs(d);
        }
        if (step < 1e-10) return true;
    }
    return false;   // did not settle even with the penalty
}

// Rao score test on the genotype column, evaluated under the NULL fit (covariates only, genotype
// coefficient fixed at zero).
//
// Why not the Wald test here. Wald divides the estimate by its own standard error, and for a rare
// variant in an unbalanced case/control study the fit approaches separation: |beta| grows, its standard
// error grows faster, and the statistic collapses toward zero. Measured on the LPA cohort (492 cases,
// 5213 controls) every feature below minor frequency 0.01 failed a uniformity check under permutation
// -- lambda_GC 0.80, worst KS p 3e-12 -- while every feature above it passed. The score test never fits
// the alternative, so it does not have a standard error to inflate, and it is the standard remedy.
//
// The genotype column is dropped from the design to form the null, which also means each feature gets
// its own null fit -- correct rather than wasteful, since each has its own complete-case sample set.
// ---- Saddlepoint approximation for the binary score statistic (Dey et al. fastSPA; SAIGE) --------
//
// The score S = sum_i Gt_i (y_i - mu_i) is a sum of INDEPENDENT bounded terms, so its cumulant
// generating function is available in closed form -- no approximation needed:
//
//   K(t)   = sum_i [ log(1 - mu_i + mu_i e^{t Gt_i}) - t Gt_i mu_i ]
//   K'(t)  = sum_i [ Gt_i mu_i e^{t Gt_i} / (1 - mu_i + mu_i e^{t Gt_i}) - Gt_i mu_i ]
//   K''(t) = sum_i [ Gt_i^2 mu_i (1-mu_i) e^{t Gt_i} / (1 - mu_i + mu_i e^{t Gt_i})^2 ]
//
// The normal approximation matches only the first two cumulants, which is why it drifts in the far tail
// exactly when the terms are skewed -- a rare variant under case/control imbalance. The saddlepoint
// expands about the point where the tilted distribution is centred on the observed value, so it stays
// accurate out where a regional Bonferroni threshold actually sits.
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

// log(1 - mu + mu e^x), computed as a log-sum-exp so a large saddlepoint cannot overflow exp(x) and
// silently turn the whole SPA into a normal fallback.
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

// Solve K'(zeta) = s. K' is strictly increasing, so bracket then bisect with a Newton step where it is
// safe -- robust matters more than fast here, and it runs once per tail per feature.
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
    if (!std::isfinite(t)) return kNaN;
    return std::min(1.0, std::max(0.0, t));
}

// Two-sided p for the observed score. The score's distribution is NOT symmetric under case/control
// imbalance, so both tails are computed rather than one doubled.
double spa_two_sided(double s, const std::vector<double>& gt, const std::vector<double>& mu,
                     const std::vector<double>* y = nullptr, bool* exact_out = nullptr) {
    // S is a bounded sum. The two tails are evaluated INDEPENDENTLY, because for an asymmetric score
    // distribution they are in different regimes: a threshold can be outside the support on one side
    // (probability exactly zero), sitting on the boundary atom on the other, or interior on both.
    // Treating "at or beyond the boundary" as one case conflated all three -- it charged the boundary
    // atom for a threshold that is strictly outside the support, where the probability is zero, and it
    // skipped the opposite tail entirely when that tail was interior and non-negligible.
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
    // Whether the OBSERVED outcome really is the extreme configuration, checked against y rather than
    // inferred from a floating-point tolerance on the statistic. It applies only to the side the
    // observation is actually ON -- s can be extremal at one end at most, and demanding it of the
    // opposite tail rejected a legitimate boundary atom there. The opposite threshold reaching its own
    // boundary is a numeric fact about the threshold, not a claim about what was observed.
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
    }

    if (!std::isfinite(up) || !std::isfinite(dn)) return kNaN;
    const double p = up + dn;
    if (!std::isfinite(p) || p <= 0.0) return kNaN;
    // Only claim an exact p when BOTH tails were exact or provably zero; a mix is a hybrid and is
    // reported as the saddlepoint, which is what it mostly is.
    if (exact_out != nullptr) *exact_out = (up_exact && dn_exact);
    return std::min(1.0, p);
}

FitResult score_logistic(const std::vector<double>& X, const std::vector<double>& y,
                         std::size_t n, std::size_t p, double spa_cutoff = 2.0) {
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
    r.ok = std::isfinite(r.p);

    // Past the cutoff, replace the normal tail with the saddlepoint. Below it the two agree to well
    // within printing precision and the normal is far cheaper, which is the same gate SAIGE uses.
    if (spa_cutoff > 0.0 && std::fabs(r.z) > spa_cutoff) {
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
        double S = 0.0;
        for (std::size_t i = 0; i < n; ++i) S += gt[i] * (y[i] - mu[i]);
        bool exact = false;
        const double ps = spa_two_sided(S, gt, mu, &y, &exact);
        if (std::isfinite(ps) && ps > 0.0) {
            r.p = std::max(ps, 1e-300);
            r.spa = !exact;
            r.exact_tail = exact;
        }
    }
    return r;
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
    // One CHECKED factorisation serves both the solve and the log-determinant. A raw determinant
    // underflows to zero for a moderately sized q and then log() of it is meaningless; taking it from
    // the factor's diagonal is exact and also tells us the design was singular in the first place.
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
    // Checked LDLT rather than a raw inverse: a singular or near-singular design (a collinear
    // covariate, a monomorphic genotype after rotation) silently produced garbage through .inverse(),
    // where the factorisation reports it. Only the (q,q) entry of A^-1 is needed for the standard
    // error, so it comes from solving against the last unit vector.
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
    const double sigma2 = rss / static_cast<double>(n - (q + 1));
    const double var_b = sigma2 * Ainv_qq;
    if (!(var_b > 0.0) || !std::isfinite(var_b)) return r;
    r.beta = beta(q);
    r.se = std::sqrt(var_b);
    r.z = r.beta / r.se;
    // Same reasoning as the ordinary linear model: sigma^2 comes from the same residuals, so the tail
    // is Student-t on n - (q+1) degrees of freedom.
    r.p = student_p(r.z, static_cast<double>(n - (q + 1)));
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
        << "                         logistic reports a Rao SCORE test: `z`/`p` come from it, while\n"
        << "                         `log_or`/`se` stay the Wald maximum-likelihood effect size, so p is\n"
        << "                         NOT recoverable from log_or/se. The Wald test collapses for a rare\n"
        << "                         variant in an unbalanced study (near-separation inflates se); the\n"
        << "                         score test never fits the alternative, so it does not.\n"
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

    // ---- sample (genotype column) order ----
    std::vector<std::string> geno_samples;
    {
        GzLineReader r(opt.samples);
        if (!r.ok()) throw std::runtime_error("cannot open --samples: " + opt.samples);
        std::string line;
        while (r.getline(line)) { line = trim(line); if (!line.empty()) geno_samples.push_back(line); }
    }
    if (geno_samples.empty()) throw std::runtime_error("no samples in " + opt.samples);
    {
        // A repeated genotype column would silently reuse one individual's phenotype for two columns,
        // inflating n and correlating rows that are not independent observations.
        std::unordered_set<std::string> seen;
        for (const std::string& sname : geno_samples)
            if (!seen.insert(sname).second)
                throw std::runtime_error("--samples has a duplicate sample id: " + sname);
    }

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
            // Without this a repeated row silently overwrites the earlier one, so which phenotype was
            // analysed depends on file order -- a difference no output column would reveal.
            if (pheno.find(s) != pheno.end())
                throw std::runtime_error("--phenotype has a duplicate sample id: " + s);
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

    // Kinship for --model lmm / --pca: external GRM over used samples (rows/cols in --samples order).
    // panvar is local and won't build one from region genotypes (proximal contamination) -- supply a
    // genome-wide GRM. --pca appends PCs as fixed covariates; the LMM uses K as a random effect.
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
            throw std::runtime_error("--kinship matrix has fewer rows than samples (" +
                                     std::to_string(M.size()) + " < " +
                                     std::to_string(geno_samples.size()) + ")");
        // Row WIDTH was never checked, so a ragged matrix indexed past the end of a short row -- silent
        // out-of-bounds rather than an error. Check every row this run will actually touch.
        for (std::size_t c : used_cols)
            if (M[c].size() < geno_samples.size())
                throw std::runtime_error("--kinship row " + std::to_string(c + 1) + " has " +
                                         std::to_string(M[c].size()) + " entries, expected at least " +
                                         std::to_string(geno_samples.size()) + " (ragged matrix)");
        for (std::size_t a = 0; a < n_used; ++a)
            for (std::size_t b = 0; b < n_used; ++b)
                K(a, b) = M[used_cols[a]][used_cols[b]];

        // A GRM has to be finite and symmetric before anything downstream can mean what it claims: the
        // eigendecomposition below assumes self-adjointness and will happily return nonsense otherwise.
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
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K, Eigen::EigenvaluesOnly);
            if (es.info() != Eigen::Success)
                throw std::runtime_error("--kinship eigendecomposition failed; the matrix is not usable");
            const double lo = es.eigenvalues().minCoeff(), hi = es.eigenvalues().maxCoeff();
            // Scale the tolerance by max(1, |hi|) rather than hi: gating on `hi > 0` let a negative
            // DEFINITE matrix (-I has hi < 0) skip the check entirely.
            const double tol = 1e-6 * std::max(1.0, std::fabs(hi));
            for (std::size_t a = 0; a < n_used; ++a)
                if (K(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(a)) < -tol)
                    throw std::runtime_error("--kinship has a negative diagonal entry at " +
                                             std::to_string(a + 1) + "; it is not a covariance matrix");
            if (lo < -tol)
                throw std::runtime_error("--kinship is not positive semi-definite (smallest eigenvalue " +
                                         std::to_string(lo) + " against a largest of " + std::to_string(hi) +
                                         "); it is not a valid GRM");
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
        // Not quiet-suppressed: the existing GEMMA comparison checks the CORRELATION of beta and
        // -log10(p), which can look excellent while every p-value is off by a systematic factor. No
        // absolute per-feature tolerance on beta, standard error or p has ever been asserted against
        // a pinned reference, so the numbers below are not independently validated and must not be
        // described as if they were.
        std::cerr << "[associate] WARNING: the LMM is EXPERIMENTAL and numerically unvalidated. Its "
                     "only external check is a correlation against GEMMA, which cannot detect a "
                     "systematic difference in beta, SE or p. Use the linear/logistic models for "
                     "results that need to be defensible.\n";
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
        // Whether the effect estimate is usable. `ok` = the maximum-likelihood fit converged.
        // `separation` = it did not (logistic, near-complete separation), so log_or/se are absent while
        // the score test's p remains valid -- the score test never fits the alternative. Reporting the
        // p without saying the effect size is missing for a REASON is what makes that confusing.
        std::string effect_status = "ok";
        // Which tail produced `p`: "t" (Student-t, quantitative), "score" (Rao score, normal tail),
        // "score_spa" (score with the saddlepoint), "score_exact" (the score is at the edge of its
        // support, where the probability is a single Bernoulli configuration and is written down
        // exactly), "lmm". Worth naming per feature because one table mixes several.
        std::string p_method = "t";
        // Binary traits: the minor-allele carrier count split by case/control. Total MAC hides what
        // actually governs asymptotic reliability -- 1 case / 19 controls is far weaker than 10 / 10 at
        // the same MAC, and it is the imbalanced one that drives a score statistic into its bad tail.
        long mac_case = -1, mac_ctrl = -1;
    };
    std::vector<Row> rows;
    // Variant tier only: full-length mean-imputed dosage per retained row, for LD r^2 between variants.
    // Kept only in variant mode (few variants), so the k-mer substrate never holds a giant matrix.
    std::vector<Eigen::VectorXd> var_dose;
    // Which entries of var_dose were OBSERVED rather than mean-imputed. LD r^2 is happy with the imputed
    // vector, but a conditional model must not be: the marginal test uses complete cases, so a
    // conditional test on imputed rows analyses a different sample and the two p-values stop being
    // comparable -- which is precisely the comparison a COJO table invites a reader to make.
    std::vector<std::vector<char>> var_obs;
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
        std::string row_effect_status = "ok";
        std::string row_p_method = is_lmm ? "lmm" : "t";
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
            if (model == "logistic") {
                bool separated = false;
                // Keep the Wald fit's beta/se as the effect size, but report the SCORE test's p. The two
                // agree for common features and diverge exactly where the Wald one breaks -- a rare
                // variant in an unbalanced case/control study, where near-separation inflates the
                // standard error and collapses the statistic.
                const FitResult sc = score_logistic(X, yy, n, p_dim);
                // The contract is that a logistic p is ALWAYS the score test's. If the score test fails
                // there is no valid p to report: keeping the Wald one and labelling it `score` is the
                // silent fallback the contract exists to forbid, so the feature is dropped instead.
                if (!sc.ok) { fr.ok = false; }
                if (sc.ok) {
                    if (!fr.ok) {
                        // ML diverged. Recover a finite effect size with Firth's penalised likelihood
                        // rather than reporting none; the p-value stays the score test's either way.
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
                    fr.z = sc.z; fr.p = sc.p; fr.ok = true;
                }
                if (separated) row_effect_status = "separation";
                row_p_method = sc.exact_tail ? "score_exact" : (sc.spa ? "score_spa" : "score");
            }
        }
        if (!fr.ok) { ++n_dropped_fit; continue; }

        Row row;
        row.effect_status = row_effect_status;
        row.p_method = row_p_method;
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
        if (binary && g.size() == yy.size() && !g.empty()) {
            long long modal = 0; std::size_t best = 0;
            for (const auto& kv : cat) if (kv.second > best) { best = kv.second; modal = kv.first; }
            long ca = 0, co = 0;
            for (std::size_t t = 0; t < g.size(); ++t)
                if (std::llround(g[t]) != modal) { if (yy[t] > 0.5) ++ca; else ++co; }
            row.mac_case = ca; row.mac_ctrl = co;
        }
        row.n = n; row.minor_freq = minor_freq;
        row.beta = fr.beta; row.se = fr.se; row.z = fr.z; row.p = fr.p;
        rows.push_back(std::move(row));

        // Retain a full-length mean-imputed dosage for variant-tier LD-clumping (cheap: few variants).
        if (variant_mode) {
            Eigen::VectorXd vd(n_used);
            std::vector<char> obs(n_used, 1);
            if (is_lmm) {
                vd = glmm;  // already full-length, mean-imputed (the LMM rotation needs a fixed set)
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
                    if (!std::isfinite(vd(static_cast<Eigen::Index>(u)))) {
                        vd(static_cast<Eigen::Index>(u)) = mean;
                        obs[u] = 0;
                    }
            }
            var_dose.push_back(std::move(vd));
            var_obs.push_back(std::move(obs));
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
    std::size_t meff_clump = 0, meff_eigen = 0;   // clumping heuristic vs the phenotype-blind estimate
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
        // Every TESTED feature has to count toward the denominator it is corrected against. A low-AF
        // variant is barred above from ANCHORING a clump -- its r^2 is unstable, so it must not claim
        // shadows -- but it is still tested, and leaving it in no clump at all made Meff smaller than the
        // number of tests: its own Bonferroni threshold was then computed from a set it was not in.
        // Measured on the LPA cohort at --min-ac 30, 5 of 20 tested variants sat in no clump and Meff read
        // 12. Give each remaining tested feature a singleton clump: it counts, but it claims nothing.
        for (std::size_t oi = 0; oi < n_tests; ++oi) {
            const std::size_t i = byp[oi];
            if (rows[i].clump != -1 || !std::isfinite(rows[i].p)) continue;
            rows[i].clump = next_clump++;
            rows[i].is_lead = 1;
        }
        meff_clump = static_cast<std::size_t>(std::max(1, next_clump));
        // Phenotype-blind effective tests (Li & Ji 2005): eigenvalues of the genotype CORRELATION
        // matrix, Meff = sum_i [ I(lambda_i >= 1) + frac(lambda_i) ]. Clumping is seeded in p-value
        // order, so the phenotype changes how many clumps there are -- a chain A-B-C with A,C
        // uncorrelated gives one clump seeded at B and two seeded at A. This estimator never looks at
        // the phenotype, so the regional threshold it implies cannot be circular.
        if (n_tests >= 2) {
            Eigen::MatrixXd C(n_tests, n_tests);
            std::vector<double> mean(n_tests, 0.0), sd(n_tests, 0.0);
            for (std::size_t i = 0; i < n_tests; ++i) {
                const Eigen::VectorXd& v = var_dose[i];
                mean[i] = v.mean();
                sd[i] = std::sqrt((v.array() - mean[i]).square().sum() / std::max<double>(1.0, v.size() - 1));
            }
            bool usable = true;
            for (std::size_t i = 0; i < n_tests && usable; ++i) if (!(sd[i] > 0.0)) usable = false;
            if (usable) {
                for (std::size_t i = 0; i < n_tests; ++i)
                    for (std::size_t j = i; j < n_tests; ++j) {
                        const double cij = ((var_dose[i].array() - mean[i]) *
                                            (var_dose[j].array() - mean[j])).sum() /
                                           (std::max<double>(1.0, var_dose[i].size() - 1) * sd[i] * sd[j]);
                        C(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = cij;
                        C(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) = cij;
                    }
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C, Eigen::EigenvaluesOnly);
                if (es.info() == Eigen::Success) {
                    // Eigenvalues that are numerically zero must be treated as exactly zero: for k
                    // identical features the spectrum is (k, 0, ..., 0) and Li-Ji gives exactly 1, but
                    // the zeros arrive as +/-1e-16 whose fractional parts push the sum a hair above 1
                    // and a bare ceil() then answers 2. The same epsilon is why the sum is rounded up
                    // from just below rather than from the raw value.
                    // Li-Ji sums an INTEGER part and a FRACTIONAL part, so it is acutely sensitive to
                    // an eigenvalue that lands a few ulps below a whole number: for k identical
                    // features the top eigenvalue comes back as 2.9999999999999996, whose floor is 2 and
                    // whose fractional part is ~1, and the estimate reads 2 where the answer is 1. (R's
                    // own eigen() reproduces this, so it is a property of the estimator's form rather
                    // than of one implementation.) Snap near-integer eigenvalues before splitting them.
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
    } else if (n_tests > 0) {
        // Same rule as the variant tier: a tested feature must count. A feature with no bubble
        // annotation belongs to no block, and skipping it made Meff smaller than the number of tests
        // whenever only SOME features were annotated. Each unannotated feature counts as its own block.
        // With no annotation at all this reduces to meff = n_tests, which is where it started.
        std::unordered_set<std::string> blocks;
        std::size_t unannotated = 0;
        for (const Row& r : rows) {
            bool any = false;
            for (const std::string& b : split(r.bubbles, ';')) {
                const std::string t = trim(b);
                if (!t.empty() && t != ".") { blocks.insert(t); any = true; }
            }
            if (!any) ++unannotated;
        }
        if (!blocks.empty() || unannotated > 0)
            meff = std::max<std::size_t>(1, blocks.size() + unannotated);
    }

    // Variant-tier forward-stepwise conditional/joint (COJO), since r^2-clumping cannot establish
    // independence: a weak tag of a strong locus stays significant below --ld-r2. Add the smallest
    // conditional p until none clears --cojo-p, then report each conditioned on the rest, so shadows
    // collapse and true signals survive. Conditioning dosages enter as covariates with the genotype at
    // column 1, so FitResult is its conditional Wald. Linear/logistic only; LMM needs the rotation.
    std::vector<double> p_cond(rows.size(), kNaN);
    // Samples the CONDITIONAL model actually used. It is the intersection of the target's complete cases
    // with every conditioning feature's, so it can be smaller than the marginal `n` -- and then p and
    // p_conditional are computed on different samples and are not directly comparable. Reporting it is
    // the minimum; a reader can see when the comparison is clean and when it is not.
    std::vector<std::size_t> n_cond(rows.size(), 0);
    // cond_role: variant tier -> "signal" (COJO-selected) / "shadow"; feature tier -> "lead" /
    // "collinear" (same-event redundancy, not scored) / "conditioned"; "." when not computed.
    std::vector<std::string> cond_role(rows.size(), ".");
    std::size_t cojo_n_signals = 0;
    if (variant_mode && !is_lmm && rows.size() > 1) {
        // conditional p of variant i given a set of conditioning dosage vectors (empty -> marginal-style).
        auto cond_p = [&](std::size_t i, const std::vector<std::size_t>& cond, std::size_t* out_n) -> double {
            const std::size_t pc = p_dim + cond.size();        // intercept, g_i, covariates, |cond|
            // Complete cases over the target and every conditioning variant, matching what the marginal
            // test did. Under the LMM the rotation is over a fixed sample set, so everything is imputed
            // there by construction and every row is marked observed.
            std::vector<std::size_t> keep;
            keep.reserve(n_used);
            for (std::size_t u = 0; u < n_used; ++u) {
                if (!var_obs[i][u]) continue;
                bool ok = true;
                for (std::size_t c : cond) if (!var_obs[c][u]) { ok = false; break; }
                if (ok) keep.push_back(u);
            }
            const std::size_t nk = keep.size();
            if (out_n != nullptr) *out_n = nk;
            if (nk <= pc + 1) return kNaN;
            std::vector<double> X(nk * pc), yk(nk);
            for (std::size_t k = 0; k < nk; ++k) {
                const std::size_t u = keep[k];
                std::size_t col = 0;
                yk[k] = y[u];
                X[k * pc + col++] = 1.0;
                X[k * pc + col++] = var_dose[i](static_cast<Eigen::Index>(u));   // target -> column 1
                for (std::size_t j = 0; j < ncov_eff; ++j) X[k * pc + col++] = Z[u][j];
                for (std::size_t c : cond) X[k * pc + col++] = var_dose[c](static_cast<Eigen::Index>(u));
            }
            const FitResult fc = (model == "logistic") ? score_logistic(X, yk, nk, pc)
                                                       : fit_linear(X, yk, nk, pc);
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
                const double pc_i = cond_p(i, selected, nullptr);
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
            if (!cond.empty()) p_cond[i] = cond_p(i, cond, &n_cond[i]);   // empty cond (sole signal) -> NA
        }
    } else if (!variant_mode && !is_lmm && rows.size() > 1) {
        // --- feature-tier single-lead conditional p with a within-bubble collinearity guard ---
        // k-mer / node / edge dosages inside one bubble are ~collinear (they measure the same event), so
        // conditioning them on the top feature is degenerate -> flag features with r^2 > 0.95 vs the lead
        // ("collinear") instead of scoring them. Cross-bubble features get a real conditional p: a mere tag
        // of the lead's variant collapses. Two bounded streaming passes (no feature x sample matrix kept).
        // `obs` records which rows carried a real dosage, so the conditional fit below can use the same
        // complete cases the marginal test used rather than mean-imputed stand-ins.
        auto parse_full = [&](const std::vector<std::string>& f, Eigen::VectorXd& vd,
                              std::vector<char>* obs = nullptr) -> bool {
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
            if (obs != nullptr) obs->assign(n_used, 1);
            for (std::size_t u = 0; u < n_used; ++u)
                if (!std::isfinite(vd(static_cast<Eigen::Index>(u)))) {
                    vd(static_cast<Eigen::Index>(u)) = mean;
                    if (obs != nullptr) (*obs)[u] = 0;
                }
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
        std::vector<char> lead_obs;
        { GzLineReader g2(opt.genotypes); std::string line;
          while (g2.getline(line)) {
              if (line.empty()) continue;
              std::vector<std::string> f = split(line, ',');
              if (f.size() < need_cols) continue;
              if (trim(f[0]) == lead_id) { parse_full(f, gl, &lead_obs); break; }
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
                std::vector<char> obs;
                if (!parse_full(f, vd, &obs)) continue;
                if (r2_vec(vd, gl) > 0.95) { cond_role[i] = "collinear"; continue; }  // same-event redundancy
                std::vector<std::size_t> keep;
                keep.reserve(n_used);
                for (std::size_t u = 0; u < n_used; ++u)
                    if (obs[u] && (lead_obs.empty() || lead_obs[u])) keep.push_back(u);
                const std::size_t nk = keep.size();
                if (nk <= pc + 1) continue;
                std::vector<double> X(nk * pc), yk(nk);
                for (std::size_t k = 0; k < nk; ++k) {
                    const std::size_t u = keep[k];
                    yk[k] = y[u];
                    X[k * pc + 0] = 1.0;
                    X[k * pc + 1] = vd(static_cast<Eigen::Index>(u));
                    for (std::size_t j = 0; j < ncov_eff; ++j) X[k * pc + 2 + j] = Z[u][j];
                    X[k * pc + 2 + ncov_eff] = gl(static_cast<Eigen::Index>(u));
                }
                const FitResult fc = (model == "logistic") ? score_logistic(X, yk, nk, pc)
                                                           : fit_linear(X, yk, nk, pc);
                if (fc.ok) { p_cond[i] = fc.p; n_cond[i] = nk; cond_role[i] = "conditioned"; }
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
    // lambda_GC is a genome-wide diagnostic: it reads the MEDIAN chi-square, on the assumption that most
    // tests are null. panvar tests one locus, where a real signal and everything in LD with it can be
    // most of the tests -- lambda then measures the signal, not inflation. Reporting it unqualified
    // invites the opposite reading, so say which situation this run is in.
    const bool lambda_meaningful =
        n_tests >= 100 && static_cast<double>(n_sig_fdr) < 0.25 * static_cast<double>(n_tests);
    const std::string lambda_label = lambda_meaningful
        ? std::string("genomic inflation lambda_GC")
        : std::string("lambda_GC (NOT an inflation estimate here: too few tests, or too many are "
                      "signal -- it reads the median chi-square assuming most tests are null)");

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
        // p_bonf = raw 0.05/n_tests. p_bonf_meff scales by Meff, an LD-CLUMPING heuristic that is
        // seeded in p-value order and therefore depends on the phenotype -- it is not a phenotype-blind
        // effective-test count and gives no formal family-wise guarantee. q_bh is the primary control.
        out << "feature_id\tlayer\tbubbles\tnodes\tn\tn_conditional\tminor_freq\t" << effect
            << "\tse\tz\tp\tp_method\teffect_status\tmac_case\tmac_ctrl\tp_bonf\tp_bonf_meff\tq_bh\taf\tan\tlow_af\tclump\tis_lead\tgene"
               "\tp_conditional\tcond_role\n";
        for (std::size_t k = 0; k < n_tests; ++k) {
            const std::size_t i = order[k];
            const Row& r = rows[i];
            const double p_bonf = std::isfinite(r.p) ? std::min(1.0, r.p * static_cast<double>(n_tests)) : kNaN;
            const double p_bonf_meff = std::isfinite(r.p) ? std::min(1.0, r.p * static_cast<double>(meff)) : kNaN;
            out << r.id << '\t' << r.layer << '\t' << r.bubbles << '\t' << r.nodes << '\t'
                << r.n << '\t' << (n_cond[i] ? std::to_string(n_cond[i]) : std::string("."))
                << '\t' << fmt(r.minor_freq) << '\t' << fmt(r.beta) << '\t' << fmt(r.se)
                << '\t' << fmt(r.z) << '\t' << fmt(r.p) << '\t' << r.p_method << '\t' << r.effect_status
                << '\t' << (r.mac_case < 0 ? std::string(".") : std::to_string(r.mac_case))
                << '\t' << (r.mac_ctrl < 0 ? std::string(".") : std::to_string(r.mac_ctrl))
                << '\t' << fmt(p_bonf) << '\t' << fmt(p_bonf_meff)
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
            << "meff_method\t" << (meff_eigen > 0 ? "eigenvalue (Li-Ji, phenotype-blind)"
                                                   : (variant_mode ? "ld_clumping (heuristic)" : "bubbles"))
               << '\n'
            << "meff_eigen\t" << (meff_eigen > 0 ? std::to_string(meff_eigen) : std::string("NA")) << '\n'
            << "meff_ld_clumping\t" << (meff_clump > 0 ? std::to_string(meff_clump) : std::string("NA")) << '\n'
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

    cli::RunLog log("associate", opt.quiet);
    log.info("model " + model + " (" + (binary ? "binary" : "quantitative") + "), " +
             std::to_string(n_used) + " samples, " + std::to_string(ncov_eff) + " covariates" +
             (opt.pca > 0 ? " (incl. " + std::to_string(opt.pca) + " PCs)" : ""));
    log.info("unit " + std::string(variant_mode ? "variant" : "feature") + "; tested " +
             std::to_string(n_tests) + ", Meff " + std::to_string(meff) + ", dropped " +
             std::to_string(n_dropped_maf) + " (MAF) + " + std::to_string(n_dropped_fit) + " (fit)");
    if (binary) {
        // The score test is well calibrated in the body but still mildly anti-conservative in the far
        // tail for very rare features, and a regional Bonferroni threshold can sit right there. SPA
        // IS implemented and is applied past |z| > 2; it moved the measured type-I error at p<0.001
        // from 0.0025 to 0.0017 against a nominal 0.001, so roughly 1.7x rather than 2.5x. It does
        // not remove the anti-conservatism and it is not rare-variant aggregation, so the warning
        // stands -- what changed is that the remedy is no longer missing, only insufficient.
        std::size_t n_rare = 0, n_imbal = 0;
        for (const Row& r : rows) {
            if (std::isfinite(r.minor_freq) && r.minor_freq < 0.01) ++n_rare;
            if (r.mac_case >= 0 && r.mac_case < 10) ++n_imbal;
        }
        if (n_rare > 0 || n_imbal > 0)
            log.info("WARNING: " + std::to_string(n_rare) + " features below 1% minor frequency and " +
                     std::to_string(n_imbal) + " with fewer than 10 minor-allele carriers among cases. "
                     "The score test is saddlepoint-corrected, but measured type-I error in the far "
                     "tail is still about 1.7x nominal at p<0.001 for such features, and panvar has no "
                     "burden/SKAT rare-variant test -- treat these as exploratory, not as calibrated "
                     "discoveries");
    }
    log.info("significant: Bonferroni(Meff=" + std::to_string(meff) + ") " + std::to_string(n_sig_meff) +
             ", FDR<0.05 " + std::to_string(n_sig_fdr) + "; " + lambda_label + "=" + fmt(lambda_gc) +
             (variant_mode ? "; COJO signals " + std::to_string(cojo_n_signals) : ""));
    log.wrote({opt.out_prefix + ".assoc.tsv", opt.out_prefix + ".summary.tsv"});
    log.done();
    return 0;
}

}  // namespace panvar
