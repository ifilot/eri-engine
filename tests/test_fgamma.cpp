#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <vector>
#include <cmath>
#include <limits>

#include <eri/math/fgamma.h>

using eri::math::Fgamma_acc;
using eri::math::Fgamma_nr;
using eri::math::Fgamma_block_interp;

/*==============================================================================
  Test: Interpolated Boys vs accurate reference
  ==============================================================================*/
/**
 * @brief Validate the interpolation-based Boys evaluator.
 *
 * This test compares:
 *   - Fgamma_block_interp (tabulation + recurrence)
 * against:
 *   - Fgamma_acc          (high-accuracy reference)
 *
 * over a wide log-spaced range of T.
 *
 * Expected accuracy:
 *   - absolute error  ≲ 1e-12
 *   - relative error  ≲ 1e-10
 *
 * These tolerances are appropriate for:
 *   - linear interpolation in log(T)
 *   - NTABLE = 2048
 *   - double precision arithmetic
 */
TEST_CASE("Boys Fgamma interpolation matches accurate reference")
{
    constexpr int nu_max = 24;
    constexpr int NT     = 2000;

    // Initialize interpolation table once
    eri::math::init_Fgamma_interp_table(nu_max);

    std::vector<double> F_interp(nu_max + 1);
    std::vector<double> F_exact (nu_max + 1);

    double max_abs_err = 0.0;
    double max_rel_err = 0.0;

    int    worst_nu_abs = -1;
    int    worst_nu_rel = -1;
    double worst_T_abs  = 0.0;
    double worst_T_rel  = 0.0;

    // Log-spaced sweep: very small -> moderate T
    for (int i = 0; i < NT; ++i) {
        const double logT = -12.0 + 6.0 * double(i) / double(NT - 1);
        const double T    = std::pow(10.0, logT);

        Fgamma_block_interp(nu_max, T, F_interp.data());

        for (int nu = 0; nu <= nu_max; ++nu) {
            F_exact[nu] = Fgamma_acc(nu, T);

            const double abs_err =
                std::abs(F_interp[nu] - F_exact[nu]);

            const double rel_err =
                abs_err / std::max(std::abs(F_exact[nu]), 1e-16);

            if (abs_err > max_abs_err) {
                max_abs_err  = abs_err;
                worst_nu_abs = nu;
                worst_T_abs  = T;
            }

            if (rel_err > max_rel_err) {
                max_rel_err  = rel_err;
                worst_nu_rel = nu;
                worst_T_rel  = T;
            }
        }
    }

    INFO("Max absolute error = " << max_abs_err
         << " at nu=" << worst_nu_abs
         << ", T=" << worst_T_abs);

    INFO("Max relative error = " << max_rel_err
         << " at nu=" << worst_nu_rel
         << ", T=" << worst_T_rel);

    CHECK(max_abs_err < 1e-12);
    CHECK(max_rel_err < 1e-10);
}


/*==============================================================================
  Test: Numerical Recipes Boys implementation (diagnostic only)
  ==============================================================================*/
/**
 * @brief Diagnostic evaluation of the Numerical-Recipes-based Boys function.
 *
 * This test exists solely to characterize and document the behavior of
 * Fgamma_nr.  It does NOT assert numerical agreement with the accurate
 * reference implementation.
 *
 * Observed behavior:
 *   - Fgamma_nr can exhibit very large absolute and relative errors
 *     (orders of magnitude) for moderate ν and T.
 *   - Accuracy is highly non-uniform and domain-dependent.
 *   - The implementation is mathematically correct but numerically
 *     ill-conditioned for Boys-function evaluation.
 *
 * Purpose:
 *   - ensure Fgamma_nr remains finite in restricted domains
 *   - detect catastrophic failures (NaN / Inf)
 *   - preserve historical/diagnostic visibility
 *
 * IMPORTANT:
 *   This test MUST NOT enforce numeric accuracy tolerances.
 *   Fgamma_nr is NOT a valid reference for Boys functions.
 */
TEST_CASE("Boys Fgamma_nr produces finite, sane values in restricted domain")
{
    constexpr int nu_max = 12;
    constexpr int NT     = 500;

    int bad_count = 0;

    for (int i = 0; i < NT; ++i) {
        // Restrict to region where NR is least pathological
        const double logT = -3.0 + 3.0 * double(i) / double(NT - 1); // 1e-3 .. 1
        const double T    = std::pow(10.0, logT);

        for (int nu = 0; nu <= nu_max; ++nu) {
            const double v = Fgamma_nr(nu, T);

            // It should at least return a finite number
            if (!std::isfinite(v)) {
                ++bad_count;
                INFO("Non-finite NR value at nu=" << nu << ", T=" << T);
            }

            // Boys functions are strictly positive
            if (std::isfinite(v)) {
                CHECK(v >= 0.0);
            }
        }
    }

    // Allow zero failures only
    CHECK(bad_count == 0);
}

/**
 * @brief Compute Gauss-Legendre nodes and weights on [-1, 1].
 *
 * @param n   Number of quadrature points
 * @param x   Output nodes (size n)
 * @param w   Output weights (size n)
 */
inline void gauss_legendre(int n,
                           std::vector<double>& x,
                           std::vector<double>& w)
{
    constexpr double EPS = 1e-15;

    x.resize(n);
    w.resize(n);

    const int m = (n + 1) / 2;
    const double pi = std::acos(-1.0);

    for (int i = 0; i < m; ++i) {
        // Initial guess (Abramowitz & Stegun)
        double z = std::cos(pi * (i + 0.75) / (n + 0.5));
        double z1;

        // Newton iteration
        do {
            double p1 = 1.0;
            double p2 = 0.0;

            for (int j = 1; j <= n; ++j) {
                const double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }

            // p1 = P_n(z), p2 = P_{n-1}(z)
            const double pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z  = z1 - p1 / pp;
        }
        while (std::abs(z - z1) > EPS);

        // Symmetric nodes
        x[i]         = -z;
        x[n - 1 - i] =  z;

        // Weight
        double p1 = 1.0;
        double p2 = 0.0;
        for (int j = 1; j <= n; ++j) {
            const double p3 = p2;
            p2 = p1;
            p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
        }
        const double pp = n * (z * p1 - p2) / (z * z - 1.0);
        const double wval = 2.0 / ((1.0 - z * z) * pp * pp);

        w[i]         = wval;
        w[n - 1 - i] = wval;
    }
}

/**
 * @brief Numerically compute Boys function using Gauss-Legendre quadrature.
 *
 * This is intended ONLY for validation / testing.
 * It is slow but very accurate for large ν and small T.
 *
 * @param nu  Boys order (nu >= 0)
 * @param T   Boys argument (T >= 0)
 * @param N   Number of quadrature points (64-256 recommended)
 */
inline double Fgamma_legendre(int nu, double T, int N = 128)
{
    if (nu < 0 || T < 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    std::vector<double> x, w;
    gauss_legendre(N, x, w);

    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
        const double t  = 0.5 * (x[i] + 1.0);  // map to [0,1]
        const double t2 = t * t;

        // integrand: t^(2nu) * exp(-T t^2)
        const double f =
            std::pow(t, 2.0 * nu) * std::exp(-T * t2);

        sum += w[i] * f;
    }

    return 0.5 * sum;
}

TEST_CASE("Fgamma_acc matches Gauss-Legendre quadrature (large nu, small T)")
{
    struct Case { int nu; double T; };
    const std::vector<Case> cases = {
        {10, 1e-12}, {20, 1e-12}, {30, 1e-12},
        {10, 1e-8 }, {20, 1e-8 }, {30, 1e-8 },
        {10, 1e-6 }, {20, 1e-6 }, {30, 1e-6 }
    };

    constexpr int N = 128; // increase to 256 if you want margin

    double worst_abs = 0.0;
    double worst_rel = 0.0;

    for (const auto& c : cases) {
        const double ref  = eri::math::Fgamma_acc(c.nu, c.T);
        const double quad = Fgamma_legendre(c.nu, c.T, N);

        const double abs_err = std::abs(ref - quad);
        const double rel_err =
            abs_err / std::max(std::abs(ref), 1e-16);

        worst_abs = std::max(worst_abs, abs_err);
        worst_rel = std::max(worst_rel, rel_err);
    }

    INFO("Worst abs error = " << worst_abs);
    INFO("Worst rel error = " << worst_rel);
    INFO("Legendre N = " << N);

    CHECK(worst_abs < 1e-15);
    CHECK(worst_rel < 1e-14);
}

TEST_CASE("Fgamma_nr violates monotonicity for large nu and small T")
{
    constexpr double T = 1e-8;   // pathological regime
    constexpr int nu_max = 40;

    bool violation_found = false;

    double prev = eri::math::Fgamma_nr(0, T);

    for (int nu = 1; nu <= nu_max; ++nu) {
        const double curr = eri::math::Fgamma_nr(nu, T);

        INFO("nu=" << nu << ", Fgamma_nr=" << curr);

        if (!std::isfinite(curr)) {
            violation_found = true;
            break;
        }

        // Boys functions must be strictly decreasing in nu
        if (curr >= prev) {
            violation_found = true;
            break;
        }

        prev = curr;
    }

    CHECK(violation_found);
}

TEST_CASE("Fgamma_acc preserves monotonicity for large nu and small T")
{
    // Same pathological regime used to expose NR failures
    constexpr double T = 1e-8;
    constexpr int nu_max = 40;

    double prev = eri::math::Fgamma_acc(0, T);

    REQUIRE(std::isfinite(prev));
    REQUIRE(prev > 0.0);

    for (int nu = 1; nu <= nu_max; ++nu) {
        const double curr = eri::math::Fgamma_acc(nu, T);

        INFO("nu=" << nu << ", Fgamma_acc=" << curr);

        // Sanity checks
        CHECK(std::isfinite(curr));
        CHECK(curr > 0.0);

        // Core invariant: strict monotonic decrease
        CHECK(curr < prev);

        prev = curr;
    }
}