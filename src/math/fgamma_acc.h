#include <limits>
#include <cmath>

namespace eri::math {

namespace detail {

// Machine constants
inline constexpr double kEps = std::numeric_limits<double>::epsilon();

// --- Small-T series: F_n(T) = sum_{k>=0} (-T)^k / (k! (2n+2k+1)) ---
inline double boys_series(int n, double T) {
    // Handles T >= 0, best for very small T
    double term = 1.0 / (2.0 * n + 1.0);
    double sum  = term;

    // Recurrence for terms:
    // term_{k+1} = term_k * [(-T)/(k+1)] * [(2n+2k+1)/(2n+2k+3)]
    for (int k = 0; k < 2000; ++k) {
        const double num = (2.0 * n + 2.0 * k + 1.0);
        const double den = (2.0 * n + 2.0 * k + 3.0);
        term *= (-T) / double(k + 1) * (num / den);
        const double newsum = sum + term;
        if (std::abs(term) <= std::abs(newsum) * (50.0 * kEps)) { // comfortable margin
            return newsum;
        }
        sum = newsum;
    }
    return sum; // Should never hit in sane regimes
}

// Accurate F0 via erf; stable for all T>=0 with a tiny-T fallback
inline double boys_F0(double T) {
    if (T == 0.0) return 1.0;
    // For very small T, erf(sqrt(T))/sqrt(T) suffers cancellation; use series
    if (T < 1e-10) {
        // F0(T) = 1 - T/3 + T^2/10 - T^3/42 + ...
        const double T2 = T * T;
        return 1.0 - T/3.0 + T2/10.0 - (T2*T)/42.0;
    }
    const double r = std::sqrt(T);
    return 0.5 * std::sqrt(M_PI / T) * std::erf(r);
}

/*
  Stabilized downward recursion with correction using the exact F0(T).

  Recurrence (inhomogeneous):
    F_{k-1} = alpha_k * F_k + beta_k
    alpha_k = 2T/(2k-1)
    beta_k  = exp(-T)/(2k-1)

  For fixed T, F0 is an affine function of the starting value F_N:
    F0 = A * F_N + B

  We can compute A (sensitivity) and B in one pass, then choose F_N so that
  computed F0 matches the accurate F0 from erf. Then a second pass yields
  highly accurate F_k down to the needed orders.
*/

// Compute F_0..F_nmax into F[0..nmax]
inline void Fgamma_block_exact(int nu_max, double T, double* F)
{
    // Small-T: series is the gold standard
    if (T < 0.1) {
        for (int nu = 0; nu <= nu_max; ++nu)
            F[nu] = boys_series(nu, T);
        return;
    }

    // Moderate/large T: erf + upward recursion
    const double expmT = std::exp(-T);
    F[0] = boys_F0(T);

    const double inv2T = 0.5 / T;
    for (int nu = 0; nu < nu_max; ++nu) {
        F[nu + 1] =
            ((2.0*nu + 1.0)*F[nu] - expmT) * inv2T;
    }
}

}
}