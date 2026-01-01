#include "math/fgamma_interpolation.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace eri::math {

/*==============================================================================
  Tunable parameters controlling interpolation behavior
==============================================================================*/

/**
 * @brief Lower bound for interpolation domain in T.
 *
 * For T < T_SMALL we bypass interpolation entirely and fall back to the
 * accurate reference implementation.  This avoids numerical issues near
 * the removable singularity at T → 0.
 */
constexpr double T_SMALL = 1e-6;

/**
 * @brief Upper bound for interpolation domain in T.
 *
 * For T > T_LARGE the Boys function is smooth and recurrence is stable;
 * direct evaluation via the accurate reference is sufficiently fast.
 */
constexpr double T_LARGE = 30.0;

/**
 * @brief Number of tabulation points in log(T).
 *
 * The table is uniformly spaced in log(T) on [T_SMALL, T_LARGE].
 * Values around 1024–2048 typically provide ~1e-12 accuracy with
 * linear interpolation.
 */
constexpr int NTABLE = 2048;


/*==============================================================================
  Internal table state (translation-unit local)
==============================================================================*/

namespace {

int nu_tab_max = -1;     ///< Maximum order ν stored in the table
double logT_min = 0.0;  ///< log(T_SMALL)
double logT_max = 0.0;  ///< log(T_LARGE)
double inv_dlogT = 0.0; ///< 1 / ΔlogT

/**
 * @brief Interpolation table for scaled Boys functions.
 *
 * We tabulate:
 *
 *   H_ν(T) = T^{ν+1/2} · F_ν(T)
 *
 * rather than F_ν(T) itself.  This removes the dominant power-law
 * behavior and yields a much smoother function in log(T), making
 * linear interpolation accurate.
 *
 * Storage layout:
 *   Htab[ ν * NTABLE + i ],  i = log-grid index
 */
std::vector<double> Htab;

/** Compute flat array index for table lookup */
inline int idx(int nu, int i) noexcept {
    return nu * NTABLE + i;
}

} // anonymous namespace


/*==============================================================================
  Table initialization
==============================================================================*/

/**
 * @brief Initialize the Boys-function interpolation table.
 *
 * This function must be called once before any use of Fgamma_interp()
 * or Fgamma_block_interp().
 *
 * The table is constructed using the *accurate* Boys implementation
 * (Fgamma_acc), not the Numerical-Recipes variant.  Interpolation is
 * therefore anchored to a high-accuracy reference.
 *
 * @param nu_table_max  Maximum Boys order ν to precompute.
 */
void init_Fgamma_interp_table(int nu_table_max)
{
    nu_tab_max = std::max(0, nu_table_max);
    Htab.assign((nu_tab_max + 1) * NTABLE, 0.0);

    logT_min  = std::log(T_SMALL);
    logT_max  = std::log(T_LARGE);
    inv_dlogT = (NTABLE - 1) / (logT_max - logT_min);

    for (int i = 0; i < NTABLE; ++i) {
        const double logT =
            logT_min + (logT_max - logT_min) * double(i) / double(NTABLE - 1);

        const double T     = std::exp(logT);
        const double sqrtT = std::sqrt(T);

        // Build T^(ν+1/2) incrementally: sqrt(T) * T^ν
        double Tpow = sqrtT;

        for (int nu = 0; nu <= nu_tab_max; ++nu) {
            const double Fnu = eri::math::Fgamma_acc(nu, T);
            Htab[idx(nu, i)] = Tpow * Fnu;
            Tpow *= T;
        }
    }
}


/*==============================================================================
  Internal interpolation helper
==============================================================================*/

/**
 * @brief Interpolate H_ν(T) = T^{ν+1/2} F_ν(T) in log(T).
 *
 * Assumes:
 *  - T ∈ [T_SMALL, T_LARGE]
 *  - 0 ≤ ν ≤ nu_tab_max
 *  - init_Fgamma_interp_table() has been called
 */
static inline double Hgamma_interp(int nu, double T) noexcept
{
    const double x = (std::log(T) - logT_min) * inv_dlogT;

    int i = static_cast<int>(x);
    if (i < 0)           i = 0;
    if (i > NTABLE - 2)  i = NTABLE - 2;

    const double w  = x - double(i);
    const double h0 = Htab[idx(nu, i)];
    const double h1 = Htab[idx(nu, i + 1)];

    return (1.0 - w) * h0 + w * h1;
}


/*==============================================================================
  Public interpolated Boys function
==============================================================================*/

/**
 * @brief Approximate Boys function F_ν(T) using interpolation.
 *
 * Strategy:
 *  - T < T_SMALL      → accurate reference implementation
 *  - T ∈ table range  → interpolate scaled function + rescale
 *  - T > T_LARGE      → accurate reference implementation
 *  - ν > table range  → accurate reference implementation
 *
 * @param nu  Boys order ν ≥ 0
 * @param T   Boys argument T ≥ 0
 */
double Fgamma_interp(int nu, double T)
{
    if (T < 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    // Small-T and out-of-range fallbacks
    if (T < T_SMALL || nu_tab_max < 0 || nu > nu_tab_max)
        return eri::math::Fgamma_acc(nu, T);

    // Interpolation regime
    if (T <= T_LARGE) {
        const double H = Hgamma_interp(nu, T);

        // Reconstruct T^(ν+1/2) without pow()
        const double sqrtT = std::sqrt(T);
        double Tnu = 1.0;
        for (int k = 0; k < nu; ++k) Tnu *= T;

        return H / (sqrtT * Tnu);
    }

    // Large-T fallback
    return eri::math::Fgamma_acc(nu, T);
}


/*==============================================================================
  Block evaluation with stable recurrence
==============================================================================*/

/**
 * @brief Compute Boys functions F_0(T) … F_νmax(T) efficiently.
 *
 * Uses:
 *  - interpolation to obtain a single anchor value
 *  - stable upward or downward recurrence depending on T
 *
 * This is the intended high-performance interface for ERI evaluation.
 *
 * @param nu_max  Maximum Boys order
 * @param T       Boys argument
 * @param F       Output array of size ≥ nu_max+1
 */
void Fgamma_block_interp(int nu_max, double T, double* F)
{
    // Very small T: use accurate implementation directly
    if (T < T_SMALL) {
        for (int nu = 0; nu <= nu_max; ++nu)
            F[nu] = eri::math::Fgamma_acc(nu, T);
        return;
    }

    const double expmT = std::exp(-T);

    // Choose stable recurrence direction
    if (T > 2.0 * nu_max) {
        // -------- Upward recurrence --------
        F[0] = Fgamma_interp(0, T);
        const double inv2T = 1.0 / (2.0 * T);

        for (int nu = 0; nu < nu_max; ++nu) {
            F[nu + 1] =
                ((2.0 * nu + 1.0) * F[nu] - expmT) * inv2T;
        }
    }
    else {
        // -------- Downward recurrence --------
        F[nu_max] = Fgamma_interp(nu_max, T);

        for (int nu = nu_max; nu > 0; --nu) {
            F[nu - 1] =
                (2.0 * T * F[nu] + expmT) / (2.0 * nu - 1.0);
        }
    }
}

} // namespace eri::math
