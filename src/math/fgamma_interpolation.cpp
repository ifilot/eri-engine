#include "math/fgamma_interpolation.h"

namespace eri::math {

constexpr double T_SMALL = 1e-6;
constexpr double T_LARGE = 30.0;
constexpr int    NTABLE  = 2048;

// ---------- Internal table state (private) ----------
namespace {

int nu_tab_max = -1;
double logT_min = 0.0, logT_max = 0.0, inv_dlogT = 0.0;
std::vector<double> Htab;

inline int idx(int nu, int i)
{
    return nu * NTABLE + i;
}

} // anonymous namespace

// ---------- Initialization ----------
void init_Fgamma_interp_table(int nu_table_max_in)
{
    nu_tab_max = std::max(0, nu_table_max_in);
    Htab.assign((nu_tab_max + 1) * NTABLE, 0.0);

    logT_min  = std::log(T_SMALL);
    logT_max  = std::log(T_LARGE);
    inv_dlogT = (NTABLE - 1) / (logT_max - logT_min);

    for (int i = 0; i < NTABLE; ++i) {
        const double logT = logT_min + (logT_max - logT_min) * (double)i / (double)(NTABLE - 1);
        const double T     = std::exp(logT);
        const double sqrtT = std::sqrt(T);

        double Tpow = sqrtT; // T^(nu+1/2), nu=0
        for (int nu = 0; nu <= nu_tab_max; ++nu) {
            const double Fnu = eri::math::Fgamma_acc(nu, T);
            Htab[idx(nu, i)] = Tpow * Fnu;
            Tpow *= T;
        }
    }
}

// ---------- Internal interpolation ----------
static inline double Hgamma_interp(int nu, double T)
{
    const double x = (std::log(T) - logT_min) * inv_dlogT;
    int i = static_cast<int>(x);
    if (i < 0) {
        i = 0;
    }
    if (i > NTABLE - 2) {
        i = NTABLE - 2;
    }

    const double w = x - i;
    const double h0 = Htab[idx(nu, i)];
    const double h1 = Htab[idx(nu, i + 1)];
    return (1.0 - w) * h0 + w * h1;
}

// ---------- Public interpolated Fgamma ----------
double Fgamma_interp(int nu, double T)
{
    if (T < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (T < T_SMALL) {
        return eri::math::Fgamma_acc(nu, T);
    }

    if (nu_tab_max < 0 || nu > nu_tab_max) {
        return eri::math::Fgamma_acc(nu, T);
    }

    if (T <= T_LARGE) {
        const double H = Hgamma_interp(nu, T);

        const double sqrtT = std::sqrt(T);
        double Tnu = 1.0;
        for (int k = 0; k < nu; ++k) Tnu *= T;

        return H / (sqrtT * Tnu);
    }

    return eri::math::Fgamma_acc(nu, T);
}

// ---------- Block evaluator ----------
void Fgamma_block_interp(int nu_max, double T, double* F)
{
    if (T < T_SMALL) {
        for (int nu = 0; nu <= nu_max; ++nu) {
            F[nu] = eri::math::Fgamma_acc(nu, T);
        }
        return;
    }

    const double expmT = std::exp(-T);

    if (T > 2.0 * nu_max) {
        // -------- UPWARD --------
        F[0] = Fgamma_interp(0, T);
        const double inv2T = 1.0 / (2.0 * T);

        for (int nu = 0; nu < nu_max; ++nu) {
            F[nu + 1] = ((2.0 * nu + 1.0) * F[nu] - expmT) * inv2T;
        }
    }
    else {
        // -------- DOWNWARD --------
        F[nu_max] = Fgamma_interp(nu_max, T);

        for (int nu = nu_max; nu > 0; --nu) {
            F[nu - 1] = (2.0 * T * F[nu] + expmT) / (2.0 * nu - 1.0);
        }
    }
}

} // namespace eri::math
