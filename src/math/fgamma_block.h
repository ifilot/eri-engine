#pragma once
#include <cmath>
#include <vector>

#include "math/fgamma.h"

namespace eri::math {

// Assumed to already exist:
// inline double Fgamma(int nu, double T);

inline void Fgamma_block(int nu_max, double T, double* F) {
    constexpr double T_SMALL = 1e-8;
    
    // Exact T=0 limit
    if (T < T_SMALL) {
        for (int nu = 0; nu <= nu_max; ++nu) {
            F[nu] = 1.0 / (2.0*nu + 1.0);
        }
        return;
    }

    const double expmT = std::exp(-T);

    if (T > 2.0 * nu_max) {  // upward recursion
        F[0] = eri::math::Fgamma(0, T);

        const double inv2T = 1.0 / (2.0 * T);

        for (int nu = 0; nu < nu_max; ++nu) {
            F[nu + 1] = ((2.0 * nu + 1.0) * F[nu] - expmT) * inv2T;
        }
    }
    else {                  // downward recursion
        F[nu_max] = eri::math::Fgamma(nu_max, T);

        for (int nu = nu_max; nu > 0; --nu) {
            F[nu - 1] = (2.0 * T * F[nu] + expmT) / (2.0 * nu - 1.0);
        }
    }
}

} // namespace eri::math
