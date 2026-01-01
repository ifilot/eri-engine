#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <vector>
#include <cmath>
#include <limits>

#include <eri/math/fgamma.h>

using eri::math::Fgamma_acc;
using eri::math::Fgamma_block_interp;

TEST_CASE("Boys Fgamma interpolation matches exact values") {
    constexpr int nu_max = 24;
    constexpr int NT = 2000;

    // Initialize interpolation table once
    eri::math::init_Fgamma_interp_table(nu_max);

    std::vector<double> F_interp(nu_max + 1);
    std::vector<double> F_exact(nu_max + 1);

    double max_abs_err = 0.0;
    double max_rel_err = 0.0;

    int worst_nu_abs = -1;
    int worst_nu_rel = -1;
    double worst_T_abs = 0.0;
    double worst_T_rel = 0.0;

    // Log-spaced T sweep: very small -> moderate
    for (int i = 0; i < NT; ++i) {
        const double logT = -12.0 + 6.0 * double(i) / double(NT - 1);
        const double T = std::pow(10.0, logT);

        Fgamma_block_interp(nu_max, T, F_interp.data());
        for (int nu = 0; nu <= nu_max; ++nu) {
            F_exact[nu] = Fgamma_acc(nu, T);

            const double abs_err =
                std::abs(F_interp[nu] - F_exact[nu]);

            const double rel_err =
                abs_err / std::max(std::abs(F_exact[nu]), 1e-16);

            if (abs_err > max_abs_err) {
                max_abs_err = abs_err;
                worst_nu_abs = nu;
                worst_T_abs = T;
            }

            if (rel_err > max_rel_err) {
                max_rel_err = rel_err;
                worst_nu_rel = nu;
                worst_T_rel = T;
            }
        }
    }

    INFO("Max absolute error = " << max_abs_err
         << " at nu=" << worst_nu_abs
         << ", T=" << worst_T_abs);

    INFO("Max relative error = " << max_rel_err
         << " at nu=" << worst_nu_rel
         << ", T=" << worst_T_rel);

    // Tolerances chosen for linear interpolation + recurrence
    CHECK(max_abs_err < 1e-12);
    CHECK(max_rel_err < 1e-10);
}
