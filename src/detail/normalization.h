#pragma once

#include <cmath>
#include "math/double_factorial.h"

namespace eri::detail {

double primitive_norm(int lx, int ly, int lz, double alpha) noexcept {
    const int L = lx + ly + lz;

    const double nom =
        std::pow(2.0, 2.0 * L + 1.5) *
        std::pow(alpha, L + 1.5);

    const double denom =
        eri::math::double_factorial(2 * lx - 1) *
        eri::math::double_factorial(2 * ly - 1) *
        eri::math::double_factorial(2 * lz - 1) *
        std::pow(M_PI, 1.5);

    return std::sqrt(nom / denom);
}

} // end namespace eri::detail