#pragma once

#include <array>

namespace eri::math {

inline std::array<double,3> gaussian_product_center(
    double a, const std::array<double,3>& A,
    double b, const std::array<double,3>& B) noexcept {
    const double g = a + b;
    return {
        (a*A[0] + b*B[0]) / g,
        (a*A[1] + b*B[1]) / g,
        (a*A[2] + b*B[2]) / g
    };
}

} // namespace eri::math
