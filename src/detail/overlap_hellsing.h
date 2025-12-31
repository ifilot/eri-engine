#pragma once

#include <cmath>
#include "math/factorial.h"
#include "math/gaussian_product.h"
#include "math/binomial_prefactor.h"

namespace eri::detail {

inline double overlap_1d_hellsing(int l1, int l2, double a1, double a2, double x, double gamma) {
    double pre = eri::math::factorial(l1) * eri::math::factorial(l2) / std::pow(gamma, l1 + l2);
    if (l1 & 1) pre *= -1.0;

    double sm = 0.0;

    for (int i1 = 0; i1 <= l1/2; ++i1) {
        for (int i2 = 0; i2 <= l2/2; ++i2) {
            const int j = l1 + l2 - 2*i1 - 2*i2;

            for (int m = 0; m <= j/2; ++m) {
                const double t1 =
                    std::pow(-1.0, m) *
                    eri::math::factorial(j) *
                    std::pow(a1, l2 - i1 - 2*i2 - m) *
                    std::pow(a2, l1 - 2*i1 - i2 - m) /
                    std::pow(4.0, i1 + i2 + m) /
                    eri::math::factorial(i1) / eri::math::factorial(i2) / eri::math::factorial(m);

                const double t2 =
                    std::pow(gamma, 2*(i1+i2) + m) *
                    std::pow(x, j - 2*m) /
                    eri::math::factorial(l1 - 2*i1) /
                    eri::math::factorial(l2 - 2*i2) /
                    eri::math::factorial(j - 2*m);

                sm += t1 * t2;
            }
        }
    }

    return pre * sm;
}

inline double overlap_primitive_hellsing(
    int lx1, int ly1, int lz1, const std::array<double,3>& A, double a,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double b) {
    const double dx = A[0] - B[0];
    const double dy = A[1] - B[1];
    const double dz = A[2] - B[2];
    const double AB2 = dx*dx + dy*dy + dz*dz;

    const double gamma = a + b;

    const double pre = std::pow(M_PI / gamma, 1.5) * std::exp(-a * b * AB2 / gamma);

    const double wx = overlap_1d_hellsing(lx1, lx2, a, b, A[0] - B[0], gamma);
    const double wy = overlap_1d_hellsing(ly1, ly2, a, b, A[1] - B[1], gamma);
    const double wz = overlap_1d_hellsing(lz1, lz2, a, b, A[2] - B[2], gamma);

    return pre * wx * wy * wz;
}

} // namespace eri::detail