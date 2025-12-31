#pragma once

#include <cmath>
#include "math/gaussian_product.h"
#include "math/binomial_prefactor.h"
#include "math/double_factorial.h"

namespace eri::detail {

inline double overlap_1d_huzinaga(int l1, int l2, double x1, double x2, double gamma) {
    double sm = 0.0;
    const int imax = (l1 + l2) / 2;

    for (int i = 0; i <= imax; ++i) {
        const double df = (i == 0) ? 1.0 : eri::math::double_factorial(2 * i - 1);
        sm += eri::math::binomial_prefactor(2*i, l1, l2, x1, x2) * df / std::pow(2 * gamma, i);
    }

    return sm;
}

inline double overlap_primitive_huzinaga(
    int lx1, int ly1, int lz1, const std::array<double,3>& A, double a,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double b)
{
    const double dx = A[0] - B[0];
    const double dy = A[1] - B[1];
    const double dz = A[2] - B[2];
    const double AB2 = dx*dx + dy*dy + dz*dz;

    const double gamma = a + b;
    const auto P = eri::math::gaussian_product_center(a, A, b, B);

    const double pre = std::pow(M_PI / gamma, 1.5) * std::exp(-a * b * AB2 / gamma);

    const double wx = overlap_1d_huzinaga(lx1, lx2, P[0] - A[0], P[0] - B[0], gamma);
    const double wy = overlap_1d_huzinaga(ly1, ly2, P[1] - A[1], P[1] - B[1], gamma);
    const double wz = overlap_1d_huzinaga(lz1, lz2, P[2] - A[2], P[2] - B[2], gamma);

    return pre * wx * wy * wz;
}

}