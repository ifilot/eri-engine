#pragma once

#include <array>
#include <vector>
#include <cmath>

#include "math/gaussian_product.h"
#include "math/binomial_prefactor.h"
#include "math/fgamma.h"
#include "math/factorial.h"

namespace eri::detail {

inline double A_term_huzinaga(
    int i, int r, int u,
    int l1, int l2,
    double pa, double pb,
    double cp,
    double gamma)
{
    const int m = i - 2*r - 2*u;   // power of cp
    if (m < 0) return 0.0;

    const double binom = eri::math::binomial_prefactor(i, l1, l2, pa, pb);

    const double sign_i = (i % 2 == 0) ? 1.0 : -1.0;
    const double sign_u = (u % 2 == 0) ? 1.0 : -1.0;

    const double fi = eri::math::factorial(i);
    const double fr = eri::math::factorial(r);
    const double fu = eri::math::factorial(u);
    const double fm = eri::math::factorial(m);

    const double cp_pow = (m == 0) ? 1.0 : std::pow(cp, m);
    const double g_pow  = std::pow(0.25 / gamma, r + u);

    return sign_i * binom * sign_u * fi * cp_pow * g_pow / (fr * fu * fm);
}

inline std::vector<double> A_array_huzinaga(
    int l1, int l2,
    double pa, double pb,
    double cp,
    double gamma)
{
    const int imax = l1 + l2 + 1;
    std::vector<double> arr(imax, 0.0);

    for (int i = 0; i < imax; ++i) {
        for (int r = 0; r <= i/2; ++r) {
            for (int u = 0; u <= (i - 2*r)/2; ++u) {
                const int iI = i - 2*r - u;
                if (iI >= 0 && iI < imax) {
                    arr[iI] += A_term_huzinaga(i, r, u, l1, l2, pa, pb, cp, gamma);
                }
            }
        }
    }

    return arr;
}

inline double nuclear_primitive_huzinaga(
    int lx1, int ly1, int lz1, const std::array<double,3>& A, double a,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double b,
    const std::array<double,3>& C)
{
    const double dx = A[0] - B[0];
    const double dy = A[1] - B[1];
    const double dz = A[2] - B[2];
    const double AB2 = dx*dx + dy*dy + dz*dz;

    const double gamma = a + b;
    const auto P = eri::math::gaussian_product_center(a, A, b, B);

    const double cpx = P[0] - C[0];
    const double cpy = P[1] - C[1];
    const double cpz = P[2] - C[2];
    const double RCP2 = cpx*cpx + cpy*cpy + cpz*cpz;

    // Build A arrays
    const auto Ax = A_array_huzinaga(lx1, lx2, P[0] - A[0], P[0] - B[0], P[0] - C[0], gamma);
    const auto Ay = A_array_huzinaga(ly1, ly2, P[1] - A[1], P[1] - B[1], P[1] - C[1], gamma);
    const auto Az = A_array_huzinaga(lz1, lz2, P[2] - A[2], P[2] - B[2], P[2] - C[2], gamma);

    // Boys values up to nu_max
    const int nu_max = (lx1+ly1+lz1) + (lx2+ly2+lz2);
    std::vector<double> F(nu_max + 1);
    const double T = gamma * RCP2;
    for (int n = 0; n <= nu_max; ++n) {
        F[n] = eri::math::Fgamma(n, T);
    }

    // Triple sum
    double s = 0.0;
    for (int i = 0; i < (int)Ax.size(); ++i)
        for (int j = 0; j < (int)Ay.size(); ++j)
            for (int k = 0; k < (int)Az.size(); ++k) {
                const int n = i + j + k;
                if (n <= nu_max)
                    s += Ax[i] * Ay[j] * Az[k] * F[n];
            }

    const double pref = -2.0 * M_PI / gamma * std::exp(-a * b * AB2 / gamma);
    return pref * s;
}

} // namespace eri::detail
