#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <algorithm>

#include "math/gaussian_product.h"
#include "math/binomial_prefactor.h"
#include "math/fgamma_block.h"
#include "math/factorial.h"
#include "math/ipow.h"
#include "math/sign_pow.h"

namespace eri::detail {

inline double fact_ratio2(int a, int b) {
    const int am2b = a - 2*b;
    if (a < 0 || b < 0 || am2b < 0) return 0.0;
    return eri::math::factorial(a) / (eri::math::factorial(b) * eri::math::factorial(am2b));
}

inline double B0_huzinaga(int i, int r, double g) {
    const double fr = fact_ratio2(i, r);
    if (fr == 0.0) return 0.0;
    return fr * eri::math::inpow(4.0 * g, r - i);
}

inline double fB_huzinaga(
    int i, int l1, int l2,
    double p, double a, double b,
    int r, double g)
{
    const double binom = eri::math::binomial_prefactor(i, l1, l2, p - a, p - b);
    if (binom == 0.0) return 0.0;
    return binom * B0_huzinaga(i, r, g);
}

inline double B_term_huzinaga(
    int i1, int i2, int r1, int r2, int u,
    int l1, int l2, int l3, int l4,
    double px, double ax, double bx,
    double qx, double cx, double dx,
    double gamma1, double gamma2, double delta)
{
    const int t = i1 + i2 - 2 * (r1 + r2);
    const int m = t - 2*u;
    if (m < 0) return 0.0;

    const int ed = t - u;
    if (ed < 0) return 0.0;

    const double f1 = fB_huzinaga(i1, l1, l2, px, ax, bx, r1, gamma1);
    if (f1 == 0.0) return 0.0;

    const double f2 = fB_huzinaga(i2, l3, l4, qx, cx, dx, r2, gamma2);
    if (f2 == 0.0) return 0.0;

    const double s_i2 = eri::math::sign_pow(i2);
    const double s_u  = eri::math::sign_pow(u);

    const double fr2 = fact_ratio2(t, u);
    if (fr2 == 0.0) return 0.0;

    const double qp = (qx - px);
    const double qp_pow = (m == 0) ? 1.0 : eri::math::ipow(qp, m);
    const double d_pow  = (ed == 0) ? 1.0 : eri::math::ipow(delta, ed);

    return f1 * s_i2 * f2 * s_u * fr2 * qp_pow / d_pow;
}

inline std::vector<double> B_array_huzinaga(
    int l1, int l2, int l3, int l4,
    double p, double a, double b,
    double q, double c, double d,
    double gamma1, double gamma2, double delta)
{
    const int imax = l1 + l2 + l3 + l4 + 1;
    std::vector<double> arr(imax, 0.0);

    for (int i1 = 0; i1 <= l1 + l2; ++i1) {
        for (int i2 = 0; i2 <= l3 + l4; ++i2) {
            for (int r1 = 0; r1 <= i1/2; ++r1) {
                for (int r2 = 0; r2 <= i2/2; ++r2) {
                    const int t = i1 + i2 - 2*(r1 + r2);
                    if (t < 0) continue;

                    // u in [0, floor(t/2)]
                    for (int u = 0; u <= t/2; ++u) {
                        const int i = t - u;
                        if (i >= 0 && i < imax) {
                            arr[i] += B_term_huzinaga(
                                i1, i2, r1, r2, u,
                                l1, l2, l3, l4,
                                p, a, b, q, c, d,
                                gamma1, gamma2, delta
                            );
                        }
                    }
                }
            }
        }
    }

    return arr;
}

inline double eri_primitive_huzinaga(
    int lx1, int ly1, int lz1, const std::array<double,3>& A, double a,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double b,
    int lx3, int ly3, int lz3, const std::array<double,3>& C, double c,
    int lx4, int ly4, int lz4, const std::array<double,3>& D, double d)
{
    const double dxAB = A[0] - B[0];
    const double dyAB = A[1] - B[1];
    const double dzAB = A[2] - B[2];
    const double rab2 = dxAB*dxAB + dyAB*dyAB + dzAB*dzAB;

    const double dxCD = C[0] - D[0];
    const double dyCD = C[1] - D[1];
    const double dzCD = C[2] - D[2];
    const double rcd2 = dxCD*dxCD + dyCD*dyCD + dzCD*dzCD;

    const double gamma1 = a + b;
    const double gamma2 = c + d;

    const auto P = eri::math::gaussian_product_center(a, A, b, B);
    const auto Q = eri::math::gaussian_product_center(c, C, d, D);

    const double dxPQ = P[0] - Q[0];
    const double dyPQ = P[1] - Q[1];
    const double dzPQ = P[2] - Q[2];
    const double rpq2 = dxPQ*dxPQ + dyPQ*dyPQ + dzPQ*dzPQ;

    const double delta = 0.25 * (1.0/gamma1 + 1.0/gamma2);

    // build B-arrays for each Cartesian component
    const auto Bx = B_array_huzinaga(lx1, lx2, lx3, lx4, P[0], A[0], B[0], Q[0], C[0], D[0], gamma1, gamma2, delta);
    const auto By = B_array_huzinaga(ly1, ly2, ly3, ly4, P[1], A[1], B[1], Q[1], C[1], D[1], gamma1, gamma2, delta);
    const auto Bz = B_array_huzinaga(lz1, lz2, lz3, lz4, P[2], A[2], B[2], Q[2], C[2], D[2], gamma1, gamma2, delta);

    // precompute Fgamma values up to nu_max
    const int nu_max = (lx1+ly1+lz1) + (lx2+ly2+lz2) + (lx3+ly3+lz3) + (lx4+ly4+lz4);
    std::vector<double> F(nu_max + 1);
    const double T = 0.25 * rpq2 / delta;

    // build Fgamma block
    if (nu_max == 0) {  // fast: no recurrence needed
         F[0] = eri::math::Fgamma(0, T);
    } else {
        eri::math::Fgamma_block(nu_max, T, F.data());
    }

    // evaluate triple-sum
    double s = 0.0;
    for (int i = 0; i < (int)Bx.size(); ++i) {
        for (int j = 0; j < (int)By.size(); ++j) {
            const int ij = i + j;
            for (int k = 0; k < (int)Bz.size(); ++k) {
                const int n = ij + k;
                if (n <= nu_max) {
                    s += Bx[i] * By[j] * Bz[k] * F[n];
                }
            }
        }
    }

    // universal pre-factor
    constexpr double pi25 = (M_PI * M_PI) * std::sqrt(M_PI);
    const double pref =
        2.0 * pi25 / (gamma1 * gamma2 * std::sqrt(gamma1 + gamma2)) *
        std::exp(-a * b * rab2 / gamma1) *
        std::exp(-c * d * rcd2 / gamma2);

    return pref * s;
}

} // namespace eri::detail
