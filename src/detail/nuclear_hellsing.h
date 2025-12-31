#pragma once

#include <array>
#include <vector>
#include <cmath>

#include "math/gaussian_product.h"
#include "math/fgamma.h"
#include "math/factorial.h"
#include "math/ipow.h"
#include "math/sign_pow.h"

namespace eri::detail {

inline std::pair<
    std::vector<double>,
    std::vector<std::pair<int,int>>
> A_array_hellsing(
    int l1, int l2,
    double alpha1, double alpha2,
    double x,
    double gamma,
    double pcx) {

    std::vector<double> arr;
    std::vector<std::pair<int,int>> mu_u;

    const double pre = eri::math::sign_pow(l1 + l2) * eri::math::factorial(l1) * eri::math::factorial(l2);

    for (int i1 = 0; i1 <= l1/2; ++i1) {
        for (int i2 = 0; i2 <= l2/2; ++i2) {
            for (int o1 = 0; o1 <= l1 - 2*i1; ++o1) {
                for (int o2 = 0; o2 <= l2 - 2*i2; ++o2) {
                    for (int r = 0; r <= (o1 + o2)/2; ++r) {

                        const double t1 =
                            eri::math::sign_pow(o2 + r) *
                            eri::math::factorial(o1 + o2) /
                            ( eri::math::ipow(4.0, i1 + i2 + r) *
                              eri::math::factorial(i1) * eri::math::factorial(i2) *
                              eri::math::factorial(o1) * eri::math::factorial(o2) * eri::math::factorial(r) );

                        const double t2 =
                            eri::math::ipow(alpha1, o2 - i1 - r) *
                            eri::math::ipow(alpha2, o1 - i2 - r) *
                            eri::math::ipow(x, o1 + o2 - 2*r) /
                            ( eri::math::factorial(l1 - 2*i1 - o1) *
                              eri::math::factorial(l2 - 2*i2 - o2) *
                              eri::math::factorial(o1 + o2 - 2*r) );

                        const int mu_x =
                            l1 + l2 - 2*(i1 + i2) - (o1 + o2);

                        for (int u = 0; u <= mu_x/2; ++u) {
                            const double t3 =
                                eri::math::sign_pow(u) *
                                eri::math::factorial(mu_x) *
                                eri::math::ipow(pcx, mu_x - 2*u) /
                                ( eri::math::ipow(4.0, u) *
                                  eri::math::factorial(u) *
                                  eri::math::factorial(mu_x - 2*u) *
                                  eri::math::ipow(gamma, o1 + o2 - r + u) );

                            arr.push_back(pre * t1 * t2 * t3);
                            mu_u.emplace_back(mu_x, u);
                        }
                    }
                }
            }
        }
    }

    return {arr, mu_u};
}

inline double nuclear_primitive_hellsing(
    int lx1, int ly1, int lz1, const std::array<double,3>& A, double alpha1,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double alpha2,
    const std::array<double,3>& C) {

    const double gamma = alpha1 + alpha2;
    const auto P = eri::math::gaussian_product_center(alpha1, A, alpha2, B);

    const double dx = A[0] - B[0];
    const double dy = A[1] - B[1];
    const double dz = A[2] - B[2];
    const double rab2 = dx*dx + dy*dy + dz*dz;

    const double cpx = P[0] - C[0];
    const double cpy = P[1] - C[1];
    const double cpz = P[2] - C[2];
    const double rcp2 = cpx*cpx + cpy*cpy + cpz*cpz;

    // Build A-arrays
    auto [Ax, mx] = A_array_hellsing(lx1, lx2, alpha1, alpha2, dx, gamma, cpx);
    auto [Ay, my] = A_array_hellsing(ly1, ly2, alpha1, alpha2, dy, gamma, cpy);
    auto [Az, mz] = A_array_hellsing(lz1, lz2, alpha1, alpha2, dz, gamma, cpz);

    // Maximum Boys order
    const int nu_max =
        mx[0].first + my[0].first + mz[0].first -
        (mx[0].second + my[0].second + mz[0].second);

    std::vector<double> F(nu_max + 1);
    for (int n = 0; n <= nu_max; ++n)
        F[n] = eri::math::Fgamma(n, gamma * rcp2);

    double s = 0.0;

    for (std::size_t i = 0; i < Ax.size(); ++i)
        for (std::size_t j = 0; j < Ay.size(); ++j)
            for (std::size_t k = 0; k < Az.size(); ++k) {
                const int nu =
                    mx[i].first + my[j].first + mz[k].first -
                    (mx[i].second + my[j].second + mz[k].second);

                s += Ax[i] * Ay[j] * Az[k] * F[nu];
            }

    const double pref = -2.0 * M_PI / gamma * std::exp(-alpha1 * alpha2 * rab2 / gamma);

    return pref * s;
}

} // namespace eri::detail
