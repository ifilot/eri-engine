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

// Each entry corresponds to one term in the Hellsing B-array expansion
struct HellsingBTerm {
    double c;  // coefficient
    int mu;    // mu
    int u;     // u
};

inline std::vector<HellsingBTerm> B_array_hellsing(
    int l1, int l2, int l3, int l4,
    double a1, double a2, double a3, double a4,
    double ax, double bx, double cx, double dx,
    double px, double qx,
    double g1, double g2)
{
    // pre1 = (-1)^(l1+l2) * l1! * l2! / g1^(l1+l2)
    // pre2 = l3! * l4! / g2^(l3+l4)
    const double pre1 = eri::math::sign_pow(l1 + l2)
                      * eri::math::factorial(l1) * eri::math::factorial(l2)
                      / eri::math::ipow(g1, l1 + l2);

    const double pre2 = eri::math::factorial(l3) * eri::math::factorial(l4)
                      / eri::math::ipow(g2, l3 + l4);

    const double eta = (g1 * g2) / (g1 + g2);

    std::vector<HellsingBTerm> arr;
    // Reserve is optional; keeping it modest avoids over-reserving for small l.
    // arr.reserve(...);

    // Python:
    // for i1 in range(l1//2+1):
    for (int i1 = 0; i1 <= l1 / 2; ++i1) {
        for (int i2 = 0; i2 <= l2 / 2; ++i2) {
            for (int o1 = 0; o1 <= l1 - 2 * i1; ++o1) {
                for (int o2 = 0; o2 <= l2 - 2 * i2; ++o2) {
                    const int oo12 = o1 + o2;
                    for (int r1 = 0; r1 <= oo12 / 2; ++r1) {
                        // t11 = (-1)^(o2+r1) * (o1+o2)! / (4^(i1+i2+r1) i1! i2! o1! o2! r1!)
                        const double t11 = eri::math::sign_pow(o2 + r1)
                                         * eri::math::factorial(oo12)
                                         / eri::math::ipow(4.0, i1 + i2 + r1)
                                         / eri::math::factorial(i1)
                                         / eri::math::factorial(i2)
                                         / eri::math::factorial(o1)
                                         / eri::math::factorial(o2)
                                         / eri::math::factorial(r1);

                        // t12 = a1^(o2-i1-r1) * a2^(o1-i2-r1) * g1^(2*(i1+i2)+r1)
                        //       * (ax-bx)^(o1+o2-2*r1) /
                        //       ( (l1-2*i1-o1)! (l2-2*i2-o2)! (o1+o2-2*r1)! )
                        const int e_a1 = o2 - i1 - r1;
                        const int e_a2 = o1 - i2 - r1;
                        const int e_g1 = 2 * (i1 + i2) + r1;
                        const int e_x1 = oo12 - 2 * r1;

                        // Negative exponents should not occur in valid ranges, but keep robust:
                        if (e_a1 < 0 || e_a2 < 0 || e_x1 < 0) continue;

                        const double num12 =
                            eri::math::ipow(a1, e_a1) *
                            eri::math::ipow(a2, e_a2) *
                            eri::math::ipow(g1, e_g1) *
                            eri::math::ipow(ax - bx, e_x1);

                        const int f1 = l1 - 2 * i1 - o1;
                        const int f2 = l2 - 2 * i2 - o2;
                        const int f3 = e_x1;

                        if (f1 < 0 || f2 < 0 || f3 < 0) continue;

                        const double den12 =
                            eri::math::factorial(f1) *
                            eri::math::factorial(f2) *
                            eri::math::factorial(f3);

                        const double t12 = num12 / den12;

                        for (int i3 = 0; i3 <= l3 / 2; ++i3) {
                            for (int i4 = 0; i4 <= l4 / 2; ++i4) {
                                for (int o3 = 0; o3 <= l3 - 2 * i3; ++o3) {
                                    for (int o4 = 0; o4 <= l4 - 2 * i4; ++o4) {
                                        const int oo34 = o3 + o4;
                                        for (int r2 = 0; r2 <= oo34 / 2; ++r2) {
                                            // t21 = (-1)^(o3+r2) * (o3+o4)! / (4^(i3+i4+r2) i3! i4! o3! o4! r2!)
                                            const double t21 = eri::math::sign_pow(o3 + r2)
                                                             * eri::math::factorial(oo34)
                                                             / eri::math::ipow(4.0, i3 + i4 + r2)
                                                             / eri::math::factorial(i3)
                                                             / eri::math::factorial(i4)
                                                             / eri::math::factorial(o3)
                                                             / eri::math::factorial(o4)
                                                             / eri::math::factorial(r2);

                                            // t22 = a3^(o4-i3-r2) * a4^(o3-i4-r2) * g2^(2*(i3+i4)+r2)
                                            //       * (cx-dx)^(o3+o4-2*r2) /
                                            //       ( (l3-2*i3-o3)! (l4-2*i4-o4)! (o3+o4-2*r2)! )
                                            const int e_a3 = o4 - i3 - r2;
                                            const int e_a4 = o3 - i4 - r2;
                                            const int e_g2 = 2 * (i3 + i4) + r2;
                                            const int e_x2 = oo34 - 2 * r2;

                                            if (e_a3 < 0 || e_a4 < 0 || e_x2 < 0) continue;

                                            const double num22 =
                                                eri::math::ipow(a3, e_a3) *
                                                eri::math::ipow(a4, e_a4) *
                                                eri::math::ipow(g2, e_g2) *
                                                eri::math::ipow(cx - dx, e_x2);

                                            const int g1f = l3 - 2 * i3 - o3;
                                            const int g2f = l4 - 2 * i4 - o4;
                                            const int g3f = e_x2;

                                            if (g1f < 0 || g2f < 0 || g3f < 0) continue;

                                            const double den22 =
                                                eri::math::factorial(g1f) *
                                                eri::math::factorial(g2f) *
                                                eri::math::factorial(g3f);

                                            const double t22 = num22 / den22;

                                            // mu = l1+l2+l3+l4 - 2*(i1+i2+i3+i4) - (o1+o2+o3+o4)
                                            const int mu =
                                                (l1 + l2 + l3 + l4)
                                                - 2 * (i1 + i2 + i3 + i4)
                                                - (o1 + o2 + o3 + o4);

                                            if (mu < 0) continue;

                                            // for u in range(mu//2+1):
                                            for (int u = 0; u <= mu / 2; ++u) {
                                                // t3 = (-1)^u * mu! * eta^(mu-u) * (px-qx)^(mu-2u) / (4^u u! (mu-2u)!)
                                                const int e_eta = mu - u;
                                                const int e_pq  = mu - 2 * u;
                                                if (e_pq < 0) continue;

                                                const double t3 =
                                                    eri::math::sign_pow(u) *
                                                    eri::math::factorial(mu) *
                                                    eri::math::ipow(eta, e_eta) *
                                                    eri::math::ipow(px - qx, e_pq) /
                                                    eri::math::ipow(4.0, u) /
                                                    eri::math::factorial(u) /
                                                    eri::math::factorial(e_pq);

                                                const double coeff = pre1 * pre2 * t11 * t12 * t21 * t22 * t3;
                                                arr.push_back(HellsingBTerm{coeff, mu, u});
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return arr;
}


inline double eri_primitive_hellsing(
    int lx1, int ly1, int lz1, const std::array<double,3>& A, double a1,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double a2,
    int lx3, int ly3, int lz3, const std::array<double,3>& C, double a3,
    int lx4, int ly4, int lz4, const std::array<double,3>& D, double a4)
{
    const double dxAB = A[0] - B[0];
    const double dyAB = A[1] - B[1];
    const double dzAB = A[2] - B[2];
    const double rab2 = dxAB*dxAB + dyAB*dyAB + dzAB*dzAB;

    const double dxCD = C[0] - D[0];
    const double dyCD = C[1] - D[1];
    const double dzCD = C[2] - D[2];
    const double rcd2 = dxCD*dxCD + dyCD*dyCD + dzCD*dzCD;

    const auto P = eri::math::gaussian_product_center(a1, A, a2, B);
    const auto Q = eri::math::gaussian_product_center(a3, C, a4, D);

    const double dxPQ = P[0] - Q[0];
    const double dyPQ = P[1] - Q[1];
    const double dzPQ = P[2] - Q[2];
    const double rpq2 = dxPQ*dxPQ + dyPQ*dyPQ + dzPQ*dzPQ;

    const double gamma1 = a1 + a2;
    const double gamma2 = a3 + a4;
    const double eta    = (gamma1 * gamma2) / (gamma1 + gamma2);

    // Build B-arrays for each Cartesian component (x, y, z)
    const auto Bx = B_array_hellsing(
        lx1, lx2, lx3, lx4,
        a1, a2, a3, a4,
        A[0], B[0], C[0], D[0],
        P[0], Q[0],
        gamma1, gamma2);

    const auto By = B_array_hellsing(
        ly1, ly2, ly3, ly4,
        a1, a2, a3, a4,
        A[1], B[1], C[1], D[1],
        P[1], Q[1],
        gamma1, gamma2);

    const auto Bz = B_array_hellsing(
        lz1, lz2, lz3, lz4,
        a1, a2, a3, a4,
        A[2], B[2], C[2], D[2],
        P[2], Q[2],
        gamma1, gamma2);

    // Precompute Fgamma values up to nu_max
    const int nu_max =
        (lx1 + ly1 + lz1) +
        (lx2 + ly2 + lz2) +
        (lx3 + ly3 + lz3) +
        (lx4 + ly4 + lz4);

    std::vector<double> F(nu_max + 1);
    const double T = eta * rpq2;
    for (int nu = 0; nu <= nu_max; ++nu) {
        F[nu] = eri::math::Fgamma(nu, T);
    }

    // Triple sum:
    // nu = (mu_x + mu_y + mu_z) - (u_x + u_y + u_z)
    double s = 0.0;
    for (const auto& tx : Bx) {
        for (const auto& ty : By) {
            const int mu_xy = tx.mu + ty.mu;
            const int u_xy  = tx.u  + ty.u;
            const double c_xy = tx.c * ty.c;

            for (const auto& tz : Bz) {
                const int nu = (mu_xy + tz.mu) - (u_xy + tz.u);
                // In theory nu is in [0, nu_max] given the construction; keep safe:
                if (nu >= 0 && nu <= nu_max) {
                    s += c_xy * tz.c * F[nu];
                }
            }
        }
    }

    // Universal pre-factor (same as Huzinaga version, but using eta in Fgamma arg)
    constexpr double pi25 = (M_PI * M_PI) * std::sqrt(M_PI);
    const double pref =
        2.0 * pi25 / (gamma1 * gamma2 * std::sqrt(gamma1 + gamma2)) *
        std::exp(-a1 * a2 * rab2 / gamma1) *
        std::exp(-a3 * a4 * rcd2 / gamma2);

    return pref * s;
}

} // namespace eri::detail
