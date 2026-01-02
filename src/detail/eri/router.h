#pragma once

#include <array>
#include <cmath>
#include <eri/cog/hellsing_cache.h>

#include "math/gaussian_product.h"
#include "math/fgamma_interpolation.h"
#include "math/ipow.h"

namespace eri::detail {

inline double eri_primitive_router(
    const eri::cog::HellsingCacheTable1D& kernels,
    int lx1, int ly1, int lz1, const std::array<double,3>& A, double a1,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double a2,
    int lx3, int ly3, int lz3, const std::array<double,3>& C, double a3,
    int lx4, int ly4, int lz4, const std::array<double,3>& D, double a4)
{
    // calculate angular momenta
    const int L1 = lx1+ly1+lz1;
    const int L2 = lx2+ly2+lz2;
    const int L3 = lx3+ly3+lz3;
    const int L4 = lx4+ly4+lz4;
    const int nu_max = L1 + L2 + L3 + L4;

    // ---------------- geometry ----------------
    const auto P = eri::math::gaussian_product_center(a1,A,a2,B);
    const auto Q = eri::math::gaussian_product_center(a3,C,a4,D);

    const double dxAB = A[0]-B[0];
    const double dyAB = A[1]-B[1];
    const double dzAB = A[2]-B[2];
    const double rab2 = dxAB*dxAB + dyAB*dyAB + dzAB*dzAB;

    const double dxCD = C[0]-D[0];
    const double dyCD = C[1]-D[1];
    const double dzCD = C[2]-D[2];
    const double rcd2 = dxCD*dxCD + dyCD*dyCD + dzCD*dzCD;

    const double dx = P[0]-Q[0];
    const double dy = P[1]-Q[1];
    const double dz = P[2]-Q[2];
    const double rpq2 = dx*dx + dy*dy + dz*dz;

    // ---------------- scalars ----------------
    const double g1 = a1 + a2;
    const double g2 = a3 + a4;
    const double eta = g1*g2/(g1+g2);

    // ---------------- prefactor ----------------
    constexpr double pi25 = (M_PI*M_PI)*std::sqrt(M_PI);
    const double pref = 2.0*pi25/(g1*g2*std::sqrt(g1+g2)) *
                        std::exp(-a1*a2*rab2/g1) *
                        std::exp(-a3*a4*rcd2/g2);

    // evaluate (ss|ss)
    if(nu_max == 0) {
        return pref * eri::math::Fgamma_interp(nu_max, eta*rpq2);
    }

    // ---------------- Fgamma ----------------
    std::vector<double> F(nu_max+1);
    eri::math::Fgamma_block_interp(nu_max, eta*rpq2, F.data());

    // fall through to Hellsing precomputed kernels
    const auto& Kx = kernels.get(lx1,lx2,lx3,lx4);
    const auto& Ky = kernels.get(ly1,ly2,ly3,ly4);
    const auto& Kz = kernels.get(lz1,lz2,lz3,lz4);

    // ---------------- precompute 1D polynomials ----------------
    auto build_poly = [&](const eri::cog::HellsingCache1D& K,
                          double AB, double CD, double PQ,
                          std::vector<double>& poly,
                          std::vector<int>& nu)
    {
        const std::size_t n = K.scalar.size();
        poly.resize(n);
        nu.resize(n);

        for (std::size_t i=0; i<n; ++i) {
            const auto& p = K.powers[i];
            poly[i] =
                K.scalar[i] *
                eri::math::inpow(a1,p[0]) *
                eri::math::inpow(a2,p[1]) *
                eri::math::inpow(a3,p[2]) *
                eri::math::inpow(a4,p[3]) *
                eri::math::inpow(g1,p[4]) *
                eri::math::inpow(g2,p[5]) *
                eri::math::inpow(AB,p[6]) *
                eri::math::inpow(CD,p[7]) *
                eri::math::inpow(eta,p[8]) *
                eri::math::inpow(PQ,p[9]);

            nu[i] = K.mu[i] - K.u[i];
        }
    };

    std::vector<double> polx, poly, polz;
    std::vector<int>    nux,  nuy,  nuz;

    build_poly(Kx, A[0]-B[0], C[0]-D[0], dx, polx, nux);
    build_poly(Ky, A[1]-B[1], C[1]-D[1], dy, poly, nuy);
    build_poly(Kz, A[2]-B[2], C[2]-D[2], dz, polz, nuz);

    // ---------------- triple contraction ----------------
    double s = 0.0;
    for (std::size_t i=0; i<polx.size(); ++i) {
        const double pi = polx[i];
        const int nui = nux[i];
        for (std::size_t j=0; j<poly.size(); ++j) {
            const double pij = pi * poly[j];
            const int nuij = nui + nuy[j];
            for (std::size_t k=0; k<polz.size(); ++k) {
                s += pij * polz[k] * F[nuij + nuz[k]];
            }
        }
    }

    return pref * s;
}

} // namespace eri::detail
