#pragma once

#include <array>
#include <cmath>
#include <eri/cog/hellsing_cache.h>

#include "math/gaussian_product.h"
#include "math/fgamma_interpolation.h"
#include "math/ipow.h"

namespace eri::detail {

inline double eri_primitive_hellsing_cached(const eri::cog::HellsingCacheTable& kernels,

    int lx1, int ly1, int lz1, const std::array<double,3>& A, double a1,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double a2,
    int lx3, int ly3, int lz3, const std::array<double,3>& C, double a3,
    int lx4, int ly4, int lz4, const std::array<double,3>& D, double a4)
{
    const auto P = eri::math::gaussian_product_center(a1,A,a2,B);
    const auto Q = eri::math::gaussian_product_center(a3,C,a4,D);

    const double g1 = a1 + a2;
    const double g2 = a3 + a4;
    const double eta = g1*g2/(g1+g2);

    const double dx = P[0]-Q[0];
    const double dy = P[1]-Q[1];
    const double dz = P[2]-Q[2];

    const double T = eta * (dx*dx + dy*dy + dz*dz);

    const int nu_max =
        lx1+ly1+lz1 + lx2+ly2+lz2 +
        lx3+ly3+lz3 + lx4+ly4+lz4;

    std::vector<double> F(nu_max+1);
    eri::math::Fgamma_block_interp(nu_max, T, F.data());

    const auto& Kx = kernels.get(lx1,lx2,lx3,lx4);
    const auto& Ky = kernels.get(ly1,ly2,ly3,ly4);
    const auto& Kz = kernels.get(lz1,lz2,lz3,lz4);

    const double AxBx = A[0]-B[0];
    const double CxDx = C[0]-D[0];
    const double AyBy = A[1]-B[1];
    const double CyDy = C[1]-D[1];
    const double AzBz = A[2]-B[2];
    const double CzDz = C[2]-D[2];

    double sx = 0.0, sy = 0.0, sz = 0.0;

    auto eval_1d = [&](const eri::cog::HellsingCache1D& K, double AB, double CD, double PQ) {
        double s = 0.0;
        for(std::size_t i=0;i<K.scalar.size();++i)
        {
            const auto& p = K.powers[i];
            s += K.scalar[i] *
                 eri::math::ipow(a1,p[0]) *
                 eri::math::ipow(a2,p[1]) *
                 eri::math::ipow(a3,p[2]) *
                 eri::math::ipow(a4,p[3]) *
                 eri::math::ipow(g1,p[4]) *
                 eri::math::ipow(g2,p[5]) *
                 eri::math::ipow(AB,p[6]) *
                 eri::math::ipow(CD,p[7]) *
                 eri::math::ipow(eta,p[8]) *
                 eri::math::ipow(PQ,p[9]);
        }
        return s;
    };

    sx = eval_1d(Kx, AxBx, CxDx, dx);
    sy = eval_1d(Ky, AyBy, CyDy, dy);
    sz = eval_1d(Kz, AzBz, CzDz, dz);

    constexpr double pi25 = (M_PI*M_PI)*std::sqrt(M_PI);
    const double pref = 2.0*pi25/(g1*g2*std::sqrt(g1+g2));

    return pref * sx * sy * sz;
}

} // namespace eri::detail
