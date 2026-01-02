#pragma once

#include "detail/eri/router.h"
#include "detail/eri/hellsing_cached.h"

namespace eri::two_electron::primitives {

struct PrimitiveRouter {
    const eri::cog::HellsingCacheTable1D& cache;

    explicit PrimitiveRouter(
        const eri::cog::HellsingCacheTable1D& c
    ) : cache(c) {}

    double operator()(
        int lx1, int ly1, int lz1, const std::array<double,3>& A, double a1,
        int lx2, int ly2, int lz2, const std::array<double,3>& B, double a2,
        int lx3, int ly3, int lz3, const std::array<double,3>& C, double a3,
        int lx4, int ly4, int lz4, const std::array<double,3>& D, double a4
    ) const {
        return eri::detail::eri_primitive_router(
            cache,
            lx1, ly1, lz1, A, a1,
            lx2, ly2, lz2, B, a2,
            lx3, ly3, lz3, C, a3,
            lx4, ly4, lz4, D, a4
        );
    }
};

}