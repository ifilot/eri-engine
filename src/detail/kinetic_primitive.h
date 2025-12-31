#pragma once

#include <array>

#include "detail/overlap_primitive.h"

namespace eri::detail {

inline double kinetic_primitive(
    int lx1, int ly1, int lz1, const std::array<double,3>& A, double a,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double b,
    eri::enums::KineticMethod method) {
    const int L2 = lx2 + ly2 + lz2;

    // Explicit, safe mapping
    const eri::enums::OverlapMethod overlap_method =
        (method == eri::enums::KineticMethod::Huzinaga)
            ? eri::enums::OverlapMethod::Huzinaga
            : eri::enums::OverlapMethod::Hellsing;

    const auto S = [&](int dx, int dy, int dz) -> double {
        const int nx = lx2 + dx, ny = ly2 + dy, nz = lz2 + dz;
        if (nx < 0 || ny < 0 || nz < 0) return 0.0;
        return overlap_primitive(
            lx1, ly1, lz1, A, a,
            nx,  ny,  nz,  B, b,
            overlap_method
        );
    };

    const double t0 = b * (2.0 * L2 + 3.0) * S(0,0,0);

    const double t1 = -2.0 * b * b * (S(2,0,0) + S(0,2,0) + S(0,0,2));

    const double t2 = -0.5 * (
        (lx2 >= 2 ? lx2 * (lx2 - 1) * S(-2,0,0) : 0.0) +
        (ly2 >= 2 ? ly2 * (ly2 - 1) * S(0,-2,0) : 0.0) +
        (lz2 >= 2 ? lz2 * (lz2 - 1) * S(0,0,-2) : 0.0)
    );

    return t0 + t1 + t2;
}

} // namespace eri::detail