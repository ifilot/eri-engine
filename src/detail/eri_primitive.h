#pragma once

#include <eri/enums.h>
#include <array>

#include "detail/eri_hellsing.h"
#include "detail/eri_huzinaga.h"

namespace eri::detail {

inline double eri_primitive(
    int lx1, int ly1, int lz1, const std::array<double,3>& A, double a,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double b,
    int lx3, int ly3, int lz3, const std::array<double,3>& C, double c,
    int lx4, int ly4, int lz4, const std::array<double,3>& D, double d,
    eri::enums::ERIMethod method) {
    switch (method) {
        case eri::enums::ERIMethod::Huzinaga:
            return eri::detail::eri_primitive_huzinaga(
                lx1, ly1, lz1, A, a,
                lx2, ly2, lz2, B, b,
                lx3, ly3, lz3, C, c,
                lx4, ly4, lz4, D, d
            );

        case eri::enums::ERIMethod::Hellsing:
            return eri::detail::eri_primitive_hellsing(
                lx1, ly1, lz1, A, a,
                lx2, ly2, lz2, B, b,
                lx3, ly3, lz3, C, c,
                lx4, ly4, lz4, D, d
            );

        default:
            __builtin_unreachable();
    }
}

} // namespace eri::detail