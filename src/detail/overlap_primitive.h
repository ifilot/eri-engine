#pragma once

#include <eri/enums.h>
#include <array>

#include "detail/overlap_hellsing.h"
#include "detail/overlap_huzinaga.h"

namespace eri::detail {

inline double overlap_primitive(
    int lx1, int ly1, int lz1, const std::array<double,3>& A, double a,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double b,
    eri::enums::OverlapMethod method) {
    switch (method) {
        case eri::enums::OverlapMethod::Huzinaga:
            return eri::detail::overlap_primitive_huzinaga(
                lx1, ly1, lz1, A, a,
                lx2, ly2, lz2, B, b
            );

        case eri::enums::OverlapMethod::Hellsing:
            return eri::detail::overlap_primitive_hellsing(
                lx1, ly1, lz1, A, a,
                lx2, ly2, lz2, B, b
            );

        default:
            __builtin_unreachable();
    }
}

} // namespace eri::detail