#pragma once

#include <array>
#include <eri/enums.h>

#include "detail/nuclear_huzinaga.h"
#include "detail/nuclear_hellsing.h"

namespace eri::detail {

inline double nuclear_primitive(
    int lx1, int ly1, int lz1, const std::array<double,3>& A, double a,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double b,
    const std::array<double,3>& C,
    eri::enums::NuclearMethod method) {
    switch (method) {
        case eri::enums::NuclearMethod::Huzinaga:
            return nuclear_primitive_huzinaga(
                lx1, ly1, lz1, A, a,
                lx2, ly2, lz2, B, b,
                C
            );
        case eri::enums::NuclearMethod::Hellsing:
            // placeholder until implemented
            return nuclear_primitive_hellsing(
                lx1, ly1, lz1, A, a,
                lx2, ly2, lz2, B, b,
                C
            );
        default:
            __builtin_unreachable();
    }
}

} // namespace eri::detail