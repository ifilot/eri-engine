#pragma once

#include <array>
#include <stdexcept>

namespace eri::math {

inline double double_factorial(int n) noexcept {
    // covers (-1)!! and 0!!
    if (n <= 0) {
        return 1.0;
    }

    // Precomputed values for n!! up to 21!!
    // This safely covers angular momenta up to l = 11
    // since we need (2l - 1)!!
    static constexpr std::array<double, 22> table = {
        /*  0 */ 1.0,
        /*  1 */ 1.0,
        /*  2 */ 2.0,
        /*  3 */ 3.0,
        /*  4 */ 8.0,
        /*  5 */ 15.0,
        /*  6 */ 48.0,
        /*  7 */ 105.0,
        /*  8 */ 384.0,
        /*  9 */ 945.0,
        /* 10 */ 3840.0,
        /* 11 */ 10395.0,
        /* 12 */ 46080.0,
        /* 13 */ 135135.0,
        /* 14 */ 645120.0,
        /* 15 */ 2027025.0,
        /* 16 */ 10321920.0,
        /* 17 */ 34459425.0,
        /* 18 */ 185794560.0,
        /* 19 */ 654729075.0,
        /* 20 */ 3715891200.0,
        /* 21 */ 13749310575.0
    };

    if (n < static_cast<int>(table.size())) {
        return table[n];
    }

    // Fallback for larger n (rare in practice)
    double result = 1.0;
    for (int i = n; i > 0; i -= 2) {
        result *= static_cast<double>(i);
    }

    return result;
}

} // namespace eri::math
