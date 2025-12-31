#pragma once

#include <cmath>

namespace eri::math {

/**
 * @brief Factorial n! as double, optimized for small n.
 *
 * - Exact for 0 <= n <= 19
 * - Returns 1.0 for n < 0 (useful for angular-momentum formulas)
 * - Falls back to tgamma for larger n (rare in practice)
 *
 * Intended for Gaussian integral kernels (l typically <= 6).
 */
inline double factorial(int n) noexcept {
    static constexpr double table[] = {
        1.0,                        // 0!
        1.0,                        // 1!
        2.0,                        // 2!
        6.0,                        // 3!
        24.0,                       // 4!
        120.0,                      // 5!
        720.0,                      // 6!
        5040.0,                     // 7!
        40320.0,                    // 8!
        362880.0,                   // 9!
        3628800.0,                  // 10!
        39916800.0,                 // 11!
        479001600.0,                // 12!
        6227020800.0,               // 13!
        87178291200.0,              // 14!
        1307674368000.0,            // 15!
        20922789888000.0,           // 16!
        355687428096000.0,          // 17!
        6402373705728000.0,         // 18!
        121645100408832000.0        // 19!
    };

    if (n < 0)
        return 1.0;

    if (n < static_cast<int>(sizeof(table) / sizeof(table[0])))
        return table[n];

    // Fallback: not expected for typical Gaussian integrals
    return std::tgamma(static_cast<double>(n) + 1.0);
}

} // namespace eri::math
