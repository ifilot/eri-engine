#pragma once

namespace eri::math {

/**
 * @brief  Integer-based exponentation
 * @note   Only evaluates positive powers
 * @param  base: base
 * @param  exp: power value
 * @retval 
 */
inline double ipow(double base, int exp) {
    double result = 1.0;
    while (exp > 0) {
        if (exp & 1)
            result *= base;
        base *= base;
        exp >>= 1;
    }
    return result;
}

/**
 * @brief  Integer-based exponentation
 * @note   Also accepts negative powers
 * @param  base: base
 * @param  exp: power value
 * @retval 
 */
inline double inpow(double base, int exp) {
    if (exp < 0)
        return 1.0 / ipow(base, -exp);
    return ipow(base, exp);
}

} // namespace eri::math
