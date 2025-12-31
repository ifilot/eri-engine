#pragma once

namespace eri::math {

inline double ipow(double base, int exp) {
    if (exp < 0)
        return 1.0 / ipow(base, -exp);

    double result = 1.0;
    while (exp > 0) {
        if (exp & 1)
            result *= base;
        base *= base;
        exp >>= 1;
    }
    return result;
}

} // namespace eri::math
