#pragma once

#include <cmath>
#include <algorithm>

namespace eri::math {

inline double binom(int n, int k) noexcept {
    if (k < 0 || k > n) return 0.0;
    k = std::min(k, n - k);
    double r = 1.0;
    for (int i = 1; i <= k; ++i) {
        r *= (n - (k - i));
        r /= i;
    }
    return r;
}

inline double binomial_prefactor(int k, int l1, int l2, double x1, double x2) noexcept {
    // sum_{t=0..k} C(l1,t) C(l2,k-t) x1^(l1-t) x2^(l2-(k-t))
    double s = 0.0;
    const int tmin = std::max(0, k - l2);
    const int tmax = std::min(k, l1);

    for (int t = tmin; t <= tmax; ++t) {
        s += binom(l1, t) * binom(l2, k - t) *
             std::pow(x1, l1 - t) *
             std::pow(x2, l2 - (k - t));
    }
    return s;
}

} // namespace eri::math
