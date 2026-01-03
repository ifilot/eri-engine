#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "math/double_factorial.h"
#include "math/ipow.h"

namespace eri::math {

namespace detail {

double approximate_error(int m, double t) {
    const double expT = std::exp(-t);

    // a0 = 1
    double a_prev = 1.0;
    double sum = 0.0;

    for (int i = 1; i < 1000; ++i) {
        // a_i = ((2m - 2i + 1) / (2t)) * a_{i-1}
        const double a = (2.0 * m - 2.0 * i + 1.0) / (2.0 * t) * a_prev;
        sum += a;

        // stop when terms start growing in magnitude
        if (std::abs(a) > std::abs(a_prev)) break;

        a_prev = a;
    }

    return expT / (2.0 * t) * sum;
}

std::pair<double, double> large_t_approximation(int m, double t) {
    const double df   = eri::math::double_factorial(2 * m - 1);
    const double pow2 = static_cast<double>(eri::math::ipow(2, m + 1));

    const double F = df / pow2 * std::sqrt(M_PI / eri::math::ipow(t, 2*m + 1));

    const double delta_F = approximate_error(m, t);

    return {delta_F, F}; // [delta, approx]
}

double small_t_approximation(int m, double t) {
    const double expT = std::exp(-t);
    
    double a = 1.0;
    double sum = 1.0;

    for (int i = 1; i < 10000; ++i) {
        a *= (2.0 * t) / (2.0 * m + 2.0 * i + 1.0);
        sum += a;

        if (std::abs(a) < 1e-16 * std::abs(sum)) {
            break;
        }
    }

    return expT * sum / (2.0 * m + 1.0);
}

} // namespace detail

double boys(int m, double t) {
    if (t == 0.0) {
        return 1.0 / (2.0 * m + 1.0);
    }

    // Small-t regime: asymptotic expansion is invalid
    if (t < m + 0.5) {
        return eri::math::detail::small_t_approximation(m, t);
    }

    // Try large-t asymptotic approximation
    auto [delta, approx] =
        eri::math::detail::large_t_approximation(m, t);

    // Accept asymptotic if error estimate is sufficiently small
    if (delta < 1e-15 * std::abs(approx)) {
        return approx;
    }

    // Fallback to convergent small-t series
    return eri::math::detail::small_t_approximation(m, t);
}

} // namespace eri math