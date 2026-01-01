#include "math/fgamma_acc.h"
#include "math/fgamma_nr.h"

#include <vector>

namespace eri::math {

double Fgamma_acc(int n, double T) {
    if (n < 0) return std::numeric_limits<double>::quiet_NaN();
    double out = 0.0;
    // For scalar, do a tiny block up to n
    std::vector<double> F(n + 1);
    detail::Fgamma_block_exact(n, T, F.data());
    out = F[n];
    return out;
}

double Fgamma_nr(int a, double x) {
    x = std::fabs(x);
    if (x < 1e-7) x = 1e-7;

    const double val = detail::gamm_inc(a + 0.5, x);
    if (std::isnan(val)) return val;

    return 0.5 * std::pow(x, -a - 0.5) * val;
}

} // namespace eri::math
