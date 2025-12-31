#include <eri/one_electron/overlap.h>

#include <cmath>
#include <array>

namespace {

inline double factorial_small(int n) noexcept {
    // Enough for typical l <= ~10 use; increase if needed.
    static constexpr double f[] = {
        1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0,
        3628800.0, 39916800.0, 479001600.0, 6227020800.0, 87178291200.0,
        1307674368000.0, 20922789888000.0, 355687428096000.0,
        6402373705728000.0, 121645100408832000.0
    };
    if (n < 0) return 1.0;
    if (n < static_cast<int>(sizeof(f)/sizeof(f[0]))) return f[n];
    // fallback (rare)
    return std::tgamma(static_cast<double>(n) + 1.0);
}

inline double overlap_1d_hellsing(int l1, int l2, double a1, double a2, double x, double gamma) {
    double pre = factorial_small(l1) * factorial_small(l2) / std::pow(gamma, l1 + l2);
    if (l1 & 1) pre *= -1.0;

    double sm = 0.0;

    for (int i1 = 0; i1 <= l1/2; ++i1) {
        for (int i2 = 0; i2 <= l2/2; ++i2) {
            const int j = l1 + l2 - 2*i1 - 2*i2;

            for (int m = 0; m <= j/2; ++m) {
                const double t1 =
                    std::pow(-1.0, m) *
                    factorial_small(j) *
                    std::pow(a1, l2 - i1 - 2*i2 - m) *
                    std::pow(a2, l1 - 2*i1 - i2 - m) /
                    std::pow(4.0, i1 + i2 + m) /
                    factorial_small(i1) / factorial_small(i2) / factorial_small(m);

                const double t2 =
                    std::pow(gamma, 2*(i1+i2) + m) *
                    std::pow(x, j - 2*m) /
                    factorial_small(l1 - 2*i1) /
                    factorial_small(l2 - 2*i2) /
                    factorial_small(j - 2*m);

                sm += t1 * t2;
            }
        }
    }

    return pre * sm;
}

inline double overlap_primitive_hellsing(
    int lx1, int ly1, int lz1, const std::array<double,3>& A, double a,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double b) {
    const double dx = A[0] - B[0];
    const double dy = A[1] - B[1];
    const double dz = A[2] - B[2];
    const double AB2 = dx*dx + dy*dy + dz*dz;

    const double gamma = a + b;

    const double pre = std::pow(M_PI / gamma, 1.5) * std::exp(-a * b * AB2 / gamma);

    // Hellsing uses x = A - B (componentwise)
    const double wx = overlap_1d_hellsing(lx1, lx2, a, b, A[0] - B[0], gamma);
    const double wy = overlap_1d_hellsing(ly1, ly2, a, b, A[1] - B[1], gamma);
    const double wz = overlap_1d_hellsing(lz1, lz2, a, b, A[2] - B[2], gamma);

    return pre * wx * wy * wz;
}

} // namespace

namespace eri::one_electron {

double overlap_hellsing(const eri::basis::CGF& a, const eri::basis::CGF& b) {
    double S = 0.0;

    const auto& A = a.ctr();
    const auto& B = b.ctr();

    for (std::size_t p = 0; p < a.exp().size(); ++p) {
        const double ap = a.exp()[p];
        const double cp = a.coef()[p] * a.norm()[p];

        for (std::size_t q = 0; q < b.exp().size(); ++q) {
            const double bq = b.exp()[q];
            const double dq = b.coef()[q] * b.norm()[q];

            const double s_pq = overlap_primitive_hellsing(
                a.lx(), a.ly(), a.lz(), A, ap,
                b.lx(), b.ly(), b.lz(), B, bq
            );

            S += cp * dq * s_pq;
        }
    }

    return S;
}

} // namespace eri::one_electron
