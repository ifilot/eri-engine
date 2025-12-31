#include <eri/one_electron/overlap.h>

#include "math/gaussian_product.h"
#include "math/binomial_prefactor.h"
#include "math/double_factorial.h"

#include <cmath>

namespace {

inline double overlap_1d_huzinaga(int l1, int l2, double x1, double x2, double gamma) {
    double sm = 0.0;
    const int imax = (l1 + l2) / 2;

    for (int i = 0; i <= imax; ++i) {
        const int k = 2 * i;
        const double df = (i == 0) ? 1.0 : eri::math::double_factorial(2*i - 1);

        sm += eri::math::binomial_prefactor(k, l1, l2, x1, x2) *
              df / std::pow(2.0 * gamma, i);
    }
    return sm;
}

inline double overlap_primitive_huzinaga(
    int lx1, int ly1, int lz1, const std::array<double,3>& A, double a,
    int lx2, int ly2, int lz2, const std::array<double,3>& B, double b)
{
    const double dx = A[0] - B[0];
    const double dy = A[1] - B[1];
    const double dz = A[2] - B[2];
    const double AB2 = dx*dx + dy*dy + dz*dz;

    const double gamma = a + b;
    const auto P = eri::math::gaussian_product_center(a, A, b, B);

    const double pre = std::pow(M_PI / gamma, 1.5) * std::exp(-a * b * AB2 / gamma);

    const double wx = overlap_1d_huzinaga(lx1, lx2, P[0] - A[0], P[0] - B[0], gamma);
    const double wy = overlap_1d_huzinaga(ly1, ly2, P[1] - A[1], P[1] - B[1], gamma);
    const double wz = overlap_1d_huzinaga(lz1, lz2, P[2] - A[2], P[2] - B[2], gamma);

    return pre * wx * wy * wz;
}

} // namespace

namespace eri::one_electron {

double overlap_huzinaga(const eri::basis::CGF& a, const eri::basis::CGF& b) {
    // Contracted overlap: sum_{p,q} (c_p N_p)(d_q M_q) * S_primitive(p,q)
    double S = 0.0;

    const auto& A = a.ctr();
    const auto& B = b.ctr();

    for (std::size_t p = 0; p < a.exp().size(); ++p) {
        const double ap = a.exp()[p];
        const double cp = a.coef()[p] * a.norm()[p];

        for (std::size_t q = 0; q < b.exp().size(); ++q) {
            const double bq = b.exp()[q];
            const double dq = b.coef()[q] * b.norm()[q];

            const double s_pq = overlap_primitive_huzinaga(
                a.lx(), a.ly(), a.lz(), A, ap,
                b.lx(), b.ly(), b.lz(), B, bq
            );

            S += cp * dq * s_pq;
        }
    }

    return S;
}

} // namespace eri::one_electron