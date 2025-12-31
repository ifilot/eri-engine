#include <eri/one_electron/overlap.h>

#include "detail/overlap_hellsing.h"

#include <cmath>
#include <array>

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

            const double s_pq = eri::detail::overlap_primitive_hellsing(
                a.lx(), a.ly(), a.lz(), A, ap,
                b.lx(), b.ly(), b.lz(), B, bq
            );

            S += cp * dq * s_pq;
        }
    }

    return S;
}

} // namespace eri::one_electron
