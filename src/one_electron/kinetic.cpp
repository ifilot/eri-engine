#include <eri/one_electron/overlap.h>
#include <eri/one_electron/kinetic.h>

#include "detail/kinetic_primitive.h"

#include <array>

namespace eri::one_electron {

double kinetic(const eri::basis::CGF& a,
               const eri::basis::CGF& b,
               eri::enums::KineticMethod method) {
    double T = 0.0;

    const auto& A = a.ctr();
    const auto& B = b.ctr();

    for (std::size_t p = 0; p < a.exp().size(); ++p) {
        const double ap = a.exp()[p];
        const double cp = a.coef()[p] * a.norm()[p];

        for (std::size_t q = 0; q < b.exp().size(); ++q) {
            const double bq = b.exp()[q];
            const double dq = b.coef()[q] * b.norm()[q];

            const double t_pq = eri::detail::kinetic_primitive(
                a.lx(), a.ly(), a.lz(), A, ap,
                b.lx(), b.ly(), b.lz(), B, bq,
                method
            );

            T += cp * dq * t_pq;
        }
    }

    return T;
}

} // namespace eri::one_electron
