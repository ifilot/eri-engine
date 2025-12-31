#include <eri/one_electron/nuclear.h>

#include "detail/nuclear_primitive.h"

namespace eri::one_electron {

double nuclear(const eri::basis::CGF& a,
               const eri::basis::CGF& b,
               const std::array<double, 3>& C,
               eri::enums::NuclearMethod method) {
    const auto& A = a.ctr();
    const auto& B = b.ctr();

    const auto& ea = a.exp();
    const auto& eb = b.exp();

    const auto& ca = a.coef();
    const auto& cb = b.coef();

    const auto& na = a.norm();
    const auto& nb = b.norm();

    double sum = 0.0;

    for (std::size_t i = 0; i < ea.size(); ++i) {
        for (std::size_t j = 0; j < eb.size(); ++j) {
            const double prim = eri::detail::nuclear_primitive(
                a.lx(), a.ly(), a.lz(), A, ea[i],
                b.lx(), b.ly(), b.lz(), B, eb[j],
                C,
                method
            );

            sum += (ca[i] * na[i]) * (cb[j] * nb[j]) * prim;
        }
    }

    return sum;
}

} // namespace eri::one_electron
