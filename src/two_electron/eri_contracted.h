#pragma once

#include <eri/basis/cgf.h>

namespace eri::two_electron::primitives {

template <typename PrimitiveEvaluator>
double eri_contracted(
    const eri::basis::CGF& a,
    const eri::basis::CGF& b,
    const eri::basis::CGF& c,
    const eri::basis::CGF& d,
    const PrimitiveEvaluator& prim_eval
) {
    double T = 0.0;

    const auto& A = a.ctr();
    const auto& B = b.ctr();
    const auto& C = c.ctr();
    const auto& D = d.ctr();

    for (std::size_t p = 0; p < a.exp().size(); ++p) {
        const double ap = a.exp()[p];
        const double ac = a.coef()[p] * a.norm()[p];

        for (std::size_t q = 0; q < b.exp().size(); ++q) {
            const double bp = b.exp()[q];
            const double bc = b.coef()[q] * b.norm()[q];

            for (std::size_t r = 0; r < c.exp().size(); ++r) {
                const double cp = c.exp()[r];
                const double cc = c.coef()[r] * c.norm()[r];

                for (std::size_t s = 0; s < d.exp().size(); ++s) {
                    const double dp = d.exp()[s];
                    const double dc = d.coef()[s] * d.norm()[s];

                    const double t = prim_eval(
                        a.lx(), a.ly(), a.lz(), A, ap,
                        b.lx(), b.ly(), b.lz(), B, bp,
                        c.lx(), c.ly(), c.lz(), C, cp,
                        d.lx(), d.ly(), d.lz(), D, dp
                    );

                    T += ac * bc * cc * dc * t;
                }
            }
        }
    }

    return T;
}

}