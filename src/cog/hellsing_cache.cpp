#include <eri/cog/hellsing_cache.h>

#include "math/factorial.h"
#include "math/sign_pow.h"
#include "math/ipow.h"

namespace eri::cog {

std::size_t HellsingCacheTable::index(int l1,int l2,int l3,int l4,int lmax) {
    const std::size_t n = lmax + 1;
    return ((l1*n + l2)*n + l3)*n + l4;
}

HellsingCacheTable::HellsingCacheTable(int lmax) : lmax_(lmax) {
    const std::size_t n = (lmax+1)*(lmax+1)*(lmax+1)*(lmax+1);
    table_.resize(n);

    for(int l1=0;l1<=lmax;++l1)
        for(int l2=0;l2<=lmax;++l2)
            for(int l3=0;l3<=lmax;++l3)
                for(int l4=0;l4<=lmax;++l4)
                    table_[index(l1,l2,l3,l4,lmax)] = build_kernel_1d(l1,l2,l3,l4);
}

const HellsingCache1D& HellsingCacheTable::get(int l1,int l2,int l3,int l4) const {
    return table_[index(l1,l2,l3,l4,lmax_)];
}

/**
 * Direct translation of Python _calculate_coefficients(),
 * but WITHOUT geometry or exponent dependence.
 */
HellsingCache1D HellsingCacheTable::build_kernel_1d(int l1,int l2,int l3,int l4) {
    HellsingCache1D K;

    const double pre1 =
        eri::math::sign_pow(l1+l2) *
        eri::math::factorial(l1) *
        eri::math::factorial(l2);

    const double pre2 =
        eri::math::factorial(l3) *
        eri::math::factorial(l4);

    for(int i1=0;i1<=l1/2;++i1)
        for(int i2=0;i2<=l2/2;++i2)
            for(int o1=0;o1<=l1-2*i1;++o1)
                for(int o2=0;o2<=l2-2*i2;++o2)
                    for(int r1=0;r1<=(o1+o2)/2;++r1)
                    {
                        const double t11 =
                            eri::math::sign_pow(o2+r1) *
                            eri::math::factorial(o1+o2) /
                            (eri::math::ipow(4.0,i1+i2+r1) *
                             eri::math::factorial(i1) *
                             eri::math::factorial(i2) *
                             eri::math::factorial(o1) *
                             eri::math::factorial(o2) *
                             eri::math::factorial(r1));

                        const double d12 =
                            eri::math::factorial(l1-2*i1-o1) *
                            eri::math::factorial(l2-2*i2-o2) *
                            eri::math::factorial(o1+o2-2*r1);

                        for(int i3=0;i3<=l3/2;++i3)
                            for(int i4=0;i4<=l4/2;++i4)
                                for(int o3=0;o3<=l3-2*i3;++o3)
                                    for(int o4=0;o4<=l4-2*i4;++o4)
                                        for(int r2=0;r2<=(o3+o4)/2;++r2)
                                        {
                                            const double t21 =
                                                eri::math::sign_pow(o3+r2) *
                                                eri::math::factorial(o3+o4) /
                                                (eri::math::ipow(4.0,i3+i4+r2) *
                                                 eri::math::factorial(i3) *
                                                 eri::math::factorial(i4) *
                                                 eri::math::factorial(o3) *
                                                 eri::math::factorial(o4) *
                                                 eri::math::factorial(r2));

                                            const double d22 =
                                                eri::math::factorial(l3-2*i3-o3) *
                                                eri::math::factorial(l4-2*i4-o4) *
                                                eri::math::factorial(o3+o4-2*r2);

                                            const int mu =
                                                l1+l2+l3+l4
                                                -2*(i1+i2+i3+i4)
                                                -(o1+o2+o3+o4);

                                            for(int u=0;u<=mu/2;++u)
                                            {
                                                const double t3 =
                                                    eri::math::sign_pow(u) *
                                                    eri::math::factorial(mu) /
                                                    (eri::math::ipow(4.0,u) *
                                                     eri::math::factorial(u) *
                                                     eri::math::factorial(mu-2*u));

                                                K.scalar.push_back(
                                                    pre1 * pre2 * t11 / d12 * t21 / d22 * t3
                                                );

                                                K.powers.push_back({
                                                    o2-i1-r1,        // a1
                                                    o1-i2-r1,        // a2
                                                    o4-i3-r2,        // a3
                                                    o3-i4-r2,        // a4
                                                    2*(i1+i2)+r1,    // g1
                                                    2*(i3+i4)+r2,    // g2
                                                    o1+o2-2*r1,      // x1
                                                    o3+o4-2*r2,      // x2
                                                    mu-u,            // eta
                                                    mu-2*u           // pq
                                                });

                                                K.mu.push_back(mu);
                                                K.u.push_back(u);
                                            }
                                        }
                    }

    return K;
}

} // namespace eri::cog
