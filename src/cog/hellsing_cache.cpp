#include <eri/cog/hellsing_cache.h>

#include "math/factorial.h"
#include "math/sign_pow.h"
#include "math/ipow.h"

namespace eri::cog {

/**
 * @brief Compute flat-table index for a given angular momentum quartet.
 *
 * The kernel table is stored as a flat array over all
 * 0 ≤ l1,l2,l3,l4 ≤ lmax.
 */
std::size_t
HellsingCacheTable1D::index(int l1, int l2, int l3, int l4, int lmax)
{
    const std::size_t n = static_cast<std::size_t>(lmax) + 1;
    return ((l1 * n + l2) * n + l3) * n + l4;
}


/**
 * @brief Construct the full Hellsing kernel table up to angular momentum lmax.
 *
 * All one-dimensional kernels for every (l1,l2,l3,l4) combination
 * are generated eagerly.  After construction, the table is immutable
 * and safe for concurrent read-only access.
 */
HellsingCacheTable1D::HellsingCacheTable1D(int lmax)
: lmax_(lmax) {
    const std::size_t n =
        static_cast<std::size_t>(lmax + 1) *
        static_cast<std::size_t>(lmax + 1) *
        static_cast<std::size_t>(lmax + 1) *
        static_cast<std::size_t>(lmax + 1);

    table_.resize(n);

    for (int l1 = 0; l1 <= lmax_; ++l1)
        for (int l2 = 0; l2 <= lmax_; ++l2)
            for (int l3 = 0; l3 <= lmax_; ++l3)
                for (int l4 = 0; l4 <= lmax_; ++l4)
                    table_[index(l1, l2, l3, l4, lmax_)] =
                        build_kernel_1d(l1, l2, l3, l4);
}


/**
 * @brief Retrieve the cached one-dimensional kernel for (l1,l2,l3,l4).
 */
const HellsingCache1D&
HellsingCacheTable1D::get(int l1, int l2, int l3, int l4) const
{
    return table_[index(l1, l2, l3, l4, lmax_)];
}


/**
 * @brief Build the one-dimensional Hellsing kernel for fixed (l1,l2,l3,l4).
 *
 * This routine is a direct, term-by-term translation of the Python
 * `_calculate_coefficients()` method used in the reference implementation.
 *
 * IMPORTANT:
 *   - Only integer combinatorics and scalar prefactors are generated here.
 *   - No geometry-, exponent-, or Boys-function dependence appears here.
 *   - The resulting kernel may contain *negative* integer exponents
 *     (notably for gamma1 and gamma2).
 */
HellsingCache1D
HellsingCacheTable1D::build_kernel_1d(int l1, int l2, int l3, int l4)
{
    HellsingCache1D K;

    // Overall angular-momentum prefactors
    const double pre1 =
        eri::math::sign_pow(l1 + l2) *
        eri::math::factorial(l1) *
        eri::math::factorial(l2);

    const double pre2 =
        eri::math::factorial(l3) *
        eri::math::factorial(l4);

    // Main combinatorial expansion
    for (int i1 = 0; i1 <= l1 / 2; ++i1)
        for (int i2 = 0; i2 <= l2 / 2; ++i2)
            for (int o1 = 0; o1 <= l1 - 2 * i1; ++o1)
                for (int o2 = 0; o2 <= l2 - 2 * i2; ++o2)
                    for (int r1 = 0; r1 <= (o1 + o2) / 2; ++r1)
                    {
                        // First-pair combinatorial factor (numerator)
                        const double t11 =
                            eri::math::sign_pow(o2 + r1) *
                            eri::math::factorial(o1 + o2) /
                            (eri::math::inpow(4.0, i1 + i2 + r1) *
                             eri::math::factorial(i1) *
                             eri::math::factorial(i2) *
                             eri::math::factorial(o1) *
                             eri::math::factorial(o2) *
                             eri::math::factorial(r1));

                        // First-pair combinatorial factor (denominator)
                        const double d12 =
                            eri::math::factorial(l1 - 2 * i1 - o1) *
                            eri::math::factorial(l2 - 2 * i2 - o2) *
                            eri::math::factorial(o1 + o2 - 2 * r1);

                        for (int i3 = 0; i3 <= l3 / 2; ++i3)
                            for (int i4 = 0; i4 <= l4 / 2; ++i4)
                                for (int o3 = 0; o3 <= l3 - 2 * i3; ++o3)
                                    for (int o4 = 0; o4 <= l4 - 2 * i4; ++o4)
                                        for (int r2 = 0; r2 <= (o3 + o4) / 2; ++r2)
                                        {
                                            // Second-pair combinatorial factor (numerator)
                                            const double t21 =
                                                eri::math::sign_pow(o3 + r2) *
                                                eri::math::factorial(o3 + o4) /
                                                (eri::math::inpow(4.0, i3 + i4 + r2) *
                                                 eri::math::factorial(i3) *
                                                 eri::math::factorial(i4) *
                                                 eri::math::factorial(o3) *
                                                 eri::math::factorial(o4) *
                                                 eri::math::factorial(r2));

                                            // Second-pair combinatorial factor (denominator)
                                            const double d22 =
                                                eri::math::factorial(l3 - 2 * i3 - o3) *
                                                eri::math::factorial(l4 - 2 * i4 - o4) *
                                                eri::math::factorial(o3 + o4 - 2 * r2);

                                            // Total Boys-function order μ
                                            const int mu =
                                                l1 + l2 + l3 + l4
                                                - 2 * (i1 + i2 + i3 + i4)
                                                - (o1 + o2 + o3 + o4);

                                            for (int u = 0; u <= mu / 2; ++u)
                                            {
                                                // Final combinatorial factor
                                                const double t3 =
                                                    eri::math::sign_pow(u) *
                                                    eri::math::factorial(mu) /
                                                    (eri::math::inpow(4.0, u) *
                                                     eri::math::factorial(u) *
                                                     eri::math::factorial(mu - 2 * u));

                                                // Scalar prefactor for this term
                                                K.scalar.push_back(
                                                    pre1 * pre2 *
                                                    t11 / d12 *
                                                    t21 / d22 *
                                                    t3
                                                );

                                                // Integer exponent pattern (see header for meaning)
                                                K.powers.push_back({
                                                    o2 - i1 - r1,                 // a1
                                                    o1 - i2 - r1,                 // a2
                                                    o4 - i3 - r2,                 // a3
                                                    o3 - i4 - r2,                 // a4
                                                    2 * (i1 + i2) + r1 - (l1 + l2), // gamma1
                                                    2 * (i3 + i4) + r2 - (l3 + l4), // gamma2
                                                    o1 + o2 - 2 * r1,             // x1
                                                    o3 + o4 - 2 * r2,             // x2
                                                    mu - u,                       // eta
                                                    mu - 2 * u                    // PQ
                                                });

                                                K.mu.push_back(mu);
                                                K.u.push_back(u);
                                            }
                                        }
                    }

    return K;
}

} // namespace eri::cog
