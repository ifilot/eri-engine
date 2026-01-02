#pragma once

#include <vector>
#include <array>

namespace eri::cog {

/**
 * One-dimensional Hellsing kernel for fixed (l1,l2,l3,l4).
 * Contains only integer structure and scalar constants.
 *
 * This is SAFE to share across threads (read-only).
 */
struct HellsingCache1D {
    std::vector<double> scalar;                 // constant prefactors
    std::vector<std::array<int,10>> powers;     // exponent patterns
    std::vector<int> mu;                        // mu indices
    std::vector<int> u;                         // u indices
};

/**
 * Kernel table for all (l1,l2,l3,l4) up to lmax.
 * Built once, then read-only.
 */
class HellsingCacheTable {
public:
    explicit HellsingCacheTable(int lmax);

    const HellsingCache1D& get(int l1, int l2, int l3, int l4) const;

    int lmax() const noexcept { return lmax_; }

private:
    int lmax_; std::vector<HellsingCache1D> table_;

    static std::size_t index(int l1, int l2, int l3, int l4, int lmax);

    static HellsingCache1D build_kernel_1d(int l1, int l2, int l3, int l4);
};

} // namespace eri::detail
