#pragma once

#include <vector>
#include <array>
#include <cstddef>

namespace eri::cog {

/**
 * @brief One-dimensional Hellsing kernel for a fixed angular momentum quartet
 *        (l1, l2, l3, l4) along a single Cartesian axis.
 *
 * This structure contains *only* integer combinatorics and scalar prefactors.
 * It is completely independent of:
 *   - Gaussian exponents
 *   - Nuclear coordinates
 *   - Basis function centers
 *
 * As such, instances of this type are:
 *   - constructed once (during cache initialization),
 *   - immutable after construction,
 *   - safe to share between threads without synchronization.
 *
 * The kernel encodes terms of the form:
 *
 *   scalar[i] *
 *   a1^{p[0]} a2^{p[1]} a3^{p[2]} a4^{p[3]} *
 *   gamma1^{p[4]} gamma2^{p[5]} *
 *   (A-B)^{p[6]} (C-D)^{p[7]} *
 *   eta^{p[8]} (P-Q)^{p[9]}
 *
 * together with the associated Boys-function index
 *
 *   nu = mu[i] - u[i]
 */
struct HellsingCache1D
{
    /// Constant (floating-point) prefactors for each kernel term
    std::vector<double> scalar;

    /**
     * Integer exponent patterns for each kernel term.
     *
     * powers[i][k] corresponds to the exponent of the k-th base:
     *
     *   [0]  a1      (Gaussian exponent of primitive 1)
     *   [1]  a2      (Gaussian exponent of primitive 2)
     *   [2]  a3      (Gaussian exponent of primitive 3)
     *   [3]  a4      (Gaussian exponent of primitive 4)
     *   [4]  gamma1  = a1 + a2   (shifted by −(l1+l2))
     *   [5]  gamma2  = a3 + a4   (shifted by −(l3+l4))
     *   [6]  x1      = A − B     (Cartesian displacement, 1st pair)
     *   [7]  x2      = C − D     (Cartesian displacement, 2nd pair)
     *   [8]  eta     = gamma1*gamma2 / (gamma1+gamma2)
     *   [9]  PQ      = P − Q     (Gaussian product center separation)
     *
     * NOTE:
     *   Exponents may be *negative* (notably gamma1/gamma2),
     *   so the element type must be signed.
     */
    std::vector<std::array<int, 10>> powers;

    /// Total angular momentum order for each term (μ)
    std::vector<int> mu;

    /// Auxiliary summation index for each term (u)
    std::vector<int> u;
};


/**
 * @brief Table of all one-dimensional Hellsing kernels for angular momenta
 *        0 ≤ l1,l2,l3,l4 ≤ lmax.
 *
 * This table is constructed once (typically at the start of an ERI build),
 * and thereafter accessed read-only from multiple threads.
 *
 * The storage layout is a flat array indexed by the 4-tuple (l1,l2,l3,l4).
 */
class HellsingCacheTable1D
{
public:
    /**
     * @brief Construct the kernel table up to a maximum angular momentum.
     *
     * @param lmax Maximum angular momentum quantum number appearing
     *             in the basis set.
     */
    explicit HellsingCacheTable1D(int lmax);

    /**
     * @brief Retrieve the one-dimensional kernel for a given angular
     *        momentum quartet.
     *
     * @param l1 Angular momentum of primitive 1
     * @param l2 Angular momentum of primitive 2
     * @param l3 Angular momentum of primitive 3
     * @param l4 Angular momentum of primitive 4
     *
     * @return Reference to an immutable kernel object.
     */
    const HellsingCache1D& get(int l1, int l2, int l3, int l4) const;

    /// Maximum angular momentum supported by this table
    int lmax() const noexcept { return lmax_; }

private:
    /// Maximum angular momentum quantum number
    int lmax_;

    /// Flat storage of all kernels
    std::vector<HellsingCache1D> table_;

    /**
     * @brief Map a 4-tuple (l1,l2,l3,l4) into a flat array index.
     */
    static std::size_t index(int l1, int l2, int l3, int l4, int lmax);

    /**
     * @brief Build a one-dimensional kernel for fixed (l1,l2,l3,l4).
     *
     * This is a direct translation of the Python
     * `_calculate_coefficients()` routine, excluding any geometry-
     * or exponent-dependent quantities.
     */
    static HellsingCache1D build_kernel_1d(int l1, int l2, int l3, int l4);
};

} // namespace eri::cog
