#pragma once

#include <eri/math/dense_square_matrix.h>

namespace eri::math {

template <typename Operator>
DenseSquareMatrix build_symmetric_matrix(const basis::BasisSet& basis)
{
    const std::size_t n = basis.size();
    DenseSquareMatrix M(n);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            const double v = Operator::eval(basis[i], basis[j]);
            M(i, j) = v;
            M(j, i) = v;
        }
    }
    return M;
}

template <typename Operator>
DenseSquareMatrix build_symmetric_matrix(const basis::BasisSet& basis, const std::array<double, 3>& nuc, double chg)
{
    const std::size_t n = basis.size();
    DenseSquareMatrix M(n);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            const double v = chg * Operator::eval(basis[i], basis[j], nuc);
            M(i, j) = v;
            M(j, i) = v;
        }
    }
    return M;
}

} // namespace eri::math