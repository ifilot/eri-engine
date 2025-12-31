#include <eri/basis/cgf.h>
#include <eri/one_electron/overlap.h>

#include "detail/normalization.h"

/**
 * @brief  Construct CGF
 * @note   assumes all primitives have the same (lx, ly, lz and ctr) values
 * @param  lx: angular momentum in x-direction
 * @param  ly: angular momentum in y-direction
 * @param  lz: angular momentum in z-direction
 * @param  ctr: CGF center
 * @param  exp: alpha values for the primitives
 * @param  coef: expansion coefficients of the primitives
 * @retval 
 */
eri::basis::CGF::CGF(int lx, int ly, int lz,
                     std::array<double,3> ctr,
                     std::vector<double> exp,
                     std::vector<double> coef)
  : lx_(lx), 
    ly_(ly), 
    lz_(lz),
    center_(ctr),
    exp_(std::move(exp)),
    coef_(std::move(coef)) 
{
    if (lx < 0 || ly < 0 || lz < 0) {
        throw std::invalid_argument("Angular momentum must be non-negative");
    }

    if (exp_.empty()) {
        throw std::invalid_argument("CGF must contain at least one primitive");
    }

    if (exp_.size() != coef_.size()) {
        throw std::invalid_argument("Exponent/coefficient size mismatch");
    }

    this->normalize();
}

/**
 * @brief  Normalize the primitives of the CGF 
 * @note   Called upon class construction
 * @retval None
 */
void eri::basis::CGF::normalize() {
    const std::size_t n = this->exp_.size();

    if (n == 0) {
        throw std::runtime_error("CGF has no primitives");
    }

    // allocate normalization constants
    this->norm_.resize(n);

    // compute primitive normalization constants
    for (std::size_t p = 0; p < n; ++p) {
        this->norm_[p] = eri::detail::primitive_norm(
            this->lx_,
            this->ly_,
            this->lz_,
            this->exp_[p]
        );
    }
}

/**
 * @brief  Calculate the self-overlap
 * @retval Self-overlap
 */
double eri::basis::CGF::overlap() const {
    return eri::one_electron::overlap(*this, *this);
}