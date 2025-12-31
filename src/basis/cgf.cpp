#include "../../include/eri/basis/cgf.h"
#include "../../include/eri/one_electron/overlap.h"

#include "detail/normalization.h"

eri::basis::CGF::CGF(int lx, int ly, int lz,
              std::array<double,3> ctr,
              std::vector<double> exp,
              std::vector<double> coef)
  : lx_(lx), ly_(ly), lz_(lz),
    center_(ctr),
    exp_(std::move(exp)),
    coef_(std::move(coef)) {
    
    // check validity
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


double eri::basis::CGF::overlap() const {
    return eri::one_electron::overlap(*this, *this);
}