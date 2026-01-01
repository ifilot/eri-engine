#pragma once

#include <vector>
#include <cstddef>
#include <cassert>
#include <algorithm>

#include <eri/two_electron/eri_index.h>

namespace eri::math {

class SymmetricERITensor {
private:
    std::size_t n_;
    std::size_t npairs_;
    std::vector<double> data_;

public:
    explicit SymmetricERITensor(std::size_t n)
        : n_(n),
          npairs_(n * (n + 1) / 2),
          data_(npairs_ * (npairs_ + 1) / 2)
    {}

    std::size_t dim() const noexcept {
        return n_;
    }

    std::size_t size() const noexcept {
        return data_.size();
    }

    double* data() noexcept {
        return data_.data();
    }

    const double* data() const noexcept {
        return data_.data();
    }

    double& operator()(std::size_t i,
                       std::size_t j,
                       std::size_t k,
                       std::size_t l) noexcept
    {
        assert(i < n_ && j < n_ && k < n_ && l < n_);
        const std::size_t idx = eri::two_electron::eri_index(i, j, k, l);
        return data_[idx];
    }

    double operator()(std::size_t i,
                      std::size_t j,
                      std::size_t k,
                      std::size_t l) const noexcept
    {
        assert(i < n_ && j < n_ && k < n_ && l < n_);
        const std::size_t idx = eri::two_electron::eri_index(i, j, k, l);
        return data_[idx];
    }
};

} // namespace eri::math
