#pragma once

#include <vector>
#include <cstddef>
#include <cassert>

namespace eri::math {

class DenseSquareMatrix {
private:
    std::size_t n_;
    std::vector<double> data_; // row-major

public:
    explicit DenseSquareMatrix(std::size_t n)
        : n_(n), data_(n * n) {}

    std::size_t size() const noexcept {
        return n_;
    }

    double* data() noexcept {
        return data_.data();
    }

    const double* data() const noexcept {
        return data_.data();
    }

    double& operator()(std::size_t i, std::size_t j) noexcept {
        assert(i < n_ && j < n_);
        return data_[i * n_ + j];
    }

    double operator()(std::size_t i, std::size_t j) const noexcept {
        assert(i < n_ && j < n_);
        return data_[i * n_ + j];
    }
};

} // namespace eri::math
