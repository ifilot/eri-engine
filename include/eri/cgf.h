#pragma once

#include <vector>
#include <array>

namespace eri {

class CGF {
private:
    int lx_, ly_, lz_;
    std::array<double,3> center_;
    std::vector<double> exp_;
    std::vector<double> coef_;

public:
    CGF(int lx, int ly, int lz,
        std::array<double,3> ctr,
        std::vector<double> exp,
        std::vector<double> coef);

    int lx() const noexcept { return lx_; }
    int ly() const noexcept { return ly_; }
    int lz() const noexcept { return lz_; }

    const std::array<double,3>& ctr() const noexcept { return center_; }
    const std::vector<double>& exp() const noexcept { return exp_; }
    const std::vector<double>& coef() const noexcept { return coef_; }
};

}