#pragma once

#include <eri/basis/cgf.h>
#include <eri/chem/molecule.h>

#include <vector>
#include <string>
#include <stdexcept>

namespace eri::basis {

class BasisSet {
private:
    std::vector<CGF> cgfs_;

public:
    void load_from_bse_json(
        const std::string& filename,
        const eri::chem::Molecule& mol
    );

    const std::vector<CGF>& functions() const noexcept {
        return cgfs_;
    }

    const CGF& operator[](std::size_t i) const noexcept {
        return cgfs_[i];
    }

    const CGF& at(std::size_t i) const {
        if (i >= cgfs_.size())
            throw std::out_of_range("BasisSet::at: index out of range");
        return cgfs_[i];
    }

    auto begin() const noexcept {
        return cgfs_.begin();
    }

    auto end() const noexcept {
        return cgfs_.end();
    }

    std::size_t size() const noexcept {
        return cgfs_.size();
    }
};

} // namespace eri::basis
