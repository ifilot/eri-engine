#pragma once

#include <vector>
#include <array>
#include <stdexcept>

namespace eri::chem {

struct Atom {
    int Z;
    std::array<double,3> position;
};

class Molecule {
private:
    std::vector<Atom> atoms_;

public:
    void add_atom(int Z, const std::array<double,3>& pos) {
        atoms_.push_back({Z, pos});
    }

    const std::vector<Atom>& atoms() const noexcept {
        return atoms_;
    }

    const auto& operator[](std::size_t i) const noexcept {
        return atoms_[i];
    }

    const auto& at(std::size_t i) const {
        if (i >= atoms_.size())
            throw std::out_of_range("BasisSet::at: index out of range");
        return atoms_[i];
    }

    auto begin() const noexcept {
        return atoms_.begin();
    }

    auto end() const noexcept {
        return atoms_.end();
    }

    auto cbegin() const noexcept {
        return atoms_.cbegin();
    }

    auto cend() const noexcept {
        return atoms_.cend();
    }

    auto rbegin() const noexcept {
        return atoms_.crbegin();
    }

    auto rend() const noexcept {
        return atoms_.crend();
    }

    auto crbegin() const noexcept {
        return atoms_.crbegin();
    }

    auto crend() const noexcept {
        return atoms_.crend();
    }

    std::size_t size() const noexcept {
        return atoms_.size();
    }
};

} // namespace eri::chem
