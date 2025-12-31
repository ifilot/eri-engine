#pragma once

#include <vector>
#include <array>

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
};

} // namespace eri::chem
