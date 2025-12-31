#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <eri/basis/basisset.h>
#include <eri/chem/molecule.h>
#include <eri/one_electron/overlap.h>

#include <unistd.h>
#include <limits.h>

#include <string>
#include <stdexcept>

// Resolve directory of the running executable (Linux)
static std::string executable_dir()
{
    char buffer[PATH_MAX];
    ssize_t len = ::readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);

    if (len == -1)
        throw std::runtime_error("Failed to resolve /proc/self/exe");

    buffer[len] = '\0';

    std::string path(buffer);
    return path.substr(0, path.find_last_of('/'));
}

TEST_CASE("BasisSet can be populated from BSE JSON") {

    // --- Build H2 molecule ---
    eri::chem::Molecule mol;
    mol.add_atom(1, {0.0, 0.0, 0.0});
    mol.add_atom(1, {0.0, 0.0, 1.4});

    // --- Resolve basis file relative to executable ---
    const std::string basis_file =
        executable_dir() + "/data/sto-3g.json";

    // --- Load basis ---
    eri::basis::BasisSet basis;
    CHECK_NOTHROW(basis.load_from_bse_json(basis_file, mol));

    // --- Basic sanity checks ---
    CHECK(basis.size() > 0);

    // STO-3G hydrogen has exactly 1 Cartesian function per atom
    CHECK(basis.size() == 2);

    // Check angular momentum of both CGFs
    for (const auto& cgf : basis) {
        CHECK(cgf.lx() == 0);
        CHECK(cgf.ly() == 0);
        CHECK(cgf.lz() == 0);
        CHECK(cgf.exp().size() == 3);
        CHECK(cgf.coef().size() == 3);
    }

    // check overlap
    auto a = basis[0];
    auto b = basis[1];
    const double S11 = eri::one_electron::overlap(a, a, eri::enums::OverlapMethod::Hellsing);
    const double S12 = eri::one_electron::overlap(a, b, eri::enums::OverlapMethod::Hellsing);
    const double S22 = eri::one_electron::overlap(b, b, eri::enums::OverlapMethod::Hellsing);

    CHECK(S11 == doctest::Approx(1.0));
    CHECK(S12 == doctest::Approx(0.65931845));
    CHECK(S22 == doctest::Approx(1.0));
}