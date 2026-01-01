#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <eri/chem/molecule.h>
#include <eri/basis/basisset.h>
#include <eri/two_electron/eri.h>
#include <eri/utils/executable_dir.h>

static eri::chem::Molecule make_h2()
{
    eri::chem::Molecule mol;
    mol.add_atom(1, { 0.00000, 0.00000, 0.00000 });
    mol.add_atom(1, { 0.00000, 0.00000, 1.40000 });
    return mol;
}

static eri::basis::BasisSet make_basis_h2() {
    eri::chem::Molecule mol = make_h2();
    eri::basis::BasisSet basis;
    basis.load_from_bse_json(eri::utils::executable_dir() + "/../data/basis/sto-3g.json", mol);
    return basis;
}

TEST_CASE("ERI Huzinaga (STO-3G H2)") {
    auto basis = make_basis_h2();

    double T1111 = eri::two_electron::eri(basis[0], basis[0], basis[0], basis[0]);
    double T1122 = eri::two_electron::eri(basis[0], basis[0], basis[1], basis[1]);
    double T1112 = eri::two_electron::eri(basis[0], basis[0], basis[0], basis[1]);
    double T2121 = eri::two_electron::eri(basis[1], basis[0], basis[1], basis[0]);
    double T1222 = eri::two_electron::eri(basis[0], basis[1], basis[1], basis[1]);
    double T2211 = eri::two_electron::eri(basis[1], basis[1], basis[0], basis[0]);

    const double epsilon = 1e-6;

    CHECK(T1111 == doctest::Approx(0.774605).epsilon(epsilon));
    CHECK(T1122 == doctest::Approx(0.569675).epsilon(epsilon));
    CHECK(T1112 == doctest::Approx(0.444107).epsilon(epsilon));
    CHECK(T2121 == doctest::Approx(0.297028).epsilon(epsilon));

    CHECK(T1222 == doctest::Approx(T1112).epsilon(epsilon));
    CHECK(T1122 == doctest::Approx(T2211).epsilon(epsilon));
}