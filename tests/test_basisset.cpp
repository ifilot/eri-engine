#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <eri/utils/executable_dir.h>
#include <eri/basis/basisset.h>
#include <eri/chem/molecule.h>
#include <eri/one_electron/overlap.h>

static eri::chem::Molecule make_h2o()
{
    eri::chem::Molecule mol;
    mol.add_atom(8, { 0.000000, -0.143207, 0.000000 });
    mol.add_atom(1, { 1.637236,  1.136548, 0.000000 });
    mol.add_atom(1, {-1.637236,  1.136548, 0.000000 });
    return mol;
}

static eri::basis::BasisSet make_basis_h2o() {
    eri::chem::Molecule mol = make_h2o();
    eri::basis::BasisSet basis;
    basis.load_from_bse_json(eri::utils::executable_dir() + "/../data/basis/sto-3g.json", mol);
    return basis;
}

TEST_CASE("BasisSet can be populated from BSE JSON") {

    // --- Build H2 molecule ---
    eri::chem::Molecule mol;
    mol.add_atom(1, {0.0, 0.0, 0.0});
    mol.add_atom(1, {0.0, 0.0, 1.4});

    // --- Resolve basis file relative to executable ---
    const std::string basis_file =
        eri::utils::executable_dir() + "/../data/basis/sto-3g.json";

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

    const std::vector<double> exps = {
        3.42525091,
        0.62391373,
        0.16885540
    };

    const std::vector<double> coefs = {
        0.15432897,
        0.53532814,
        0.44463454
    };

    for(std::size_t i=0; i<exps.size(); ++i) {
        CHECK(basis[0].exp()[i] == doctest::Approx(exps[i]));
        CHECK(basis[0].coef()[i] == doctest::Approx(coefs[i]));
        CHECK(basis[1].exp()[i] == doctest::Approx(exps[i]));
        CHECK(basis[1].coef()[i] == doctest::Approx(coefs[i]));
    }
}

TEST_CASE("Build basis for H2O (ordering and structure)")
{
    const auto basis = make_basis_h2o();

    // --- size ---
    REQUIRE(basis.size() == 7);

    // --- Oxygen center ---
    const std::array<double,3> O = { 0.000000, -0.143207, 0.000000 };

    CHECK(basis[0].lx() == 0);
    CHECK(basis[0].ly() == 0);
    CHECK(basis[0].lz() == 0);
    CHECK(basis[0].ctr() == O);  // O 1s

    CHECK(basis[1].lx() == 0);
    CHECK(basis[1].ly() == 0);
    CHECK(basis[1].lz() == 0);
    CHECK(basis[1].ctr() == O);  // O 2s

    CHECK(basis[2].lx() == 1);
    CHECK(basis[2].ly() == 0);
    CHECK(basis[2].lz() == 0);
    CHECK(basis[2].ctr() == O);  // O 2px

    CHECK(basis[3].lx() == 0);
    CHECK(basis[3].ly() == 1);
    CHECK(basis[3].lz() == 0);
    CHECK(basis[3].ctr() == O);  // O 2py

    CHECK(basis[4].lx() == 0);
    CHECK(basis[4].ly() == 0);
    CHECK(basis[4].lz() == 1);
    CHECK(basis[4].ctr() == O);  // O 2pz

    // --- Hydrogen centers ---
    const std::array<double,3> H1 = { 1.637236,  1.136548, 0.000000 };
    const std::array<double,3> H2 = {-1.637236,  1.136548, 0.000000 };

    CHECK(basis[5].lx() == 0);
    CHECK(basis[5].ly() == 0);
    CHECK(basis[5].lz() == 0);
    CHECK(basis[5].ctr() == H1); // H1 1s

    CHECK(basis[6].lx() == 0);
    CHECK(basis[6].ly() == 0);
    CHECK(basis[6].lz() == 0);
    CHECK(basis[6].ctr() == H2); // H2 1s

    CHECK(basis[0].exp()[0] == doctest::Approx(130.70932));
    CHECK(basis[0].exp()[1] == doctest::Approx(23.808861));
    CHECK(basis[0].exp()[2] == doctest::Approx(6.443608));

    CHECK(basis[0].coef()[0] == doctest::Approx(0.154329));
    CHECK(basis[0].coef()[1] == doctest::Approx(0.535328));
    CHECK(basis[0].coef()[2] == doctest::Approx(0.444635));

    for(std::size_t i=2; i<5; ++i) {
        CHECK(basis[i].exp()[0] == doctest::Approx(5.033151));
        CHECK(basis[i].exp()[1] == doctest::Approx(1.169596));
        CHECK(basis[i].exp()[2] == doctest::Approx(0.380389));

        CHECK(basis[i].coef()[0] == doctest::Approx(0.155916));
        CHECK(basis[i].coef()[1] == doctest::Approx(0.607684));
        CHECK(basis[i].coef()[2] == doctest::Approx(0.391957));
    }
}