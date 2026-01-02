#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <eri/chem/molecule.h>
#include <eri/basis/basisset.h>
#include <eri/two_electron/eri.h>
#include <eri/utils/executable_dir.h>
#include <eri/math/eri_tensor_builder.h>
#include <eri/two_electron/eri_index.h>
#include <eri/ops/ops.h>

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

TEST_CASE("ERI Huzinaga (STO-3G H2)") {
    auto basis = make_basis_h2();

    double T1111 = eri::two_electron::eri_huzinaga(basis[0], basis[0], basis[0], basis[0]);
    double T1122 = eri::two_electron::eri_huzinaga(basis[0], basis[0], basis[1], basis[1]);
    double T1112 = eri::two_electron::eri_huzinaga(basis[0], basis[0], basis[0], basis[1]);
    double T2121 = eri::two_electron::eri_huzinaga(basis[1], basis[0], basis[1], basis[0]);
    double T1222 = eri::two_electron::eri_huzinaga(basis[0], basis[1], basis[1], basis[1]);
    double T2211 = eri::two_electron::eri_huzinaga(basis[1], basis[1], basis[0], basis[0]);

    const double epsilon = 1e-6;

    CHECK(T1111 == doctest::Approx(0.774605).epsilon(epsilon));
    CHECK(T1122 == doctest::Approx(0.569675).epsilon(epsilon));
    CHECK(T1112 == doctest::Approx(0.444107).epsilon(epsilon));
    CHECK(T2121 == doctest::Approx(0.297028).epsilon(epsilon));

    CHECK(T1222 == doctest::Approx(T1112).epsilon(epsilon));
    CHECK(T1122 == doctest::Approx(T2211).epsilon(epsilon));
}

TEST_CASE("ERI Huzinaga (STO-3G H2O - selected)") {
    auto basis = make_basis_h2o();
    const double epsilon = 1e-6;
    CHECK(eri::two_electron::eri_huzinaga(basis[0], basis[0], basis[0], basis[0]) == doctest::Approx(4.785069).epsilon(epsilon));
    CHECK(eri::two_electron::eri_huzinaga(basis[6], basis[6], basis[6], basis[6]) == doctest::Approx(0.774605).epsilon(epsilon));
    CHECK(eri::two_electron::eri_huzinaga(basis[1], basis[2], basis[3], basis[4]) == doctest::Approx(0).epsilon(epsilon));
    CHECK(eri::two_electron::eri_huzinaga(basis[2], basis[2], basis[5], basis[5]) == doctest::Approx(0.469864).epsilon(epsilon));
    CHECK(eri::two_electron::eri_huzinaga(basis[2], basis[2], basis[6], basis[6]) == doctest::Approx(0.469864).epsilon(epsilon));
}

TEST_CASE("ERI Huzinaga (STO-3G H2O)") {
    auto basis = make_basis_h2o();
    const std::size_t N = basis.size();
    auto tetensor = eri::math::build_eri_tensor<eri::ops::two_electron::ERIHuzinaga>(basis);

    const std::string fname = eri::utils::executable_dir() + "/h2o_eri.txt";

    std::ifstream in(fname);
    REQUIRE(in.good());

    constexpr double tol = 1e-6;

    double ref;
    std::size_t count = 0;

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            for (std::size_t k = 0; k < N; ++k) {
                for (std::size_t l = 0; l < N; ++l) {
                    in >> ref;
                    const double val = tetensor(i,j,k,l);
                    CHECK(val == doctest::Approx(ref).epsilon(tol));
                    ++count;
                }
            }
        }
    }

    INFO("Verified unique ERIs: " << count);
}

TEST_CASE("ERI Hellsing (STO-3G H2O)") {
    auto basis = make_basis_h2o();
    const std::size_t N = basis.size();
    auto tetensor = eri::math::build_eri_tensor<eri::ops::two_electron::ERIHellsing>(basis);

    const std::string fname = eri::utils::executable_dir() + "/h2o_eri.txt";

    std::ifstream in(fname);
    REQUIRE(in.good());

    constexpr double tol = 1e-6;

    double ref;
    std::size_t count = 0;

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            for (std::size_t k = 0; k < N; ++k) {
                for (std::size_t l = 0; l < N; ++l) {
                    in >> ref;
                    const double val = tetensor(i,j,k,l);
                    CHECK(val == doctest::Approx(ref).epsilon(tol));
                    ++count;
                }
            }
        }
    }

    INFO("Verified unique ERIs: " << count);
}