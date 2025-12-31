#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <eri/basis/basisset.h>
#include <eri/chem/molecule.h>
#include <eri/math/matrix_builder.h>
#include <eri/ops/ops.h>
#include <eri/utils/executable_dir.h>

#include <cmath>

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

const std::array<std::array<double, 7>, 7> h2o_kinetic = {{
    //  O(1s)        O(2s)        O(2px)       O(2py)       O(2pz)       H1(1s)        H2(1s)
    {{ 29.00320000, -0.16801100,  0.00000000,  5.61892e-17, 0.00000000, -0.00841086, -0.00841086 }},
    {{ -0.16801100,  0.80812800,  0.00000000, -2.26480e-17, 0.00000000,  0.07063230,  0.07063230 }},
    {{  0.00000000,  0.00000000,  2.52873000,  0.00000000, 0.00000000,  0.14722400, -0.14722400 }},
    {{  5.61892e-17,-2.26480e-17, 0.00000000,  2.52873000, 0.00000000,  0.11507900,  0.11507900 }},
    {{  0.00000000,  0.00000000,  0.00000000,  0.00000000, 2.52873000,  0.00000000,  0.00000000 }},
    {{ -0.00841086,  0.07063230,  0.14722400,  0.11507900, 0.00000000,  0.76003200, -0.00394818 }},
    {{ -0.00841086,  0.07063230, -0.14722400,  0.11507900, 0.00000000, -0.00394818,  0.76003200 }}
}};

TEST_CASE("Kinetic energy matrix for H2O / STO-3G (Huzinaga)")
{
    auto basis = make_basis_h2o();

    const std::size_t n = basis.size();
    REQUIRE(n == 7); // STO-3G H2O has 7 Cartesian basis functions

    const auto T = eri::math::build_symmetric_matrix<eri::ops::KineticHuzinaga>(basis);

    for (std::size_t i = 0; i < n; ++i) {
        CHECK(std::isfinite(T(i,i)));
        CHECK(T(i,i) > 0.0);          // kinetic energy is positive
        for (std::size_t j = 0; j < n; ++j) {
            CHECK(T(i,j) == doctest::Approx(T(j,i)));
        }
    }

    constexpr double eps = 1e-5;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            CHECK_MESSAGE(
                T(i,j) == doctest::Approx(h2o_kinetic[i][j]).epsilon(eps),
                "Mismatch at (" << i << "," << j << "): "
                << T(i,j) << " vs " << h2o_kinetic[i][j]
            );
        }
    }
}

TEST_CASE("Kinetic energy matrix for H2O / STO-3G (Hellsing)")
{
    auto basis = make_basis_h2o();

    const std::size_t n = basis.size();
    REQUIRE(n == 7); // STO-3G H2O has 7 Cartesian basis functions

    const auto T = eri::math::build_symmetric_matrix<eri::ops::KineticHellsing>(basis);

    for (std::size_t i = 0; i < n; ++i) {
        CHECK(std::isfinite(T(i,i)));
        CHECK(T(i,i) > 0.0);          // kinetic energy is positive
        for (std::size_t j = 0; j < n; ++j) {
            CHECK(T(i,j) == doctest::Approx(T(j,i)));
        }
    }

    constexpr double eps = 1e-5;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            CHECK_MESSAGE(
                T(i,j) == doctest::Approx(h2o_kinetic[i][j]).epsilon(eps),
                "Mismatch at (" << i << "," << j << "): "
                << T(i,j) << " vs " << h2o_kinetic[i][j]
            );
        }
    }
}
