#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <eri/basis/basisset.h>
#include <eri/chem/molecule.h>
#include <eri/math/matrix_builder.h>
#include <eri/ops/ops.h>
#include <eri/utils/executable_dir.h>

static eri::chem::Molecule make_h2o()
{
    eri::chem::Molecule mol;
    mol.add_atom(8, { 0.000000, -0.143207, 0.000000 });
    mol.add_atom(1, { 1.637236,  1.136548, 0.000000 });
    mol.add_atom(1, {-1.637236,  1.136548, 0.000000 });
    return mol;
}

static eri::basis::BasisSet make_basis(const eri::chem::Molecule& mol) {
    eri::basis::BasisSet basis;
    basis.load_from_bse_json(eri::utils::executable_dir() + "/../data/basis/sto-3g.json", mol);
    return basis;
}

const std::array<std::array<double, 7>, 7> h2o_nuclear = {{
    //  O(1s)              O(2s)              O(2px)             O(2py)             O(2pz)             H1(1s)             H2(1s)
    {{ -6.06184779432341827e+01, -7.18308970982290873e+00,  0.00000000000000000e+00, -2.63908815287175561e-17,  0.00000000000000000e+00, -1.19519,              -1.19519              }},
    {{ -7.18308970982290873e+00, -9.05191849590399400e+00,  0.00000000000000000e+00,  0.00000000000000000e+00,  0.00000000000000000e+00, -2.52536,              -2.52536              }},
    {{  0.00000000000000000e+00,  0.00000000000000000e+00, -8.97945049381568339e+00,  0.00000000000000000e+00,  0.00000000000000000e+00, -1.47018,               1.47018              }},
    {{ -2.63908815287175561e-17,  0.00000000000000000e+00,  0.00000000000000000e+00, -8.97945049381568339e+00,  0.00000000000000000e+00, -1.14917,              -1.14917              }},
    {{  0.00000000000000000e+00,  0.00000000000000000e+00,  0.00000000000000000e+00,  0.00000000000000000e+00, -8.97945049381568339e+00,  0.00000000000000000e+00,  0.00000000000000000e+00 }},
    {{ -1.19519,                -2.52536,                -1.47018,                -1.14917,                 0.00000000000000000e+00, -3.76984,              -0.851032             }},
    {{ -1.19519,                -2.52536,                 1.47018,                -1.14917,                 0.00000000000000000e+00, -0.851032,             -3.76984              }}
}};

TEST_CASE("Nuclear attraction matrix for H2O / STO-3G (Huzinaga)") {
    auto mol = make_h2o();
    auto basis = make_basis(mol);
    const std::size_t n = basis.size();
    REQUIRE(n == 7);

    const auto S = eri::math::build_symmetric_matrix<eri::ops::NuclearHuzinaga>(basis, mol[0].position, mol[0].Z);
    constexpr double eps = 1e-5;

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            CHECK_MESSAGE(
                S(i,j) == doctest::Approx(h2o_nuclear[i][j]).epsilon(eps),
                "Mismatch at (" << i << "," << j << "): "
                << S(i,j) << " vs " << h2o_nuclear[i][j]
            );
        }
    }
}

TEST_CASE("Nuclear attraction matrix for H2O / STO-3G (Hellsing)") {
    auto mol = make_h2o();
    auto basis = make_basis(mol);
    const std::size_t n = basis.size();
    REQUIRE(n == 7);

    const auto S = eri::math::build_symmetric_matrix<eri::ops::NuclearHellsing>(basis, mol[0].position, mol[0].Z);
    constexpr double eps = 1e-5;

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            CHECK_MESSAGE(
                S(i,j) == doctest::Approx(h2o_nuclear[i][j]).epsilon(eps),
                "Mismatch at (" << i << "," << j << "): "
                << S(i,j) << " vs " << h2o_nuclear[i][j]
            );
        }
    }
}