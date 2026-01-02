#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <eri/chem/molecule.h>
#include <eri/basis/basisset.h>
#include <eri/one_electron/overlap.h>
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

static eri::basis::BasisSet make_basis_h2o() {
    eri::chem::Molecule mol = make_h2o();
    eri::basis::BasisSet basis;
    basis.load_from_bse_json(eri::utils::executable_dir() + "/../data/basis/sto-3g.json", mol);
    return basis;
}

TEST_CASE("Overlap Huzinaga vs Hellsing agree (STO-3G 1s)") {
    eri::basis::CGF a(0,0,0, {0,0,0},
                     {3.42525091, 0.62391373, 0.16885540},
                     {0.15432897, 0.53532814, 0.44463454});
    eri::basis::CGF b = a;

    const double S_h = eri::one_electron::overlap(a, b, eri::enums::OverlapMethod::Huzinaga);
    const double S_s = eri::one_electron::overlap(a, b, eri::enums::OverlapMethod::Hellsing);

    CHECK(S_h == doctest::Approx(S_s).epsilon(1e-12));
}

TEST_CASE("Overlap Huzinaga H2") {
    eri::basis::CGF a(0,0,0, {0,0,0},
                     {3.42525091, 0.62391373, 0.16885540},
                     {0.15432897, 0.53532814, 0.44463454});
    eri::basis::CGF b(0,0,0, {0,0,1.4},
                     {3.42525091, 0.62391373, 0.16885540},
                     {0.15432897, 0.53532814, 0.44463454});

    const double S11 = eri::one_electron::overlap(a, a, eri::enums::OverlapMethod::Huzinaga);
    const double S12 = eri::one_electron::overlap(a, b, eri::enums::OverlapMethod::Huzinaga);
    const double S22 = eri::one_electron::overlap(b, b, eri::enums::OverlapMethod::Huzinaga);

    CHECK(S11 == doctest::Approx(1.0));
    CHECK(S12 == doctest::Approx(0.65931845));
    CHECK(S22 == doctest::Approx(1.0));
}

TEST_CASE("Overlap Hellsing H2") {
    eri::basis::CGF a(0,0,0, {0,0,0},
                     {3.42525091, 0.62391373, 0.16885540},
                     {0.15432897, 0.53532814, 0.44463454});
    eri::basis::CGF b(0,0,0, {0,0,1.4},
                     {3.42525091, 0.62391373, 0.16885540},
                     {0.15432897, 0.53532814, 0.44463454});

    const double S11 = eri::one_electron::overlap(a, a, eri::enums::OverlapMethod::Hellsing);
    const double S12 = eri::one_electron::overlap(a, b, eri::enums::OverlapMethod::Hellsing);
    const double S22 = eri::one_electron::overlap(b, b, eri::enums::OverlapMethod::Hellsing);

    CHECK(S11 == doctest::Approx(1.0));
    CHECK(S12 == doctest::Approx(0.65931845));
    CHECK(S22 == doctest::Approx(1.0));
}

TEST_CASE("Self-overlap p-type self-overlap (STO-3G Li 2p Huzinaga)") {
    eri::basis::CGF px(
        1, 0, 0,
        {0.0, 0.0, 0.0},
        {0.636290, 0.147860, 0.048089},
        {0.155916, 0.607684, 0.391957}
    );

    const double S = eri::one_electron::overlap(px, px, eri::enums::OverlapMethod::Huzinaga);

    CHECK(S == doctest::Approx(1.0).epsilon(1e-6));
}

TEST_CASE("Self-overlap p-type self-overlap (STO-3G Li 2p Hellsing)") {
    eri::basis::CGF px(
        1, 0, 0,
        {0.0, 0.0, 0.0},
        {0.636290, 0.147860, 0.048089},
        {0.155916, 0.607684, 0.391957}
    );

    const double S = eri::one_electron::overlap(px, px, eri::enums::OverlapMethod::Hellsing);

    CHECK(S == doctest::Approx(1.0).epsilon(1e-6));
}

TEST_CASE("Overlap p-type orthogonality (p_x | p_y)") {
    eri::basis::CGF px(
        1, 0, 0,
        {0.0, 0.0, 0.0},
        {0.636290, 0.147860, 0.048089},
        {0.155916, 0.607684, 0.391957}
    );

    eri::basis::CGF py(
        0, 1, 0,
        {0.0, 0.0, 0.0},
        {0.636290, 0.147860, 0.048089},
        {0.155916, 0.607684, 0.391957}
    );

    const double S =
        eri::one_electron::overlap(px, py, eri::enums::OverlapMethod::Huzinaga);

    CHECK(std::abs(S) < 1e-12);
}

const std::array<std::array<double, 7>, 7> h2o_overlap = {{
    //  O(1s)        O(2s)        O(2px)       O(2py)       O(2pz)       H1(1s)        H2(1s)
    {{ 1.00000000,  0.23670400,  0.00000000,  8.60875e-18, 0.00000000,  0.0384367,   0.0384367 }},
    {{ 0.23670400,  1.00000000,  0.00000000, -2.46745e-17,0.00000000,  0.386337,    0.386337  }},
    {{ 0.00000000,  0.00000000,  1.00000000,  0.00000000, 0.00000000,  0.268493,   -0.268493  }},
    {{ 8.60875e-18,-2.46745e-17, 0.00000000,  1.00000000, 0.00000000,  0.209869,    0.209869  }},
    {{ 0.00000000,  0.00000000,  0.00000000,  0.00000000, 1.00000000,  0.00000000,  0.00000000}},
    {{ 0.0384367,   0.386337,    0.268493,    0.209869,    0.00000000,  1.00000000,  0.181995  }},
    {{ 0.0384367,   0.386337,   -0.268493,    0.209869,    0.00000000,  0.181995,    1.00000000}}
}};

TEST_CASE("Full overlap matrix for H2O / STO-3G (Huzinaga)") {
    auto basis = make_basis_h2o();
    const std::size_t n = basis.size();
    REQUIRE(n == 7);

    const auto S = eri::math::build_symmetric_matrix<eri::ops::one_electron::OverlapHuzinaga>(basis);
    constexpr double eps = 1e-5;

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            CHECK_MESSAGE(
                S(i,j) == doctest::Approx(h2o_overlap[i][j]).epsilon(eps),
                "Mismatch at (" << i << "," << j << "): "
                << S(i,j) << " vs " << h2o_overlap[i][j]
            );
        }
    }
}

TEST_CASE("Full overlap matrix for H2O / STO-3G (Hellsing)") {
    auto basis = make_basis_h2o();
    const std::size_t n = basis.size();
    REQUIRE(n == 7);

    const auto S = eri::math::build_symmetric_matrix<eri::ops::one_electron::OverlapHellsing>(basis);
    constexpr double eps = 1e-5;

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            CHECK_MESSAGE(
                S(i,j) == doctest::Approx(h2o_overlap[i][j]).epsilon(eps),
                "Mismatch at (" << i << "," << j << "): "
                << S(i,j) << " vs " << h2o_overlap[i][j]
            );
        }
    }
}