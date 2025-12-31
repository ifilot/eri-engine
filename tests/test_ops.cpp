#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <eri/basis/basisset.h>
#include <eri/chem/molecule.h>
#include <eri/math/dense_square_matrix.h>
#include <eri/math/matrix_builder.h>
#include <eri/ops/ops.h>
#include <eri/utils/executable_dir.h>

// --- small, robust test system ---
static eri::chem::Molecule make_h2()
{
    eri::chem::Molecule mol;
    mol.add_atom(1, {0.0, 0.0, 0.0});
    mol.add_atom(1, {0.0, 0.0, 1.4});
    return mol;
}

static eri::basis::BasisSet make_basis()
{
    eri::basis::BasisSet basis;
    auto mol = make_h2();

    // adjust path if you resolve relative to executable
    basis.load_from_bse_json(eri::utils::executable_dir() + "/../data/basis/sto-3g.json", mol);
    return basis;
}

TEST_CASE("Overlap matrix via operator policies")
{
    const auto basis = make_basis();
    const std::size_t n = basis.size();

    REQUIRE(n > 0);

    const auto S_def =
        eri::math::build_symmetric_matrix<eri::ops::Overlap>(basis);

    const auto S_huz =
        eri::math::build_symmetric_matrix<eri::ops::OverlapHuzinaga>(basis);

    const auto S_hel =
        eri::math::build_symmetric_matrix<eri::ops::OverlapHellsing>(basis);

    // --- basic invariants ---
    CHECK(S_def.size() == n);
    CHECK(S_huz.size() == n);
    CHECK(S_hel.size() == n);

    for (std::size_t i = 0; i < n; ++i) {
        CHECK(S_def(i,i) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(S_huz(i,i) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(S_hel(i,i) == doctest::Approx(1.0).epsilon(1e-6));

        for (std::size_t j = 0; j < n; ++j) {
            CHECK(S_def(i,j) == doctest::Approx(S_def(j,i)));
            CHECK(S_huz(i,j) == doctest::Approx(S_huz(j,i)));
            CHECK(S_hel(i,j) == doctest::Approx(S_hel(j,i)));
        }
    }

    // --- default policy must match explicit Huzinaga ---
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            CHECK(
                S_def(i,j)
                ==
                doctest::Approx(S_huz(i,j)).epsilon(1e-12)
            );
        }

    // --- Huzinaga vs Hellsing should agree numerically ---
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            CHECK(
                S_huz(i,j)
                ==
                doctest::Approx(S_hel(i,j)).epsilon(1e-8)
            );
        }
}

TEST_CASE("Kinetic matrix via operator policies")
{
    const auto basis = make_basis();
    const std::size_t n = basis.size();

    REQUIRE(n > 0);

    const auto T_def =
        eri::math::build_symmetric_matrix<eri::ops::Kinetic>(basis);

    const auto T_huz =
        eri::math::build_symmetric_matrix<eri::ops::KineticHuzinaga>(basis);

    const auto T_hel =
        eri::math::build_symmetric_matrix<eri::ops::KineticHellsing>(basis);

    // --- symmetry & positivity ---
    for (std::size_t i = 0; i < n; ++i) {
        CHECK(T_def(i,i) > 0.0);
        CHECK(T_huz(i,i) > 0.0);
        CHECK(T_hel(i,i) > 0.0);

        for (std::size_t j = 0; j < n; ++j) {
            CHECK(T_def(i,j) == doctest::Approx(T_def(j,i)));
            CHECK(T_huz(i,j) == doctest::Approx(T_huz(j,i)));
            CHECK(T_hel(i,j) == doctest::Approx(T_hel(j,i)));
        }
    }

    // --- default policy must match explicit Huzinaga ---
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            CHECK(
                T_def(i,j)
                ==
                doctest::Approx(T_huz(i,j)).epsilon(1e-12)
            );
        }
}
