#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <eri/basis/cgf.h>
#include <eri/one_electron/nuclear.h>
#include <eri/enums.h>

TEST_CASE("Nuclear attraction self-integral s-type STO-3G is negative")
{
    // STO-3G hydrogen 1s CGF
    eri::basis::CGF a(
        0, 0, 0,
        {0.0, 0.0, 0.0},
        {3.42525091, 0.62391373, 0.16885540},
        {0.15432897, 0.53532814, 0.44463454}
    );

    // nucleus at same center (hydrogen)
    const std::array<double,3> C{0.0, 0.0, 0.0};

    const double V = eri::one_electron::nuclear(
        a, a, C,
        eri::enums::NuclearMethod::Huzinaga
    );

    CHECK(std::isfinite(V));
    CHECK(V < 0.0);
}

TEST_CASE("Nuclear attraction Huzinaga vs Hellsing agree (CGF level)")
{
    eri::basis::CGF a(
        0, 0, 0,
        {0.0, 0.0, 0.0},
        {3.42525091, 0.62391373, 0.16885540},
        {0.15432897, 0.53532814, 0.44463454}
    );

    // place nucleus slightly off-center to avoid trivial symmetry
    const std::array<double,3> C{0.0, 0.0, 0.2};

    const double V_huz = eri::one_electron::nuclear(
        a, a, C,
        eri::enums::NuclearMethod::Huzinaga
    );

    const double V_hel = eri::one_electron::nuclear(
        a, a, C,
        eri::enums::NuclearMethod::Hellsing
    );

    CHECK(std::isfinite(V_huz));
    CHECK(std::isfinite(V_hel));

    CHECK(V_huz < 0.0);
    CHECK(V_hel < 0.0);

    CHECK(V_huz == doctest::Approx(V_hel).epsilon(1e-8));
}
