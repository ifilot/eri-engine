#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <eri/one_electron/overlap.h>

TEST_CASE("Overlap Huzinaga vs Hellsing agree (STO-3G 1s)") {
    eri::basis::CGF a(0,0,0, {0,0,0},
                     {3.42525091, 0.62391373, 0.16885540},
                     {0.15432897, 0.53532814, 0.44463454});
    eri::basis::CGF b = a;

    const double S_h = eri::one_electron::overlap(a, b, eri::one_electron::OverlapMethod::Huzinaga);
    const double S_s = eri::one_electron::overlap(a, b, eri::one_electron::OverlapMethod::Hellsing);

    CHECK(S_h == doctest::Approx(S_s).epsilon(1e-12));
}

TEST_CASE("Overlap Huzinaga H2") {
    eri::basis::CGF a(0,0,0, {0,0,0},
                     {3.42525091, 0.62391373, 0.16885540},
                     {0.15432897, 0.53532814, 0.44463454});
    eri::basis::CGF b(0,0,0, {0,0,1.4},
                     {3.42525091, 0.62391373, 0.16885540},
                     {0.15432897, 0.53532814, 0.44463454});

    const double S11 = eri::one_electron::overlap(a, a, eri::one_electron::OverlapMethod::Huzinaga);
    const double S12 = eri::one_electron::overlap(a, b, eri::one_electron::OverlapMethod::Huzinaga);
    const double S22 = eri::one_electron::overlap(b, b, eri::one_electron::OverlapMethod::Huzinaga);

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

    const double S11 = eri::one_electron::overlap(a, a, eri::one_electron::OverlapMethod::Hellsing);
    const double S12 = eri::one_electron::overlap(a, b, eri::one_electron::OverlapMethod::Hellsing);
    const double S22 = eri::one_electron::overlap(b, b, eri::one_electron::OverlapMethod::Hellsing);

    CHECK(S11 == doctest::Approx(1.0));
    CHECK(S12 == doctest::Approx(0.65931845));
    CHECK(S22 == doctest::Approx(1.0));
}

TEST_CASE("Overlap p-type self-overlap (STO-3G Li 2p Huzinaga)") {
    eri::basis::CGF px(
        1, 0, 0,
        {0.0, 0.0, 0.0},
        {0.636290, 0.147860, 0.048089},
        {0.155916, 0.607684, 0.391957}
    );

    const double S = eri::one_electron::overlap(px, px, eri::one_electron::OverlapMethod::Huzinaga);

    CHECK(S == doctest::Approx(1.0).epsilon(1e-6));
}

TEST_CASE("Overlap p-type self-overlap (STO-3G Li 2p Hellsing)") {
    eri::basis::CGF px(
        1, 0, 0,
        {0.0, 0.0, 0.0},
        {0.636290, 0.147860, 0.048089},
        {0.155916, 0.607684, 0.391957}
    );

    const double S = eri::one_electron::overlap(px, px, eri::one_electron::OverlapMethod::Hellsing);

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
        eri::one_electron::overlap(px, py, eri::one_electron::OverlapMethod::Huzinaga);

    CHECK(std::abs(S) < 1e-12);
}

TEST_CASE("Overlap p-type Huzinaga vs Hellsing (displaced centers)") {
    eri::basis::CGF px_A(
        1, 0, 0,
        {0.0, 0.0, 0.0},
        {0.636290, 0.147860, 0.048089},
        {0.155916, 0.607684, 0.391957}
    );

    eri::basis::CGF px_B(
        1, 0, 0,
        {0.0, 0.0, 1.4},
        {0.636290, 0.147860, 0.048089},
        {0.155916, 0.607684, 0.391957}
    );

    const double S_h =
        eri::one_electron::overlap(px_A, px_B, eri::one_electron::OverlapMethod::Huzinaga);

    const double S_s =
        eri::one_electron::overlap(px_A, px_B, eri::one_electron::OverlapMethod::Hellsing);

    CHECK(S_h == doctest::Approx(S_s).epsilon(1e-12));
}