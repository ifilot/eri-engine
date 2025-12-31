#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"
#include "eri/basis/cgf.h"

TEST_CASE("CGF construction and basic properties") {
    eri::basis::CGF a(
        0, 0, 0,
        {0.0, 0.0, 0.0},
        {1.24, 0.45},
        {0.7, 0.3}
    );

    CHECK(a.lx() == 0);
    CHECK(a.ly() == 0);
    CHECK(a.lz() == 0);

    CHECK(a.ctr()[0] == doctest::Approx(0.0));
    CHECK(a.ctr()[1] == doctest::Approx(0.0));
    CHECK(a.ctr()[2] == doctest::Approx(0.0));

    CHECK(a.exp().size() == 2);
    CHECK(a.coef().size() == 2);
    CHECK(a.norm().size() == 2);
}

TEST_CASE("CGF construction with STO-3G (H 1s)") {
    // STO-3G parameters for Hydrogen 1s
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

    eri::basis::CGF cgf(
        0, 0, 0,
        {0.0, 0.0, 0.0},
        exps,
        coefs
    );

    // Basic properties
    CHECK(cgf.lx() == 0);
    CHECK(cgf.ly() == 0);
    CHECK(cgf.lz() == 0);

    CHECK(cgf.exp().size() == 3);
    CHECK(cgf.coef().size() == 3);
    CHECK(cgf.norm().size() == 3);

    // Coefficients must be preserved exactly
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(cgf.coef()[i] == doctest::Approx(coefs[i]));
    }

    // Primitive normalization constants must be positive and finite
    for (double N : cgf.norm()) {
        CHECK(std::isfinite(N));
        CHECK(N > 0.0);
    }

    CHECK(cgf.overlap() == doctest::Approx(1.0));
}

TEST_CASE("CGF construction with STO-3G (2p_x)") {
    // STO-3G parameters for a 2p shell
    const std::vector<double> exps = {
        2.94124940,
        0.68348310,
        0.22228990
    };

    const std::vector<double> coefs = {
        0.15591627,
        0.60768372,
        0.39195739
    };

    // p_x function: lx = 1, ly = 0, lz = 0
    eri::basis::CGF cgf(
        1, 0, 0,
        {0.0, 0.0, 0.0},
        exps,
        coefs
    );

    // Angular momentum
    CHECK(cgf.lx() == 1);
    CHECK(cgf.ly() == 0);
    CHECK(cgf.lz() == 0);

    // Sizes
    CHECK(cgf.exp().size() == 3);
    CHECK(cgf.coef().size() == 3);
    CHECK(cgf.norm().size() == 3);

    // Coefficients must be preserved (basis-set semantics)
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(cgf.coef()[i] == doctest::Approx(coefs[i]));
    }

    // Primitive normalization constants:
    // must be positive, finite, and different from s-type norms
    for (double N : cgf.norm()) {
        CHECK(std::isfinite(N));
        CHECK(N > 0.0);
    }

    CHECK(cgf.overlap() == doctest::Approx(1.0));
}