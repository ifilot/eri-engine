#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <eri/cgf.h>

TEST_CASE("CGF can be constructed") {
    eri::CGF a{
        .lx = 0, .ly = 0, .lz = 0,
        .center = {0.0, 0.0, 0.0},
        .exp = {1.24, 0.45},
        .coef = {0.7, 0.3}
    };

    CHECK(a.lx == 0);
    CHECK(a.center[0] == doctest::Approx(0.0));
}