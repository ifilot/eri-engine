#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <eri/one_electron/kinetic.h>

TEST_CASE("Kinetic self-integral s-type STO-3G is positive") {
    eri::basis::CGF a(
         0,0,0,
        {0,0,0},
        {3.42525091, 0.62391373, 0.16885540},
        {0.15432897, 0.53532814, 0.44463454}
    );

    const double T = eri::one_electron::kinetic(a, a);
    CHECK(T > 0.0);
}