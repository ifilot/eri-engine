#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include <eri/engine/engine.h>

TEST_CASE("Engine smoke test") {
    eri::Engine eng;

    CHECK_NOTHROW(eng.compute());
}