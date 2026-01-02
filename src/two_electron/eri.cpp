#include <eri/two_electron/eri.h>

#include "two_electron/eri_contracted.h"
#include "two_electron/primitives/huzinaga.h"
#include "two_electron/primitives/hellsing.h"
#include "two_electron/primitives/hellsing_cached.h"

using eri::two_electron::primitives::PrimitiveHuzinaga;
using eri::two_electron::primitives::PrimitiveHellsing;
using eri::two_electron::primitives::PrimitiveHellsingCached;

namespace eri::two_electron {

double eri_huzinaga(const CGF& a, const CGF& b, const CGF& c, const CGF& d) {
    return eri_contracted(a, b, c, d, PrimitiveHuzinaga{});
}

double eri_hellsing(const CGF& a, const CGF& b, const CGF& c, const CGF& d) {
    return eri_contracted(a, b, c, d, PrimitiveHellsing{});
}

double eri_hellsing_cached(const CGF& a, const CGF& b, const CGF& c, const CGF& d,
                                  const eri::cog::HellsingCacheTable1D& cache) {
    return eri_contracted(a, b, c, d, PrimitiveHellsingCached{cache});
}

}