#pragma once

#include <eri/enums.h>
#include <eri/basis/cgf.h>
#include <eri/cog/hellsing_cache.h>

using eri::basis::CGF;

namespace eri::two_electron {

double eri_huzinaga(const CGF& a, const CGF& b, const CGF& c, const CGF& d);

double eri_hellsing(const CGF& a, const CGF& b, const CGF& c, const CGF& d);

double eri_hellsing_cached(const CGF& a, const CGF& b, const CGF& c, const CGF& d,
                           const eri::cog::HellsingCacheTable1D& cache);

} // namespace eri::one_electron