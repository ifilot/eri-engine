#pragma once

#include <eri/one_electron/overlap.h>
#include <eri/one_electron/kinetic.h>
#include <eri/one_electron/nuclear.h>
#include <eri/two_electron/eri.h>
#include <eri/cog/hellsing_cache.h>

namespace eri::ops::two_electron {

struct ERIHuzinaga {
    explicit ERIHuzinaga(const basis::BasisSet&) {
        // nothing to do
    }

    double operator()(const basis::CGF& a,
                      const basis::CGF& b,
                      const basis::CGF& c,
                      const basis::CGF& d) const
    {
        return eri::two_electron::eri_huzinaga(a, b, c, d);
    }
};

struct ERIHellsing {
    explicit ERIHellsing(const basis::BasisSet&) {
        // nothing to do
    }

    double operator()(const basis::CGF& a,
                      const basis::CGF& b,
                      const basis::CGF& c,
                      const basis::CGF& d) const
    {
        return eri::two_electron::eri_hellsing(a, b, c, d);
    }
};

struct ERIHellsingCached {
    eri::cog::HellsingCacheTable cache;

    explicit ERIHellsingCached(const basis::BasisSet& basis)
        : cache(compute_lmax(basis))
    {}

    double operator()(const basis::CGF& a,
                      const basis::CGF& b,
                      const basis::CGF& c,
                      const basis::CGF& d) const
    {
        return eri::two_electron::eri_hellsing_cached(a, b, c, d, cache);
    }

private:
    static int compute_lmax(const basis::BasisSet& basis)
    {
        int lmax = 0;
        for (const auto& shell : basis) {
            lmax = std::max(
                lmax,
                shell.lx() + shell.ly() + shell.lz()
            );
        }
        return lmax;
    }
};

}