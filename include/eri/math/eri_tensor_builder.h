#pragma once

#include <vector>
#include <cstddef>

#include <eri/math/symmetric_eri_tensor.h>
#include <eri/two_electron/eri.h>
#include <eri/two_electron/eri_index.h>
#include <eri/math/fgamma.h>

namespace eri::math {

struct ERIJob {
    std::size_t i, j, k, l;
};

template <typename ERIEvaluator>
SymmetricERITensor build_eri_tensor(const basis::BasisSet& basis)
{
    const std::size_t n = basis.size();

    // TODO: establish proper value for interpolation table from basis set
    eri::math::init_Fgamma_interp_table(4);

    // construct evaluator ONCE
    ERIEvaluator eri_eval(basis);

    SymmetricERITensor eritensor(n);

    std::vector<ERIJob> jobs;
    jobs.reserve(n * n * n * n / 8);

    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j <= i; ++j)
            for (std::size_t k = 0; k < n; ++k)
                for (std::size_t l = 0; l <= k; ++l)
                    if (eri::two_electron::pair_index(i,j) >=
                        eri::two_electron::pair_index(k,l))
                        jobs.push_back({i,j,k,l});

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (const auto& job : jobs) {
        eritensor(job.i, job.j, job.k, job.l) =
            eri_eval(
                basis[job.i], basis[job.j],
                basis[job.k], basis[job.l]
            );
    }

    return eritensor;
}

} // namespace eri::math
