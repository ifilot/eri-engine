#pragma once

#include <vector>
#include <cstddef>

#include <eri/math/symmetric_eri_tensor.h>
#include <eri/two_electron/eri.h>
#include <eri/two_electron/eri_index.h>

namespace eri::math {

struct ERIJob {
    std::size_t i, j, k, l;
};

template <typename Operator>
SymmetricERITensor build_eri_tensor(const basis::BasisSet& basis)
{
    const std::size_t n = basis.size();
    SymmetricERITensor eri(n);

    const std::size_t npairs = n * (n + 1) / 2;

    // collect all jobs
    std::vector<ERIJob> jobs;
    jobs.reserve(npairs * (npairs + 1) / 2);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            const std::size_t ij = eri::two_electron::pair_index(i, j);

            for (std::size_t k = 0; k < n; ++k) {
                for (std::size_t l = 0; l <= k; ++l) {
                    const std::size_t kl = eri::two_electron::pair_index(k, l);

                    if (ij >= kl) {
                        jobs.push_back({i, j, k, l});
                    }
                }
            }
        }
    }

    // evaluate and populate
    for (const auto& job : jobs) {
        eri(job.i, job.j, job.k, job.l) =
            Operator::eval(basis[job.i], basis[job.j],
                           basis[job.k], basis[job.l]);
    }

    return eri;
}

} // namespace eri::math
