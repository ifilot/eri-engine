#include <benchmark/benchmark.h>

#include <eri/basis/basisset.h>
#include <eri/chem/molecule.h>
#include <eri/one_electron/overlap.h>
#include <eri/enums.h>
#include <eri/math/dense_square_matrix.h>
#include <eri/utils/executable_dir.h>

static eri::basis::BasisSet make_basis() {
    eri::chem::Molecule mol;
    // Carbon ring (Z = 6)
    mol.add_atom(6, { 1.397,  0.000, 0.0});
    mol.add_atom(6, { 0.699,  1.211, 0.0});
    mol.add_atom(6, {-0.699,  1.211, 0.0});
    mol.add_atom(6, {-1.397,  0.000, 0.0});
    mol.add_atom(6, {-0.699, -1.211, 0.0});
    mol.add_atom(6, { 0.699, -1.211, 0.0});

    // Hydrogens (Z = 1)
    mol.add_atom(1, { 2.479,  0.000, 0.0});
    mol.add_atom(1, { 1.240,  2.148, 0.0});
    mol.add_atom(1, {-1.240,  2.148, 0.0});
    mol.add_atom(1, {-2.479,  0.000, 0.0});
    mol.add_atom(1, {-1.240, -2.148, 0.0});
    mol.add_atom(1, { 1.240, -2.148, 0.0});

    eri::basis::BasisSet basis;
    basis.load_from_bse_json(eri::utils::executable_dir() + "/../data/basis/sto-3g.json", mol);
    return basis;
}

static void BM_Overlap_Huzinaga(benchmark::State& state) {
    static eri::basis::BasisSet basis = make_basis();

    for (auto _ : state) {
        eri::math::DenseSquareMatrix S(basis.size());

        for (std::size_t i = 0; i < basis.size(); ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                double v = eri::one_electron::overlap(
                    basis[i], basis[j],
                    eri::enums::OverlapMethod::Huzinaga
                );
                S(i, j) = v;
                S(j, i) = v;
            }
        }

        benchmark::DoNotOptimize(S.data());
    }
}

static void BM_Overlap_Hellsing(benchmark::State& state) {
    static eri::basis::BasisSet basis = make_basis();

    for (auto _ : state) {
        eri::math::DenseSquareMatrix S(basis.size());

        for (std::size_t i = 0; i < basis.size(); ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                double v = eri::one_electron::overlap(
                    basis[i], basis[j],
                    eri::enums::OverlapMethod::Hellsing
                );
                S(i, j) = v;
                S(j, i) = v;
            }
        }

        benchmark::DoNotOptimize(S.data());
    }
}

BENCHMARK(BM_Overlap_Huzinaga);
BENCHMARK(BM_Overlap_Hellsing);
BENCHMARK_MAIN();
