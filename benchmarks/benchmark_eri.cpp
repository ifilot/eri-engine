#include <benchmark/benchmark.h>

#include <eri/chem/molecule.h>
#include <eri/basis/basisset.h>
#include <eri/two_electron/eri.h>
#include <eri/utils/executable_dir.h>
#include <eri/math/eri_tensor_builder.h>
#include <eri/two_electron/eri_index.h>
#include <eri/ops/ops.h>

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

static void BM_ERI_Huzinaga(benchmark::State& state) {
    static eri::basis::BasisSet basis = make_basis();
    eri::math::init_Fgamma_interp_table(4);

    for (auto _ : state) {
        auto tetensor = eri::math::build_eri_tensor<eri::ops::two_electron::ERIHuzinaga>(basis);
        benchmark::DoNotOptimize(tetensor.data());
        benchmark::ClobberMemory();
    }
}

static void BM_ERI_Hellsing(benchmark::State& state) {
    static eri::basis::BasisSet basis = make_basis();
    eri::math::init_Fgamma_interp_table(4);

    for (auto _ : state) {
        auto tetensor = eri::math::build_eri_tensor<eri::ops::two_electron::ERIHellsing>(basis);
        benchmark::DoNotOptimize(tetensor.data());
        benchmark::ClobberMemory();
    }
}

static void BM_ERI_HellsingCached(benchmark::State& state) {
    static eri::basis::BasisSet basis = make_basis();
    eri::math::init_Fgamma_interp_table(4);

    for (auto _ : state) {
        auto tetensor = eri::math::build_eri_tensor<eri::ops::two_electron::ERIHellsingCached>(basis);
        benchmark::DoNotOptimize(tetensor.data());
        benchmark::ClobberMemory();
    }
}

BENCHMARK(BM_ERI_Huzinaga);
BENCHMARK(BM_ERI_Hellsing);
BENCHMARK(BM_ERI_HellsingCached);
BENCHMARK_MAIN();
