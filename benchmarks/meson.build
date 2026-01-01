benchmark_dep = dependency('benchmark', required: true)

bench_deps = [benchmark_dep]

if openmp_dep.found()
  bench_deps += openmp_dep
endif

bench_inc = include_directories(
  '..',
  '../include'
)

benchmark_sources = {
  'overlap': 'benchmark_overlap.cpp',
  'eri':     'benchmark_eri.cpp',
}

foreach name, src : benchmark_sources
  executable(
    'benchmark-' + name,
    src,
    link_with: eri_core,
    dependencies: bench_deps,
    include_directories: bench_inc,
    install: false
  )
endforeach
