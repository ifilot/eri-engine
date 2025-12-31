benchmark_dep = dependency('benchmark', required: true)

bench_inc = include_directories(
  '..',
  '../include'
)

benchmark_sources = {
  'overlap': 'benchmark_overlap.cpp',
}

foreach name, src : benchmark_sources
  executable(
    'benchmark-' + name,
    src,
    link_with: eri_core,
    dependencies: benchmark_dep,
    include_directories: bench_inc,
    install: false
  )
endforeach
