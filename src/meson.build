# general
sources = files(
    'basis/basisset.cpp',
    'basis/cgf.cpp',
    'engine/engine.cpp',
)

# integrals
sources += files(
    'one_electron/overlap.cpp',
    'one_electron/kinetic.cpp',
    'one_electron/nuclear.cpp',
    'two_electron/eri.cpp',
)

# math
sources += files(
    'math/fgamma.cpp',
    'math/fgamma_interpolation.cpp',
    'math/boys.cpp'
)

sources += files(
    'cog/hellsing_cache.cpp',
)

inc = include_directories(
    '.',
    '../include',
)

jsoncpp_dep = dependency('jsoncpp', required: true)
openmp_dep = dependency('openmp', required: false)

eri_deps = [jsoncpp_dep]
if openmp_dep.found()
  eri_deps += openmp_dep
endif

eri_core = library(
    'eri-core',
    sources,
    include_directories: inc,
    dependencies: eri_deps,
    install: true
)

install_subdir(
    'include/eri',
    install_dir: get_option('includedir')
)
