# general
sources = files(
    'basis/basisset.cpp',
    'basis/cgf.cpp',
    'engine/engine.cpp',
)

# one-electron integrals
sources += files(
    'one_electron/overlap_huzinaga.cpp',
    'one_electron/overlap_hellsing.cpp',
    'one_electron/overlap.cpp',
    'one_electron/kinetic.cpp',
    'one_electron/nuclear.cpp',
)

inc = include_directories(
    '.',
    '../include',
)

jsoncpp_dep = dependency('jsoncpp', required: true)

eri_core = library(
    'eri-core',
    sources,
    include_directories: inc,
    dependencies: [jsoncpp_dep],
    install: true
)

install_subdir(
    'include/eri',
    install_dir: get_option('includedir')
)
