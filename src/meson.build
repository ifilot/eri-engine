# general
sources = files(
    'basis/cgf.cpp',
    'engine/engine.cpp',
)

# overlap integrals
sources += files(
    'one_electron/overlap_huzinaga.cpp',
    'one_electron/overlap_hellsing.cpp',
    'one_electron/overlap.cpp',
)

inc = include_directories(
    '.',
    '../include',
)

eri_core = library(
    'eri-core',
    sources,
    include_directories: inc,
    install: true
)

install_subdir(
    'include/eri',
    install_dir: get_option('includedir')
)
