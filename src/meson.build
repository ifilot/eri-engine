sources = files(
    'engine/engine.cpp',
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
