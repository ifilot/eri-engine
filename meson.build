project(
    'eri-engine',
    'cpp',
    version: '0.1.0',
    default_options: [
        'cpp_std=c++20',
        'buildtype=release'
    ]
)

version = meson.project_version()
version_parts = version.split('.')

conf = configuration_data()
conf.set('ERI_VERSION', version)
conf.set('ERI_VERSION_MAJOR', version_parts[0])
conf.set('ERI_VERSION_MINOR', version_parts[1])
conf.set('ERI_VERSION_PATCH', version_parts[2])

configure_file(
    input: 'include/eri/version.h.in',
    output: 'version.h',
    configuration: conf
)

subdir('data/basis')
subdir('src')
subdir('tests')
subdir('benchmarks')

run_target(
  'bench-eri',
  command: [
    'sh', '-c',
    'meson test -C ' + meson.current_build_dir() +
    ' --print-errorlogs && ' +
    join_paths(meson.current_build_dir(), 'benchmarks', 'benchmark-eri') +
    ' --benchmark_repetitions=3' +
    ' --benchmark_report_aggregates_only=true' +
    ' --benchmark_time_unit=s'
  ],
  env: ['OMP_NUM_THREADS=20']
)