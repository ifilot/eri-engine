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

subdir('src')
subdir('tests')