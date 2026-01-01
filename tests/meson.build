test_inc = include_directories('..', '../include')

test_sources = {
  'engine':   'test_engine.cpp',
  'cgf':      'test_cgf.cpp',
  'overlap':  'test_overlap.cpp',
  'kinetic':  'test_kinetic.cpp',
  'nuclear':  'test_nuclear.cpp',
  'basisset': 'test_basisset.cpp',
  'ops':      'test_ops.cpp',
  'eri':      'test_eri.cpp',
}

foreach name, src : test_sources
  exe = executable(
    'test-' + name,
    src,
    link_with: eri_core,
    include_directories: test_inc
  )

  test(name, exe)
endforeach

h2o_eri_data = configure_file(
  input: 'data/h2o_eri.txt',
  output: 'h2o_eri.txt',
  copy: true
)