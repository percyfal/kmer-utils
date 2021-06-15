project('kmer-utils', ['c', 'cpp'], default_options: ['c_std=c99', 'cpp_std=c++11'])

cc = meson.get_compiler('c')

extra_c_args = [ '-Wall', '-Wextra', '-Werror', '-Wpedantic', '-W',
    '-Wmissing-prototypes', '-Wconversion', '-Wshadow',
    '-Wpointer-arith', '-Wcast-align', '-Wcast-qual',
    '-Wwrite-strings', '-Wnested-externs', '-fshort-enums',
    '-fno-common']
  

jellyfishdep = dependency('jellyfish-2.0', version: '>=2.3.0', method: 'pkg-config')

executable('kmer_count_pairs',
	   sources: 'kmer_count_pairs.cc', dependencies : jellyfishdep, install: true)
