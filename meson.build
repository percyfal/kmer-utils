project('jellyfish-utils', ['c', 'cpp'], default_options: ['c_std=c99', 'cpp_std=c++11'])

cc = meson.get_compiler('c')

extra_c_args = [
    '-Wall', '-Wextra', '-Werror', '-Wpedantic', '-W',
    '-Wmissing-prototypes', # '-Wstrict-prototypes', # does not work with R.h
    '-Wconversion', '-Wshadow', '-Wpointer-arith', '-Wcast-align',
    '-Wcast-qual', '-Wwrite-strings', '-Wnested-externs',
    '-fshort-enums', '-fno-common']
  

# pkg = import('pkgconfig')
jellyfishdep = dependency('jellyfish-2.0')
#jellyfishdep = cc.find_library('jellyfish-2.0')

executable('kmer_count_pairs',
	   sources: 'kmer_count_pairs.cc', dependencies : jellyfishdep, install: true)
