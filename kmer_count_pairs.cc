/**
 * @file   kmer_count_pairs.cc
 * @author Per Unneberg
 * @date   Thu Feb 18 19:54:02 2021
 *
 * @brief Count kmer occurrences from two jellyfish databases
 *
 * Based on
 * https://github.com/gmarcais/Jellyfish/tree/master/examples/count_in_file
 * and
 * https://github.com/tallgran/conifer-hifi-assembly/tree/master/src/kmer_compare
 *
 */

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <string>
#include <unordered_map>

#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/mer_heap.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/rectangular_binary_matrix.hpp>
#include <jellyfish/cpp_array.hpp>

namespace err = jellyfish::err;

using jellyfish::file_header;
using jellyfish::RectangularBinaryMatrix;
using jellyfish::mer_dna;
using jellyfish::cpp_array;
typedef std::unique_ptr<binary_reader>           binary_reader_ptr;
typedef std::unique_ptr<text_reader>             text_reader_ptr;
typedef jellyfish::cooperative::hash_counter<jellyfish::mer_dna>                  mer_hash_t;

struct file_info {std::ifstream is;
  file_header   header;

  file_info(const char* path) :
    is(path),
    header(is)
  { }
};


struct common_info {
  unsigned int            key_len;
  size_t                  max_reprobe_offset;
  size_t                  size;
  unsigned int            out_counter_len;
  std::string             format;
  RectangularBinaryMatrix matrix;

  common_info(RectangularBinaryMatrix&& m) : matrix(std::move(m))
  { }
};

common_info read_headers(int argc, char* input_files[], cpp_array<file_info>& files) {
  // Read first file
  files.init(0, input_files[0]);
  if(!files[0].is.good())
    err::die(err::msg() << "Failed to open input file '" << input_files[0] << "'");

  file_header& h = files[0].header;
  common_info res(h.matrix());
  res.key_len            = h.key_len();
  res.max_reprobe_offset = h.max_reprobe_offset();
  res.size               = h.size();
  res.format = h.format();
  size_t reprobes[h.max_reprobe() + 1];
  h.get_reprobes(reprobes);
  res.out_counter_len = h.counter_len();

  // Other files must match
  for(int i = 1; i < argc; i++) {
    files.init(i, input_files[i]);
    file_header& nh = files[i].header;
    if(!files[i].is.good())
      err::die(err::msg() << "Failed to open input file '" << input_files[i] << "'");
    if(res.format != nh.format())
      err::die(err::msg() << "Can't compare files with different formats (" << res.format << ", " << nh.format() << ")");
    if(res.key_len != nh.key_len())
      err::die(err::msg() << "Can't compare hashes of different key lengths (" << res.key_len << ", " << nh.key_len() << ")");
    if(res.max_reprobe_offset != nh.max_reprobe_offset())
      err::die("Can't compare hashes with different reprobing strategies");
    if(res.size != nh.size())
      err::die(err::msg() << "Can't compare hash with different size (" << res.size << ", " << nh.size() << ")");
    if(res.matrix != nh.matrix())
      err::die("Can't compare hash with different hash function");
  }

  return res;
}

template<typename reader_type>
void output_counts(cpp_array<file_info>& files, mer_hash_t& mer_hash_, char *outfile) {
  cpp_array<reader_type> readers(files.size());
  typedef jellyfish::mer_heap::heap<mer_dna, reader_type> heap_type;
  typedef typename heap_type::const_item_t heap_item;
  heap_type heap(files.size());
  std::ofstream outfile_;

  // Prime heap
  for(size_t i = 0; i < files.size(); ++i) {
    readers.init(i, files[i].is, &files[i].header);
    if(readers[i].next())
      heap.push(readers[i]);
  }

  heap_item          head      = heap.head();
  mer_dna            key;
  const int          num_files = files.size();
  const reader_type* base      = &readers[0];
  uint64_t           counts[num_files];

  std::unordered_map<uint64_t, std::unordered_map<uint64_t, uint64_t>> coverage_count;

  while(heap.is_not_empty()) {
    key = head->key_;
    memset(counts, '\0', sizeof(uint64_t) * num_files);
    // Heap consists of ordered array from both files; collect all
    // keys of given type before moving on
    do {
      counts[head->it_ - base] = head->val_;
      heap.pop();
      if(head->it_->next())
        heap.push(*head->it_);
      head = heap.head();
    } while(head->key_ == key && heap.is_not_empty());

    // Assembly counts in slot 1, read counts in slot 2
    coverage_count[counts[0]][counts[1]]++;
    mer_hash_.add(key, counts[0]);
  }

  outfile_.open(outfile);
  for (const std::pair<uint64_t, std::unordered_map<uint64_t, uint64_t>>& x : coverage_count) {
    for (const std::pair<uint64_t, uint64_t>& y : x.second) {
      outfile_ << x.first << "\t" << y.first << "\t" << y.second << std::endl;
    }
  }
}

int main(int argc, char *argv[])
{
  jellyfish::file_header header;
  header.fill_standard();
  header.set_cmdline(argc, argv);

  // Get options
  int c;
  bool saveMers = false;
  while (1) {
    int option_index = 0;
    static struct option long_options[] = {
      {"savemers",  no_argument,       0,  0 },
      {0,         0,                 0,  0 }
    };
    c = getopt_long(argc, argv, "m", long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'm':
      saveMers = true;
      break;
    case '?':
      break;
    }
  };

  // Check number of arguments
  if ((argc - optind) != 3) {err::die(err::msg() << "usage: " <<
                           argv[0] << " [-m] assembly_file read_file out_prefix" <<
                           "\n\nARGUMENTS:\n" <<
                           "\tassembly_file\t\tjellyfish database from genome assembly\n" <<
                           "\tread_file\t\tjellyfish database from short read data\n" <<
                           "\tout_prefix\t\toutput prefix\n");
  }

  // Read the header of each input files and do sanity checks.
  cpp_array<file_info> files(2);
  common_info cinfo = read_headers(2, argv + optind, files);
  mer_dna::k(cinfo.key_len / 2);

  std::unique_ptr<jellyfish::dumper_t<mer_array> > dumper;

  char outfile[1024];
  strcpy(outfile, argv[argc - 1]);
  strcat(outfile, "_mers.jf");

  mer_hash_t mer_hash(cinfo.size, cinfo.key_len, 24, 1, 126);
  dumper.reset(new binary_dumper(4, mer_hash.key_len(), 1, outfile, &header));
  dumper->one_file(true);
  mer_hash.dumper(dumper.get());

  // table output file name
  char tablefile[1024];
  strcpy(tablefile, argv[3]);
  strcat(tablefile, ".tsv");
  if(cinfo.format == binary_dumper::format)
    output_counts<binary_reader>(files, mer_hash, tablefile);
  else if (cinfo.format == text_dumper::format)
    output_counts<text_reader>(files, mer_hash, tablefile);
  else
    err::die(err::msg() << "Format '" << cinfo.format << "' not supported\n");
  if (saveMers)
    dumper->dump(mer_hash.ary());
  return 0;
}
