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

struct file_info {
  std::ifstream is;
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
void output_counts(cpp_array<file_info>& files) {
  cpp_array<reader_type> readers(files.size());
  typedef jellyfish::mer_heap::heap<mer_dna, reader_type> heap_type;
  typedef typename heap_type::const_item_t heap_item;
  heap_type heap(files.size());

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
    do {
      counts[head->it_ - base] = head->val_;
      heap.pop();
      if(head->it_->next())
        heap.push(*head->it_);
      head = heap.head();
    } while(head->key_ == key && heap.is_not_empty());
	// Assembly counts in slot 1, read counts in slot 2
	coverage_count[counts[0]][counts[1]]++;
  }
  for (const std::pair<uint64_t, std::unordered_map<uint64_t, uint64_t>>& x : coverage_count) {
	  for (const std::pair<uint64_t, uint64_t>& y : x.second) {
		  std::cout << x.first << "\t" << y.first << "\t" << y.second << std::endl;
	  }
  }
}

int main(int argc, char *argv[])
{
  // Check number of input files
    if (argc < 2) {err::die(err::msg() << "usage: " <<
							argv[0] << "assembly_file read_file" <<
							"\n\nARGUMENTS:\n" <<
							"\tassembly_file\t\tjellyfish database from genome assembly\n" <<
							"\tread_file\t\tjellyfish database from short read data\n");							
    }

  // Read the header of each input files and do sanity checks.
  cpp_array<file_info> files(argc - 1);
  common_info cinfo = read_headers(argc - 1, argv + 1, files);
  mer_dna::k(cinfo.key_len / 2);

  if(cinfo.format == binary_dumper::format)
    output_counts<binary_reader>(files);
  else if(cinfo.format == text_dumper::format)
    output_counts<text_reader>(files);
  else
    err::die(err::msg() << "Format '" << cinfo.format << "' not supported\n");

  return 0;
}
