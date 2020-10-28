/* align - Align the reads to the reference
    Copyright (C) 2020 Massimiliano Rossi
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file align.cpp
   \brief align.cpp Align the reads to the reference.
   \author Massimiliano Rossi
   \date 13/07/2020
*/

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/io.hpp>

#include <ms_pointers.hpp>

#include <malloc_count.h>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>

#include <ssw_cpp.h>

typedef std::pair<std::string, std::vector<uint8_t>> pattern_t;

// std::vector<pattern_t> read_patterns(std::string filename)
// {
//   // Open File
//   FILE *fd;
//   if ((fd = fopen(filename.c_str(), "r")) == nullptr)
//     error("open() file " + filename + " failed");

//   std::vector<pattern_t> patterns;

//   pattern_t pattern;

//   char c;
//   while (fread(&c, sizeof(char), 1, fd) == 1)
//   {
//     if (c == '>')
//     {
//       if (pattern.second.size() > 0)
//         patterns.push_back(pattern);

//       pattern.first.clear();
//       pattern.second.clear();

//       pattern.first.append(1, c);
//       while (fread(&c, sizeof(char), 1, fd) == 1 && c != '\n')
//         pattern.first.append(1, c);
//     }
//     else
//     {
//       pattern.second.push_back(c);
//       while (fread(&c, sizeof(char), 1, fd) == 1 && c != '\n')
//         pattern.second.push_back(c);
//     }
//   }

//   if (pattern.second.size() > 0)
//     patterns.push_back(pattern);

//   fclose(fd);

//   return patterns;
// }

struct Patterns
{
private:
  std::ifstream input_file;

public:
  using pattern_t = std::pair<std::string, std::vector<uint8_t>>;

  Patterns(std::string input_path) : input_file(input_path) {}

  pattern_t next();
  bool end();

  size_t n_patterns = 0;
};

Patterns::pattern_t Patterns::next()
{
  n_patterns++;
  std::string header, read;
  std::vector<uint8_t> pattern;
  // header
  std::getline(input_file, header);
  // read
  std::getline(input_file, read);

  pattern_t out;
  out.first = header;
  std::copy(read.begin(), read.end(), std::back_inserter(out.second));
  return out;
}

bool Patterns::end()
{
  return input_file.eof();
}

class aligner_t
{
public:
  using SelSd = SelectSdvec<>;
  using DagcSd = DirectAccessibleGammaCode<SelSd>;

  aligner_t(std::string filename, size_t min_len_ = 50) : min_len(min_len_)
  {
    verbose("Building the matching statistics index");
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    std::string filename_ms = filename + ".ms";

    ifstream fs_ms(filename_ms);
    ms.load(fs_ms);

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Matching statistics index construction complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    verbose("Building random access");
    t_insert_start = std::chrono::high_resolution_clock::now();

    std::string filename_slp = filename + ".slp";

    ifstream fs(filename_slp);
    ra.load(fs);

    n = ra.getLen();

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Matching statistics index construction complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
  }

  void align(pattern_t &pattern)
  {
    size_t mem_pos = 0;
    size_t mem_len = 0;
    size_t mem_idx = 0;

    auto pointers = ms.query(pattern.second);
    std::vector<size_t> lengths(pointers.size());
    size_t l = 0;
    for (size_t i = 0; i < pointers.size(); ++i)
    {
      size_t pos = pointers[i];
      while ((i + l) < pattern.second.size() && (pos + l) < n && pattern.second[i + l] == ra.charAt(pos + l))
        ++l;

      lengths[i] = l;
      l = (l == 0 ? 0 : (l - 1));

      // Update MEM
      if (lengths[i] > mem_len)
      {
        mem_len = lengths[i];
        mem_pos = pointers[i];
        mem_idx = i;
      }
    }

    // Align the read
    if (mem_len >= min_len)
    {
      char *str = (char *)malloc(400);

      int32_t maskLen = pattern.second.size() / 2;
      maskLen = maskLen < 15 ? 15 : maskLen;

      std::string query(pattern.second.begin(), pattern.second.end());
      // Extract the context from the reference
      size_t left_occ = (mem_pos > 100 ? mem_pos - 100 : 0);
      size_t len = mem_len + 200;
      ra.expandSubstr(left_occ, len, str);

      // Declares a default Aligner
      StripedSmithWaterman::Aligner aligner;
      // Declares a default filter
      StripedSmithWaterman::Filter filter;
      // Declares an alignment that stores the result
      StripedSmithWaterman::Alignment alignment;
      // Aligns the query to the ref
      aligner.Align(query.c_str(), str, len, filter, &alignment, maskLen);

      std::string ref(str, str + len);

      cout << "===== SSW result =====" << endl;
      cout << ref << endl;
      cout << query << endl;
      cout << "MEM length:\t" << mem_len << endl;
      cout << "MEM position:\t" << mem_pos << endl;
      cout << "Best Smith-Waterman score:\t" << alignment.sw_score << endl
           << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << endl
           << "Reference start:\t" << alignment.ref_begin << endl
           << "Reference end:\t" << alignment.ref_end << endl
           << "Query start:\t" << alignment.query_begin << endl
           << "Query end:\t" << alignment.query_end << endl
           << "Next-best reference end:\t" << alignment.ref_end_next_best << endl
           << "Number of mismatches:\t" << alignment.mismatches << endl
           << "Cigar: " << alignment.cigar_string << endl;
      cout << "======================" << endl;

      aligned_reads++;

      delete str;
    }
  }

  size_t get_aligned_reads()
  {
    return aligned_reads;
  }

protected:
  ms_pointers<> ms;
  SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd> ra;

  size_t min_len = 0;
  size_t aligned_reads = 0;
  size_t n = 0;
};

int
main(int argc, char *const argv[])
{

  Args args;
  parseArgs(argc, argv, args);

  verbose("Construction of the aligner");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  aligner_t aligner(args.filename, 0);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();

  Patterns patterns(args.patterns);

  while (not patterns.end())
  {
    pattern_t pattern = patterns.next();
    aligner.align(pattern);
  }

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Number of aligned reads: ", aligner.get_aligned_reads(), "/", patterns.n_patterns);

  size_t space = 0;
  if (args.memo)
  {
  }

  if (args.store)
  {
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;
}