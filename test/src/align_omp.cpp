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

#include <omp.h>

#include <libgen.h>


typedef std::pair<std::string, std::vector<uint8_t>> pattern_t;

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
    verbose("Loading the matching statistics index");
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    std::string filename_ms = filename + ".ms";

    ifstream fs_ms(filename_ms);
    ms.load(fs_ms);

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Matching statistics index construction complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    verbose("Loading random access");
    t_insert_start = std::chrono::high_resolution_clock::now();

    std::string filename_slp = filename + ".slp";

    ifstream fs(filename_slp);
    ra.load(fs);

    n = ra.getLen();

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Matching statistics index loading complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
  
    verbose("Minimum MEM length: ", min_len);

  }


  bool align(kseq_t *read, FILE* out)
  {
    size_t mem_pos = 0;
    size_t mem_len = 0;
    size_t mem_idx = 0;

    bool aligned = false;

    auto pointers = ms.query(read->seq.s, read->seq.l);
    std::vector<size_t> lengths(pointers.size());
    size_t l = 0;
    for (size_t i = 0; i < pointers.size(); ++i)
    {
      size_t pos = pointers[i];
      while ((i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
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
      uint8_t strand = 0;

      int32_t maskLen = read->seq.l / 2;
      maskLen = maskLen < 15 ? 15 : maskLen;

      size_t ext_l = (read->seq.l - mem_len)*3;

      // Extract the context from the reference
      size_t left_occ = (mem_pos > ext_l ? mem_pos - ext_l : 0);
      size_t len = mem_len + ext_l + (mem_pos > ext_l ? ext_l : ext_l - mem_pos);

      char *str = (char *)malloc(len + 10);
      ra.expandSubstr(left_occ, len, str);

      size_t min_score = 20 + 8 * log(read->seq.l);

      // Declares a default Aligner
      StripedSmithWaterman::Aligner aligner;
      // Declares a default filter
      StripedSmithWaterman::Filter filter;
      // StripedSmithWaterman::Filter filter(true, true, min_score, 32767);
      // Declares an alignment that stores the result
      StripedSmithWaterman::Alignment alignment;
      // Aligns the query to the ref
      aligner.Align(read->seq.s, str, len, filter, &alignment, maskLen);

      // Update alignment method
      alignment.ref_begin += left_occ;
      alignment.ref_end += left_occ;
      alignment.ref_end_next_best += left_occ;

      if(alignment.sw_score >= min_score)
      {
        ssw_write_sam(alignment,"human",read,strand,out);
        aligned = true;
      }

      // aligned_reads++;

      delete str;
    }
    return aligned;
  }

  size_t get_aligned_reads()
  {
    return aligned_reads;
  }


  // Adapted from SSW
  static void ssw_write_sam(StripedSmithWaterman::Alignment &a,
                            const char *ref_seq_name,
                            const kseq_t *read,
                            int8_t strand,
                            FILE* out) // 0: forward aligned ; 1: reverse complement aligned
  {
    // Sam format output
    fprintf(out, "%s\t", read->name.s);
    if (a.sw_score == 0)
      fprintf(out, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
    else
    {
      int32_t c, p;
      uint32_t mapq = -4.343 * log(1 - (double)abs(a.sw_score - a.sw_score_next_best) / (double)a.sw_score);
      mapq = (uint32_t)(mapq + 4.99);
      mapq = mapq < 254 ? mapq : 254;
      if (strand)
        fprintf(out, "16\t");
      else
        fprintf(out, "0\t");
      // TODO: Find the correct reference name.
      fprintf(out, "%s\t%d\t%d\t", ref_seq_name, a.ref_begin + 1, mapq);
      // mismatch = mark_mismatch(a.ref_begin, a.query_begin, a.query_end, ref_num, read_num, read->seq.l, &a->cigar, &a->cigarLen);
      // for (c = 0; c < a->cigarLen; ++c)
      // {
      //   char letter = cigar_int_to_op(a->cigar[c]);
      //   uint32_t length = cigar_int_to_len(a->cigar[c]);
      //   fprintf(out, "%lu%c", (unsigned long)length, letter);
      // }
      fprintf(out, "%s", a.cigar_string.c_str());
      fprintf(out, "\t*\t0\t0\t");
      fprintf(out, "%s", read->seq.s);
      fprintf(out, "\t");
      if (read->qual.s && strand)
      {
        for (p = read->qual.l - 1; p >= 0; --p)
          fprintf(out, "%c", read->qual.s[p]);
      }
      else if (read->qual.s)
        fprintf(out, "%s", read->qual.s);
      else
        fprintf(out, "*");
      fprintf(out, "\tAS:i:%d", a.sw_score);
      fprintf(out, "\tNM:i:%d\t", a.mismatches);
      if (a.sw_score_next_best > 0)
        fprintf(out, "ZS:i:%d\n", a.sw_score_next_best);
      else
        fprintf(out, "\n");
    }
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

  aligner_t aligner(args.filename, args.l);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();
  
  size_t n_reads = 0;
  size_t n_aligned_reads = 0;

  std::string base_name = basename(args.filename.data());

// #pragma omp parallel
// {
#pragma omp parallel for schedule(static)
  for(size_t i = 1; i <= 64; ++i)
  {
    size_t n_reads_p = 0;
    size_t n_aligned_reads_p = 0;

    gzFile fp;
    kseq_t *seq;
    int l;

    std::string file_path = args.patterns + "_" + std::to_string(i) + ".fa";

    // Open out file
    std::string filename = file_path + "_" + base_name + "_" +std::to_string(args.l) +  ".sam";
    FILE *sam_fd;

    if ((sam_fd = fopen(filename.c_str(), "w")) == nullptr)
      error("open() file " + filename + " failed");

    fp = gzopen(file_path.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0)
    {
      if(aligner.align(seq,sam_fd))
        n_aligned_reads_p ++;
      n_reads_p++;

      // std::cout << "\rSequenced patterns on block " << i << " : "
      //           << n_reads_p << std::flush;
    }

    verbose("Number of aligned reads block ", i, " : ", n_aligned_reads_p, "/", n_reads_p);
    kseq_destroy(seq);
    gzclose(fp);
    fclose(sam_fd);

  }

// #pragma omp critical
//   {
//     // Count the total number of reads.
//     n_reads += n_reads_p;
//     // Count the total number of aligned reads
//     n_aligned_reads += n_aligned_reads_p;
//   }
// }




  // TODO: Merge the SAM files.

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Number of aligned reads: ", n_aligned_reads, "/", n_reads);

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