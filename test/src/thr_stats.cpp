/* matching_statistics - Computes the matching statistics from BWT and Thresholds
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
   \file matching_statistics.cpp
   \brief matching_statistics.cpp Computes the matching statistics from BWT and Thresholds.
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

#define maxr(a,b) ((a) = std::max((a),(b)))
#define minr(a,b) ((a) = std::min((a), (b)))

//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
  std::string filename = "";
  size_t w = 10;             // sliding window size and its default
  bool store = false;        // store the data structure in the file
  bool memo = false;         // print the memory usage
  bool csv = false;          // print stats on stderr in csv format
  bool rle = false;          // outpt RLBWT
  std::string patterns = ""; // path to patterns file
  size_t l = 25;             // minumum MEM length
  size_t th = 1;             // number of threads
  bool is_fasta = false;     // read a fasta file
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string usage("usage: " + std::string(argv[0]) + " infile [-s store] [-m memo] [-c csv] [-p patterns] [-f fasta] [-r rle] [-t threads] [-l len]\n\n" +
                    "Computes the pfp data structures of infile, provided that infile.parse, infile.dict, and infile.occ exists.\n" +
                    "  wsize: [integer] - sliding window size (def. 10)\n" +
                    "  store: [boolean] - store the data structure in infile.pfp.ds. (def. false)\n" +
                    "   memo: [boolean] - print the data structure memory usage. (def. false)\n" +
                    "  fasta: [boolean] - the input file is a fasta file. (def. false)\n" +
                    "    rle: [boolean] - output run length encoded BWT. (def. false)\n" +
                    "pattens: [string]  - path to patterns file.\n" +
                    "    len: [integer] - minimum MEM lengt (def. 25)\n" +
                    " thread: [integer] - number of threads (def. 1)\n" +
                    "    csv: [boolean] - print the stats in csv form on strerr. (def. false)\n");

  std::string sarg;
  while ((c = getopt(argc, argv, "w:smcfl:rhp:t:")) != -1)
  {
    switch (c)
    {
    case 'w':
      sarg.assign(optarg);
      arg.w = stoi(sarg);
      break;
    case 's':
      arg.store = true;
      break;
    case 'm':
      arg.memo = true;
      break;
    case 'c':
      arg.csv = true;
      break;
    case 'r':
      arg.rle = true;
      break;
    case 'p':
      arg.patterns.assign(optarg);
      break;
    case 'l':
      sarg.assign(optarg);
      arg.l = stoi(sarg);
      break;
    case 't':
      sarg.assign(optarg);
      arg.th = stoi(sarg);
      break;
    case 'f':
      arg.is_fasta = true;
      break;
    case 'h':
      error(usage);
    case '?':
      error("Unknown option.\n", usage);
      exit(1);
    }
  }
  // the only input parameter is the file name
  if (argc == optind + 1)
  {
    arg.filename.assign(argv[optind]);
  }
  else
  {
    error("Invalid number of arguments\n", usage);
  }
}

//********** end argument options ********************


int main(int argc, char *const argv[])
{

  Args args;
  parseArgs(argc, argv, args);

  verbose("Loading the matching statistics index");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  std::string bwt_fname = args.filename + ".bwt";

  std::string bwt_heads_fname = bwt_fname + ".heads";
  std::ifstream ifs_heads(bwt_heads_fname);
  if (!ifs_heads.is_open())
    error("open() file " + bwt_heads_fname + " failed");
  std::string bwt_len_fname = bwt_fname + ".len";
  std::ifstream ifs_len(bwt_len_fname);
  if (!ifs_len.is_open())
    error("open() file " + bwt_len_fname + " failed");
  ms_rle_string_sd bwt = ms_rle_string_sd(ifs_heads, ifs_len);

  ifs_heads.close();
  ifs_len.close();

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Matching statistics index construction complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Reading thresholds from file");

  t_insert_start = std::chrono::high_resolution_clock::now();


  std::string tmp_filename = args.filename + std::string(".thr_pos");

  size_t n = uint64_t(bwt.size());
  int log_n = bitsize(uint64_t(bwt.size()));

  struct stat filestat;
  FILE *fd;

  if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
    error("open() file " + tmp_filename + " failed");

  int fn = fileno(fd);
  if (fstat(fn, &filestat) < 0)
    error("stat() file " + tmp_filename + " failed");

  if (filestat.st_size % THRBYTES != 0)
    error("invilid file " + tmp_filename);

  size_t length = filestat.st_size / THRBYTES;
  size_t threshold;

  sdsl::int_vector<> thresholds = sdsl::int_vector<>(length, 0, log_n);

  size_t pos = 0;

  size_t max_thr = 0;
  long long max_off = 0;
  long long min_off = n;

  uint8_t char_at_max = 0;

  for (size_t i = 0; i < length; ++i)
  {
    if ((fread(&threshold, THRBYTES, 1, fd)) != 1)
      error("fread() file " + tmp_filename + " failed");

    // verbose(pos, " ",bwt.head_of(i), " ", bwt.run_at(i));
    if(threshold > 0)
    {
      uint8_t c = bwt.head_of(i);
      size_t pred = bwt.select(bwt.rank(pos-1,c)-1,c); 
      size_t mid_int = (pos - pred + 1)>>1;
      assert(threshold > pred);
      if(threshold < pred)
      {
        error("Threshold before the beginning of the run", threshold, " ", pred);
      }
      threshold = threshold - pred;

      long long off = mid_int - threshold;

      if(max_thr < threshold)
      {
        char_at_max = c;
        maxr(max_thr,threshold);
      }

      maxr(max_off,off);
      minr(min_off,off);

    }



    thresholds[i] = threshold;

    pos += bwt.run_at(i);
  }

  fclose(fd);

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());


  verbose("Log n : ", log_n);
  verbose("Max threshold: ", max_thr);
  verbose("Char at max: ", char_at_max);
  int log_thr = bitsize(max_thr);
  verbose("Log threshold: ", log_thr);

  verbose("Max offset   : ", max_off);
  verbose("Min offset   : ", min_off);
  int log_off = bitsize((size_t)(max_off - min_off + 1));
  verbose("Log offset   : ", log_off);

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

  size_t space = 0;
  if (args.memo)
  {
    verbose("Thresholds size (bytes): ", space);
  }

  if (args.store)
  {
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;
}