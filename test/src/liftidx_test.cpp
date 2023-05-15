/* lifting_test - Testing levioSAM lifting
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
   \file lifting_test.cpp
   \brief lifting_test.cpp Testing levioSAM lifting.
   \author Massimiliano Rossi
   \date 19/11/2021
*/

#include <iostream>

#define VERBOSE

#include <kseq.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread);

#include <liftidx.hpp>
#include <common.hpp>

#include <malloc_count.h>


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

  verbose("Loading the lifting index");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  verbose("Index filename: " + args.filename);

  liftidx lift;
  std::ifstream in(args.filename + lift.get_file_extension());
  lift.load(in);
  in.close();

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Lifting index loading complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("SAM Header");
  verbose(lift.to_sam());

  verbose("Lifting positions");
  verbose("Lifting position ", 48126296, ": ", lift.lift(48126296));
  verbose("Lifting position ", 68126296, ": ", lift.lift(68126296));
  verbose("Lifting position ", 130126296, ": ", lift.lift(130126296));
  verbose("Lifting position ", 180126296, ": ", lift.lift(180126296));
  verbose("Lifting position ", 230126296, ": ", lift.lift(230126296));
  // for(size_t i = 48126296; i < 48126395; ++i)
  //   verbose("Lifting position ", i, ": ", lift.lift(i));

  const auto names = lift.get_names();
  size_t len = 0;
  for (size_t i = 0; i < names.size(); len+=lift.length(i)+10, ++i)
  {
    const auto pos = len;
    const auto lpos = lift.lift(pos);
    const auto lft_ref = lift.index(lpos);
    verbose("Lifting position ", pos, " in ", names[i], ": ", lpos, " corresponding to: ", lft_ref.second + 1, " in ", lft_ref.first);
    if (names[i] == "HG01358_H2_chr21")
    {
      const auto pos = len + 29312633;
      const auto lpos = lift.lift(pos);
      const auto lft_ref = lift.index(lpos);
      verbose("Lifting position ", pos, " in ", names[i], ": ", lpos, " corresponding to: ", lft_ref.second + 1, " in ", lft_ref.first);
    }
  }

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

  return 0;
}