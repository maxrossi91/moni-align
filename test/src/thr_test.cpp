/* thr_test - Test if plain and compressed version are equals
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
   \file thr_test.cpp
   \brief thr_test.cpp Test if plain and compressed version are equals.
   \author Massimiliano Rossi
   \date 17/11/2020
*/

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/io.hpp>

#include <ms_pointers.hpp>

#include <thresholds_ds.hpp>

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

  verbose("Building the thresholds from file");

  t_insert_start = std::chrono::high_resolution_clock::now();


  thr_plain<> plain(args.filename,&bwt);
  verbose("Plain thresholds construction complete");
  thr_compressed<> compressed(args.filename,&bwt);

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Checking the thresholds");

  t_insert_start = std::chrono::high_resolution_clock::now();

  for (size_t i = 0; i < bwt.number_of_runs(); ++i)
  {
    if(plain[i] != compressed[i])
      error("Plain and compressed are different in position ", i , " . Plain: ", plain[i], " Compressed: ", compressed[i]);
  }

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());


  

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