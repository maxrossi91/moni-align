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

extern "C" {
#include <xerrors.h>
}

#include <iostream>

#define VERBOSE
#define MTIME


#include <sdsl/io.hpp>

#define _REALIGN
#include <aligner_ksw2.hpp>
#include <align_reads_dispatcher.hpp>

#include <common.hpp>
#include <stacktrace.hpp>
#include <malloc_count.h>


//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
  std::string filename = "";
  size_t w = 10; // sliding window size and its default
  bool store = false; // store the data structure in the file
  bool memo  = false; // print the memory usage
  bool csv   = false; // print stats on stderr in csv format
  bool rle   = false; // outpt RLBWT
  std::string patterns = ""; // path to patterns file
  std::string mate1 = "";    // path to file with #1 mates paired with mate2.
  std::string mate2 = "";    // path to file with #2 mates paired with mate1.
  std::string output   = ""; // output file prefix
  size_t b = 512;       // size of the batch of read to be processed
  size_t l = 25; // minumum MEM length
  size_t th = 1; // number of threads
  bool is_fasta = false; // read a fasta file
  bool shaped_slp = false; // use shaped slp
  size_t ext_len = 100;      // Extension length
  // size_t top_k = 1;       // Report the top_k alignments
  bool filter_dir = true; // Use MEMs average length to filter the orientation of the reads
  double dir_thr = 50.0; // Use MEMs average length distance to filter the orientation of the reads

  bool filter_seeds = true; // Filter seed if occurs more than threshold
  size_t n_seeds_thr = 5000;   // Filter seed if occurs more than threshold
  

  // ksw2 parameters
  int8_t smatch = 2;      // Match score default
  int8_t smismatch = 4;   // Mismatch score default
  int8_t gapo = 4;        // Gap open penalty
  int8_t gapo2 = 13;      // Gap open penalty
  int8_t gape = 2;        // Gap extension penalty
  int8_t gape2 = 1;       // Gap extension penalty
  // int end_bonus = 400;    // Bonus to add at the extension score to declare the alignment

  // int w = -1;             // Band width
  // int zdrop = -1;         // Zdrop enable
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string usage("usage: " + std::string(argv[0]) + " infile [-p patterns] [-1 mate1] [-2 mate2] [-o output] [-t threads] [-b batch] [-l len] [-q shaped_slp]  [-L ext_l] [-A smatch] [-B smismatc] [-O gapo] [-E gape] [-d dir_en] [-s seeds_en] [-D dir_thr] [-S seeds_thr]\n\n" +
                    "Align the reads in the pattern against the reference index in infile.\n" +
                    "   pattens: [string]  - path to patterns file.\n" +
                    "     mate1: [string]  - path to file with #1 mates paired with mate2.\n" +
                    "     mate2: [string]  - path to file with #2 mates paired with mate1.\n" +
                    "    output: [string]  - output file prefix.\n" +
                    "   threads: [integer] - number of threads (def. " + std::to_string(arg.th) + ")\n" +
                    "     batch: [integer] - number of batches per therad pool (def. " + std::to_string(arg.b) + ")\n" + 
                    "       len: [integer] - minimum MEM lengt (def. " + std::to_string(arg.l) + ")\n" +
                    "shaped_slp: [boolean] - use shaped slp. (def. " + std::to_string(arg.shaped_slp) + ")\n" +
                    "     ext_l: [integer] - length of reference substring for extension (def. " + std::to_string(arg.ext_len) + ")\n" +
                    "   dir_dis: [boolean] - enable direction filtering (def. " + std::to_string(arg.filter_dir) + ")\n" +
                    "   dir_thr: [float]   - direction filtering threshold (def. " + std::to_string(arg.dir_thr) + ")\n" +
                    " seeds_dis: [boolean] - enable seed filtering (def. " + std::to_string(arg.filter_dir) + ")\n" +
                    " seeds_thr: [float]   - seed filtering threshold (def. " + std::to_string(arg.n_seeds_thr) + ")\n" +
                    "    smatch: [integer] - match score value (def. " + std::to_string(arg.smatch) + ")\n" +
                    " smismatch: [integer] - mismatch penalty value (def. " + std::to_string(arg.smismatch) + ")\n" +
                    "      gapo: [integer] - gap open penalty value (def. " + std::to_string(arg.gapo) + "," + std::to_string(arg.gapo2) + ")\n" +
                    "      gape: [integer] - gap extension penalty value (def. " + std::to_string(arg.gape) + "," + std::to_string(arg.gape2) + ")\n");

  std::string sarg;
  char* s;
  while ((c = getopt(argc, argv, "ql:hp:o:t:1:2:b:A:B:O:E:L:dsD:S:")) != -1)
  {
    switch (c)
    {
    case 'p':
      arg.patterns.assign(optarg);
      break;
    case 'o':
      arg.output.assign(optarg);
      break;
    case '1':
      arg.mate1.assign(optarg);
      break;
    case '2':
      arg.mate2.assign(optarg);
      break;
    case 'l':
      sarg.assign(optarg);
      arg.l = stoi(sarg);
      break;
    case 'b':
      sarg.assign(optarg);
      arg.b = stoi(sarg);
      break;
    case 't':
      sarg.assign(optarg);
      arg.th = stoi(sarg);
      break;
    case 'L':
      sarg.assign(optarg);
      arg.ext_len = stoi(sarg);
      break;
    case 'A':
      sarg.assign(optarg);
      arg.smatch = stoi(sarg);
      break;
    case 'B':
      sarg.assign(optarg);
      arg.smismatch = stoi(sarg);
      break;
    case 'd':
      arg.filter_dir = false;
      break;
    case 's':
      arg.filter_seeds = false;
      break;
    case 'D':
      sarg.assign(optarg);
      arg.dir_thr = stoi(sarg);
      break;
    case 'S':
      sarg.assign(optarg);
      arg.n_seeds_thr = stoi(sarg);
      break;
    case 'O':
      arg.gapo = arg.gapo2 = strtol(optarg, &s, 10);
      if (*s == ',') arg.gapo2 = strtol(s+1, &s, 10);
      break;
    case 'E':
      arg.gape = arg.gape2 = strtol(optarg, &s, 10);
      if (*s == ',') arg.gape2 = strtol(s+1, &s, 10);
      break;
    case 'q':
      arg.shaped_slp = true;
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


template<typename aligner_t>
typename aligner_t::config_t configurer(Args &args){
  typename aligner_t::config_t config;
  
  config.min_len    = args.l;           // Minimum MEM length
  config.ext_len    = args.ext_len;     // Extension length

  // ksw2 parameters
  config.smatch     = args.smatch;      // Match score default
  config.smismatch  = args.smismatch;   // Mismatch score default
  config.gapo       = args.gapo;        // Gap open penalty
  config.gapo2      = args.gapo2;       // Gap open penalty
  config.gape       = args.gape;        // Gap extension penalty
  config.gape2      = args.gape2;       // Gap extension penalty

  return config;
}

template<typename aligner_t>
void dispatcher(Args &args){

  verbose("Construction of the aligner");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  aligner_t aligner(args.filename, configurer<aligner_t>(args));

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();

  std::string base_name = basename(args.filename.data());


  if(args.patterns != "")
  {
    std::string sam_filename = args.patterns + "_" + base_name + "_" + std::to_string(args.l) + ".sam";
    if(args.output != "")
      sam_filename = args.output;

    verbose("Output file: ", sam_filename);

    // if (is_gzipped(args.patterns))
    // {
    //   verbose("The input is gzipped - forcing single thread alignment.");
    //   args.th = 1;
    // }
    if (args.th == 1)
      st_align<aligner_t>(&aligner, args.patterns, sam_filename, args.b);
    else
      mt_align<aligner_t>(&aligner, args.patterns, sam_filename, args.th, args.b);
  }
  else
  {
    std::string sam_filename = args.mate1 + "_" + base_name + "_" + std::to_string(args.l) + ".sam";
    if(args.output != "")
      sam_filename = args.output;

    verbose("Output file: ", sam_filename);

    // if (is_gzipped(args.mate1) or is_gzipped(args.mate2))
    // {
    //   verbose("The input is gzipped - forcing single thread alignment.");
    //   args.th = 1;
    // }
    if (args.th == 1)
      st_align<aligner_t>(&aligner, args.mate1, sam_filename, args.b, args.mate2 );
    else
      mt_align<aligner_t>(&aligner, args.mate1, sam_filename, args.th, args.b, args.mate2);
  }

  // TODO: Merge the SAM files.

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

}

int main(int argc, char *const argv[])
{

  Args args;
  parseArgs(argc, argv, args);

  enable_stacktrace();

  if(args.shaped_slp){
    dispatcher<aligner<shaped_slp_t, ms_pointers<>>>(args);
  }else{
    dispatcher<aligner<plain_slp_t, ms_pointers<>>>(args);
  }

  MTIME_TSAFE_REPORT_ALL;

  return 0;
}