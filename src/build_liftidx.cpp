/* build_seqidx - Builds the sequence index for the reference
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
   \file build_seqidx.cpp
   \brief build_seqidx.cpp Builds the sequence index for the reference.
   \author Massimiliano Rossi
   \date 07/08/2021
*/

#include <iostream>

#define VERBOSE


#include <malloc_count.h>
#include <kseq.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread);

#include <liftidx.hpp>
#include <common.hpp>

#include <filesystem>
namespace fs = std::filesystem;

inline lift::Lift get_null_lift(const size_t len)
{
    sdsl::bit_vector ibv(len);
    sdsl::bit_vector dbv(len);
    sdsl::bit_vector sbv(len);

    return lift::Lift(ibv, dbv, sbv);
}

//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
  std::string filename = "";
  std::string outpath = ""; // path where to output the file
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string usage("usage: " + std::string(argv[0]) + " infile [-o outpath]\n\n" +
                    "Computes the .idx file storing the sequence names and starting positions.\n" +
                    "outpath: [string]  - path to where to output the file.\n");

  std::string sarg;
  while ((c = getopt(argc, argv, "o:")) != -1)
  {
    switch (c)
    {
    case 'o':
      arg.outpath.assign(optarg);
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

  // Building the sequence idx

  verbose("Building the sequence index");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  gzFile fp(gzopen(args.filename.c_str(), "r"));
  if (fp == nullptr)
      error("gzopen() file " + std::string(args.filename) + " failed");

  kseq_t *seq = kseq_init(fp);

  std::vector<size_t> onset(1,0);
  std::vector<std::string> names;
  std::vector<std::pair<lift::Lift,size_t>> lifts;

  size_t u = 0;

  while (kseq_read(seq) >= 0)
  {
    lifts.push_back(std::make_pair(get_null_lift(seq->seq.l),u));

    u += seq->seq.l;
    names.push_back(std::string(seq->name.s));
    onset.push_back(u);
  }

  kseq_destroy(seq);
  gzclose(fp);

  liftidx idx(onset, names, lifts, u);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Sequence index construction complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());


  std::string outfile = "";
  if(args.outpath == "") outfile = args.filename;
  else outfile = args.outpath + fs::path(args.filename).filename().string();
  outfile += idx.get_file_extension();

  std::ofstream out(outfile);
  idx.serialize(out);

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Sequence index serialzation complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());
  return 0;
}