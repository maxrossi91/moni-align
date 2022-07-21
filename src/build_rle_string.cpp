/* build_rle_string - Builds the rle_strings
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
   \file build_rle_string.cpp
   \brief build_rle_string.cpp Builds the rle_strings.
   \author Massimiliano Rossi
   \date 20/07/2022
*/

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <malloc_count.h>

#include <ms_rle_string.hpp>
#include <ms_rle_simple_string.hpp>

//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
  std::string input = "";
  std::string outpath = ""; // path where to output the file
  bool simple = false; // build the simple rle_string
  bool rle = false; // build the simple rle_string
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string usage("usage: " + std::string(argv[0]) + " infile [-o outpath] [-s simple] [-r rle]\n\n" +
                    "Computes the rle string and dumps it to file.\n" +
                    " simple: [bool]    - build simple rle_string.\n" +
                    "    rle: [bool]    - consider run-length encoded pair of files as input.\n" +
                    "outpath: [string]  - path to where to output the file.\n");

  std::string sarg;
  while ((c = getopt(argc, argv, "o:sr")) != -1)
  {
    switch (c)
    {
    case 'o':
      arg.outpath.assign(optarg);
      break;
    case 's':
      arg.simple = true;
      break;
    case 'r':
      arg.rle = true;
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
    arg.input.assign(argv[optind]);
  }
  else
  {
    error("Invalid number of arguments\n", usage);
  }
}

//********** end argument options ********************

template<typename rle_t>
void build(const Args& args)
{
  // Building the sequence idx


  if(args.rle)
  {
    verbose("Building the rle string");
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    std::string heads_filename = args.input + ".heads";
    std::string lengths_filename = args.input + ".len";

    std::ifstream in_heads(heads_filename);
    std::ifstream in_lengths(lengths_filename);

    rle_t rle(in_heads, in_lengths);

    in_heads.close();
    in_lengths.close();

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Rle string construction complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    t_insert_start = std::chrono::high_resolution_clock::now();

    std::string out_filename = args.input + rle.get_file_extension();
    verbose("Writing the rle string to file: ", out_filename);

    std::ofstream out(out_filename);
    rle.serialize(out);

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Rle serialization complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  }else{
    
    verbose("Building the rle string");
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    std::ifstream in(args.input);
    rle_t rle(in);
    in.close();

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Rle string construction complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    
    t_insert_start = std::chrono::high_resolution_clock::now();

    std::string out_filename = args.input + rle.get_file_extension();
    verbose("Writing the rle string to file: ", out_filename);

    std::ofstream out(out_filename);
    rle.serialize(out);

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Rle serialization complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  }

}


int main(int argc, char *const argv[])
{
  Args args;
  parseArgs(argc, argv, args);

  if(args.simple) build<ms_rle_simple_string<>>(args);
  else build<ms_rle_string<>>(args);

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());
  return 0;
}