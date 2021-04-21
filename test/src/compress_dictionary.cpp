/* compress_dictionary - Computes the compressed dictionary from prefix-free parse dictionary
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
   \file compress_dictionary.cpp
   \brief compress_dictionary.cpp Computes the compressed dictionary from prefix-free parse dictionary.
   \author Massimiliano Rossi
   \date 16/09/2020
*/

#include <iostream>

#define VERBOSE

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

  // Building the r-index

  verbose("Compressing the dictionary");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  // Open output files
  std::string dicz_filename = args.filename + ".dicz";
  std::string dicz_len_filename = args.filename + ".dicz.len";

  FILE *dicz;
  FILE *dicz_len;

  if ((dicz = fopen(dicz_filename.c_str(), "w")) == nullptr)
    error("open() file " + std::string(dicz_filename) + " failed");

  if ((dicz_len = fopen(dicz_len_filename.c_str(), "w")) == nullptr)
    error("open() file " + std::string(dicz_len_filename) + " failed");

  // Open the dictionary
  std::string dict_filename = args.filename + ".dict";
  std::vector<uint8_t> dict;
  read_file(dict_filename.c_str(), dict);

  // Start processing

  
  // Generating phrase lengths
  verbose("Generating phrase lengths");
  std::vector<size_t> lengths(1,0);
  
  // Counting the number of Dollars at the beginning
  size_t i = 0, j = 0;
  while(dict[i++] == Dollar)
    j++;
  dict.erase(dict.begin(), dict.begin() + j);

  for(auto chr: dict)
  {
    // Skip the Dollars
    if(chr == EndOfDict)
      continue;

    // Hit end of phrase
    if(chr == EndOfWord)
      lengths.push_back(0);
    else
      lengths.back()++;
  }

  if (lengths.back()==0)
    lengths.pop_back();

  verbose("Found", lengths.size(), " phrases ");

  verbose("Generating phrases");
  uint8_t* ptr = dict.data(); // Beginning of the current phrase
  for(auto length: lengths)
  {
    size_t compressed_length = length - args.w;

    if ((fwrite(&compressed_length, 4, 1, dicz_len)) != 1)
      error("fwrite() file " + std::string(dicz_len_filename) + " failed");

    if ((fwrite(ptr, sizeof(uint8_t), compressed_length, dicz)) != compressed_length)
      error("fwrite() file " + std::string(dicz_filename) + " failed");

    ptr += length + 1;
  }


  fclose(dicz);
  fclose(dicz_len);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

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