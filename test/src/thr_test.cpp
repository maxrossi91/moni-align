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