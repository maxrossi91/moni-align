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

#include <moni.hpp>

#include <malloc_count.h>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>

int main(int argc, char *const argv[])
{
  using SelSd = SelectSdvec<>;
  using DagcSd = DirectAccessibleGammaCode<SelSd>;

  Args args;
  parseArgs(argc, argv, args);

  verbose("Loading the matching statistics index");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  
  ms_pointers<> ms;

  std::string filename_ms = args.filename + ms.get_file_extension();
  verbose("Index filename: " + filename_ms);

  ifstream fs_ms(filename_ms);
  ms.load(fs_ms);
  fs_ms.close();

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Matching statistics index construction complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Loading random access");
  t_insert_start = std::chrono::high_resolution_clock::now();

  std::string filename_slp = args.filename + ".slp";

  SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd> ra;

  ifstream fs(filename_slp);
  ra.load(fs);
  fs.close();

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Matching statistics index loading complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());


  sdsl::nullstream ns;

  size_t ms_size = ms.serialize(ns);
  // size_t ra_size = ra.serialize(ns);

  ms.print_stats();
  verbose("MS size (bytes): ", ms_size);
  // verbose("RA size (bytes): ", ra_size);

  verbose("Start generating the whole SA from back (PHI)");

  std::vector<size_t> sa_phi = ms.get_SA_Phi();

  verbose("Start generating the whole SA from front (PHI^-1)");

  std::vector<size_t> sa_phi_inv = ms.get_SA_Phi_inv();

  verbose("Checking if the two arrays are equal");

  std::reverse(sa_phi.begin(),sa_phi.end());
  if(sa_phi.size() != sa_phi_inv.size())
    error("Different sizes: " + std::to_string(sa_phi.size()) + 
                          " " + std::to_string(sa_phi_inv.size()));

  verbose("BWT length: " + std::to_string(sa_phi.size()));

  for(size_t i = 0; i < sa_phi.size(); ++i)
    if(sa_phi[i] != sa_phi_inv[i])
      error("Different in position " + std::to_string(i) + 
                        " values : " + std::to_string(sa_phi[i]) + 
                                " " + std::to_string(sa_phi_inv[i]));


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