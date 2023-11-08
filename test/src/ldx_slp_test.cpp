/* ldx_slp_test - Testing lifting and random access
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
   \file ldx_slp_test.cpp
   \brief ldx_slp_test.cpp Testing lifting and random access.
   \author Massimiliano Rossi
   \date 19/11/2021
*/

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>

#include <iostream>
#include <fstream>

#define VERBOSE
#include <leviosam.hpp> // Included here to avoid conflicts with common.hpp
#include <common.hpp>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>
#include <PlainSlp.hpp>
#include <FixedBitLenCode.hpp>

#include <kseq.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)
#include <liftidx.hpp>

#include <slp_definitions.hpp>

// //*********************** Catch2 listener update *******************************

// struct listener : Catch::TestEventListenerBase
// {
//     using TestEventListenerBase::TestEventListenerBase;
    
//     virtual void testCaseStarting(Catch::TestCaseInfo const& testInfo) override
//     {
//         std::cout << testInfo.tagsAsString() << " " << testInfo.name << std::endl;
//     }
// };
// CATCH_REGISTER_LISTENER(listener)
// //*********************** End Catch2 listener update ***************************

//*********************** Global variables *************************************

std::string base_name = "../data/Chr21.10";
std::size_t w = 10;

//*********************** End global variables *********************************

//*********************** Tests LevioSAM ***************************************

TEST_CASE("ldx and slp data structure test", "[ldx_slp]")
{

  // Load lifting index   
  verbose("Loading ldx data structure...");
  liftidx ldx;
  std::string filename_ldx = base_name + ldx.get_file_extension();

  std::ifstream in_ldx(filename_ldx);
  ldx.load(in_ldx);
  in_ldx.close();

  // Loding golden truth lifting
  verbose("Loading random access SLP...");

  plain_slp_t ra;
  std::string filename_slp = base_name + get_slp_file_extension<plain_slp_t>();

  std::ifstream fs(filename_slp);
  ra.load(fs);
  fs.close();


  verbose("Checking correctness of ldx data structure");

  bool check = true;

  size_t n = ra.getLen();
  auto& names = ldx.get_names();

  size_t buff_size = 600;
  char *ref = (char *)malloc(buff_size + 1);


  size_t acc = 0;
  for( size_t i = 0; i < names.size(); ++i)
  {
    verbose("Processing ", names[i]);
    size_t last = w + ((i == (names.size() -1)) ? w-1 : 0);
    size_t chr_len = ldx.length(i) - last;
    bool check = true;
    size_t k = 0;
    while (k < chr_len)
    {
      size_t len = std::min(chr_len - k, buff_size);
      ra.expandSubstr(acc + k, len, ref);
      ref[len] = 0;
      size_t j = 0;
      for (; j < len and check; ++j)
        check = check and (ref[j] > 5);
      if( not check )
        verbose("First character <=5 in ", names[i], " at position: ", k + j - 1);
      k += len;
    }
    REQUIRE(check);

    check = true;
    
    ra.expandSubstr(acc + k, last, ref);
    ref[last] = 0;
    size_t j = 0;
    for (; j < last and check; ++j)
      check = check and (ref[j] <= 5);
    if (not check)
      verbose("First character >5 in ", names[i], " at position: ", k + j - 1, j -1);
    k += last;

    REQUIRE(check);
    acc += k;
    
  }

  delete ref;
}

//*********************** End tests LevioSAM ***********************************


int main( int argc, char* argv[] )
{
  Catch::Session session; // There must be exactly one instance
  
  // Build a new parser on top of Catch2's
  using namespace Catch::Clara;
  auto cli = session.cli() |
             Opt(base_name, "index")["--base_name"]("index base name") |
             Opt(w, "int")["-W"]("specify w");

  // Now pass the new composite back to Catch2 so it uses that
  session.cli( cli );

  // Let Catch2 (using Clara) parse the command line
  int returnCode = session.applyCommandLine( argc, argv );
  if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

  session.run();
}
