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

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>

#include <iostream>

#define VERBOSE

#include <kseq.h>
#include <zlib.h>

KSEQ_DECLARE(gzFile);
#include <liftidx.hpp>


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

std::string test_dir = "../data";

//*********************** End global variables *********************************

//*********************** Tests LevioSAM ***************************************

TEST_CASE("LevioSAM Chr19 pangenome", "[LevioSAM]")
{

  // Load contigs lengths and structures  
  verbose("Loading contigs lengths and names...");
  std::ifstream in_lidx(test_dir + "/Chr21.10.lidx");
  std::vector<std::string> from_lidx_file_names;
  std::vector<size_t> from_lidx_file_lengths;
  while (not in_lidx.eof()) 
  { 
      std::string tmp_name;
      std::size_t tmp_length;
      in_lidx >> tmp_name >> tmp_length;
      if (tmp_name != "")
      {
          from_lidx_file_names.push_back(tmp_name);
          from_lidx_file_lengths.push_back(tmp_length);
      }
  }


  // Check lifting data   
  verbose("Loading ldx data structure...");
  std::ifstream in_ldx(test_dir + "/Chr21.10.ldx");
  liftidx ldx;
  ldx.load(in_ldx);
  in_ldx.close();

  // Loding golden truth lifting
  verbose("Loading LevioSAM indicies...");
  std::vector< std::string > filenames = {
    "HG00096_H1_21.lft",
    "HG00096_H2_21.lft",
    "HG00097_H1_21.lft",
    "HG00097_H2_21.lft",
    "HG00099_H1_21.lft",
    "HG00099_H2_21.lft",
    "HG00100_H1_21.lft",
    "HG00100_H2_21.lft"
  };
  // those are necessary because the construction of LevioSAM is buggy
  std::vector< size_t > limits = {
    46708362,
    46708362,
    46708095,
    46708070,
    46708362,
    46708362,
    46708362,
    46708362
  };

  std::vector< lift::LiftMap > lifts(filenames.size());

  for ( size_t i = 0; i < filenames.size(); ++i )
  {
    verbose("Loading ", test_dir + "/lifts/" + filenames[i]);
    std::ifstream in_lft(test_dir + "/lifts/" + filenames[i]);
    if ( not in_lft.is_open() ) error("Error opening file...");
    lifts[i].load(in_lft);
    in_lft.close();
  }

  verbose("Checking correctness of ldx data structure");

  bool check = true;
  size_t i = 0;
  for( i = 0; i < from_lidx_file_lengths[0] and check; ++i)
      check = check and (ldx.lift(i) == i);
  verbose("First missmatch: ", i);
  REQUIRE(check);

  for (size_t j = 1; j < from_lidx_file_lengths.size(); ++j)
  {
    const std::string& contig_name = from_lidx_file_names[0];
    size_t k = 0 ;
    auto s1_lengths = lifts[j-1].get_s1_lens();
    auto length_map = lifts[j-1].get_lengthmap();
    for( k = 0 ; k < from_lidx_file_lengths[j] and check; ++i, ++k)
      if( k < limits[j-1])
        check = check and ((ldx.lift(i)) == lifts[j-1].lift_pos(contig_name,k));
      else
        check = check and ((ldx.lift(i)) == ((k - limits[j-1] + 1) + lifts[j-1].lift_pos(contig_name,limits[j-1]-1)));
    verbose("First missmatch: ", i);
    if ( not check )
      verbose("Left: ", ldx.lift(i-1) , " Right: ", lifts[j-1].lift_pos(contig_name,k-1), " Pos: ", ldx.index(i-1).first, " ", ldx.index(i-1).second);

    REQUIRE(check);
  }
}

//*********************** End tests LevioSAM ***********************************


int main( int argc, char* argv[] )
{
  Catch::Session session; // There must be exactly one instance
  
  // Build a new parser on top of Catch2's
  using namespace Catch::Clara;
  auto cli = session.cli() |
  Opt( test_dir, "dir" ) ["--test-dir"] ("test directory");

  // Now pass the new composite back to Catch2 so it uses that
  session.cli( cli );

  // Let Catch2 (using Clara) parse the command line
  int returnCode = session.applyCommandLine( argc, argv );
  if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

  session.run();
}
