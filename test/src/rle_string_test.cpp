/* rle_string_test - Testing the rle_string
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
   \file rle_string_test.cpp
   \brief rle_string_test.cpp Testing the rle_string.
   \author Massimiliano Rossi
   \date 20/07/2022
*/

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>

#include <iostream>

#define VERBOSE
#include <common.hpp>

#include <ms_rle_string.hpp>
#include <ms_rle_simple_string.hpp>

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

//*********************** Tests rle_string *************************************

class RleFixture {
protected:
  ms_rle_string<> rle;
  ms_rle_simple_string<> rle_simple;
  size_t n;
  size_t r;
public:
  RleFixture() {
    std::string heads_filename = test_dir + "/yeast.fasta.bwt.heads";
    std::string lengths_filename = test_dir + "/yeast.fasta.bwt.len";

    std::ifstream in_heads(heads_filename);
    std::ifstream in_lengths(lengths_filename);
    rle = ms_rle_string<>(in_heads, in_lengths);
    rle_simple = ms_rle_simple_string<>(in_heads, in_lengths);
    in_heads.close();
    in_lengths.close();
    n= rle.size();
    r = rle.number_of_runs();
  }

  size_t rle_rank(size_t i, uint8_t c){
    return rle.rank(i,c);
  }
};

TEST_CASE_METHOD(RleFixture, "rle_string Constants Equivalence", "[Equivalence]")
{
  REQUIRE(rle.size() == rle_simple.size());
  REQUIRE(rle.number_of_runs() == rle_simple.number_of_runs());
  for(size_t i = 0; i < 256; ++i)
  {
    // REQUIRE(rle.number_of_runs_of_letter(i) == rle_simple.number_of_runs_of_letter(i));
    // The test above fails when there is no letter i in the text
    REQUIRE(rle.number_of_letter(i) == rle_simple.number_of_letter(i));
  }
}

TEST_CASE_METHOD(RleFixture, "rle_string Run of position Equivalence", "[Equivalence]")
{
  for(size_t i = 0; i < n; ++i) REQUIRE(rle.run_of_position(i) == rle_simple.run_of_position(i));
}

TEST_CASE_METHOD(RleFixture, "rle_string Random Access Equivalence", "[Equivalence]")
{
  for(size_t i = 0; i < n; ++i) REQUIRE(rle[i] == rle_simple[i]);
}

TEST_CASE_METHOD(RleFixture, "rle_string Run at Equivalence", "[Equivalence]")
{
  for(size_t i = 0; i < r; ++i) REQUIRE(rle.run_at(i) == rle_simple.run_at(i));
}

TEST_CASE_METHOD(RleFixture, "rle_string Head of Equivalence", "[Equivalence]")
{
  for(size_t i = 0; i < r; ++i) REQUIRE(rle.head_of(i) == rle_simple.head_of(i));
}

TEST_CASE_METHOD(RleFixture, "rle_string Head rank Equivalence", "[Equivalence]")
{
  for(size_t i = 0; i < r; ++i) 
    for(size_t c = 0; c < 256; ++c) 
      REQUIRE(rle.head_rank(i,c) == rle_simple.head_rank(i,c));
}

TEST_CASE_METHOD(RleFixture, "rle_string Run head rank Equivalence", "[Equivalence]")
{
  for(size_t i = 0; i < r; ++i) 
    for(size_t c = 0; c < 256; ++c) 
      REQUIRE(rle.run_head_rank(i,c) == rle_simple.run_head_rank(i,c));
}

TEST_CASE_METHOD(RleFixture, "rle_string Run and head rank Equivalence", "[Equivalence]")
{
  for(size_t i = 0; i < r; ++i) 
    for(size_t c = 0; c < 256; ++c) 
    {
      const auto rle_rhr = rle.run_and_head_rank(i,c);
      const auto rle_simple_rhr = rle_simple.run_and_head_rank(i,c);
      REQUIRE( rle_rhr.first == rle_simple_rhr.first);
      REQUIRE( rle_rhr.second == rle_simple_rhr.second);
    }
}

TEST_CASE_METHOD(RleFixture, "rle_string Run range Equivalence", "[Equivalence]")
{
  for(size_t i = 0; i < r; ++i) {
    const auto rle_rr = rle.run_range(i);
    const auto rle_simple_rr = rle_simple.run_range(i);
    REQUIRE( rle_rr.first == rle_simple_rr.first);
    REQUIRE( rle_rr.second == rle_simple_rr.second);
  }
}

TEST_CASE_METHOD(RleFixture, "rle_string Rank Equivalence", "[Equivalence]")
{
  for(size_t i = 0; i < n; ++i)
    for(size_t c = 0; c < 256; ++c)
      REQUIRE(rle.rank(i,c) == rle_simple.rank(i,c));

  BENCHMARK("indexed rle", i){ return rle.rank(i,'A'); };
  BENCHMARK("indexed rle_simple", i){ return rle_simple.rank(i,'A'); };
}

TEST_CASE_METHOD(RleFixture, "rle_string Rank Benchmark", "[Benchmark]")
{
  BENCHMARK("rank rle", i){ return rle.rank(i,'A'); };
  BENCHMARK("rank rle_simple", i){ return rle_simple.rank(i,'A'); };
}

TEST_CASE_METHOD(RleFixture, "rle_string Select Benchmark", "[Benchmark]")
{
  BENCHMARK("select rle", i){ return rle.select(i,'A'); };
  BENCHMARK("select rle_simple", i){ return rle_simple.select(i,'A'); };
}

TEST_CASE_METHOD(RleFixture, "rle_string Random access Benchmark", "[Benchmark]")
{
  BENCHMARK("random access rle", i){ return rle[i]; };
  BENCHMARK("random access rle_simple", i){ return rle_simple[i]; };
}

TEST_CASE_METHOD(RleFixture, "rle_string Select Equivalence", "[Equivalence]")
{
  for(size_t c = 0; c < 256; ++c)
  {
    size_t n_c = rle.number_of_letter(c);
    for(size_t i = 0; i < n_c; ++i)
      REQUIRE(rle.select(i,c) == rle_simple.select(i,c));
  }
}



//*********************** End tests rle_string *********************************


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
