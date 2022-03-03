/* r-index_test - Testing the r-index
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
   \file r-index_test.cpp
   \brief r-index_test.cpp Testing the r-index.
   \author Massimiliano Rossi
   \date 02/03/2022
*/

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>

#include <iostream>

#define VERBOSE
#include <common.hpp>

#include <moni.hpp>
#include <sam.hpp>
#include <mapq.hpp>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>
#include <PlainSlp.hpp>
#include <FixedBitLenCode.hpp>


#include <slp_definitions.hpp>
#include <chain.hpp>

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
typedef plain_slp_t slp_t;
typedef ms_pointers<> ms_t;
//*********************** End global variables *********************************

//*********************** Tests r-index ****************************************

TEST_CASE("r-index Phi", "[Phi]")
{
  ms_t ms;
  slp_t ra;
  // Load contigs lengths and structures  
  verbose("Loading r-index...");
  std::ifstream in_ms(test_dir + "/Chr21.100.thrbv.full.ms");
  ms.load(in_ms);
  in_ms.close();

  // Check r-index data   
  verbose("Loading slp data structure...");
  std::ifstream in_slp(test_dir +  "/Chr21.100.plain.slp");
  ra.load(in_slp);
  in_slp.close();

  size_t n = ra.getLen();

  size_t pos = 799258409;
  size_t len = 100;

  char *ref = (char *)malloc(len);

  ra.expandSubstr(pos, len, ref);
  std::cout << std::setw(10) << pos << ": " << std::string(ref) << std::endl;

  const auto sa_first = ms.get_first_run_sample();

  size_t count = 1;
  // Phi direction
  size_t curr = pos;
  if(curr != sa_first)
  {
    size_t next = ms.Phi(curr);
    if((n-curr) >= len and (n-next) >= len)
    {
      size_t lcp =  lceToRBounded(ra,curr,next,len);
      while(lcp >= len)
      {
        ra.expandSubstr(next, len, ref);
        std::cout << std::setw(10) << next << ": " << std::string(ref) << std::endl;

        if( ++count > 300)
          break;
        
        curr = next;
        if(curr != sa_first)
        {
          next = ms.Phi(curr);
          if((n-curr) >= len and (n-next) >= len)
            lcp  = lceToRBounded(ra,curr,next,len);
          else
            lcp = 0;
        }else lcp = 0;
      }
    }
  }
  
}

//*********************** End tests r-index ***********************************


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
