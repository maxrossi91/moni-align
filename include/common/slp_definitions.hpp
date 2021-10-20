/* slp_definitions - Definitions wrapper for SLP.
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
   \file slp_definitions.hpp
   \brief slp_definitions.hpp Definitions wrapper for SLP.
   \author Massimiliano Rossi
   \date 19/10/2021
*/

#ifndef _SLP_DEFINITIONS_HH
#define _SLP_DEFINITIONS_HH


////////////////////////////////////////////////////////////////////////////////
/// SLP definitions
////////////////////////////////////////////////////////////////////////////////

using SelSd = SelectSdvec<>;
using DagcSd = DirectAccessibleGammaCode<SelSd>;
using Fblc = FixedBitLenCode<>;

using shaped_slp_t = SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd>;
using plain_slp_t = PlainSlp<uint32_t, Fblc, Fblc>;

template< typename slp_t>
std::string get_slp_file_extension()
{
  return std::string(".slp");
}

template <>
std::string get_slp_file_extension<shaped_slp_t>()
{
  return std::string(".slp");
}

template <>
std::string get_slp_file_extension<plain_slp_t>()
{
  return std::string(".plain.slp");
}
////////////////////////////////////////////////////////////////////////////////


#endif /* end of include guard: _SLP_DEFINITIONS_HH */