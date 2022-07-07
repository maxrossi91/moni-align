/* mems - MEMs types
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
   \file mems.hpp
   \brief mems.hpp MEMs type
   \author Massimiliano Rossi
   \date 19/10/2021
*/

#ifndef _MEMS_HH
#define _MEMS_HH

typedef long long int ll;

#define MATE_1 0 // The MEM cames from mate_1
#define MATE_2 1 // The MEM cames from mate_2
#define MATE_F 0 // The MEM cames from the forward strand
#define MATE_RC 2 // The MEM cames from the reverse-complement strand

typedef struct mem_t{
    size_t pos = 0;  // Position in the reference
    size_t len = 0;  // Length
    size_t idx = 0;  // Position in the pattern
    size_t mate = 0; // Left mate (0) or Right mate (1)
    size_t rpos = 0; // Position in the read for chaining
                        // With a Forward-Reverse library
                        // If the mem is in the FWD strand it is the position of the last character in the read
                        // If the mem is in the REV strand it is the position of the first character in the read
    std::vector<size_t> occs; // List of occurrences of the MEM

    mem_t(size_t p, size_t l, size_t i)
    {
        pos = p;  // Position in the reference
        len = l;  // Length of the MEM
        idx = i;  // Position in the read
    }
    mem_t(size_t p, size_t l, size_t i, size_t m, size_t r)
    {
        pos = p;  // Position in the reference
        len = l;  // Length of the MEM
        idx = i;  // Position in the read
        mate = m; // Left mate (0) or Right mate (1)
        rpos = r; // Position in the read for chaining
    }
} mem_t;



#endif /* end of include guard: _MEMS_HH */