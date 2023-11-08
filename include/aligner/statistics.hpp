/* statistics - statistics type
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
   \file statistics.hpp
   \brief statistics.hpp STATISTICs type
   \author Massimiliano Rossi
   \date 16/07/2022
*/

#ifndef _STATISTICS_HH
#define _STATISTICS_HH

#include <common.hpp>

typedef struct statistics_t{
    size_t processed_reads = 0;
    size_t aligned_reads = 0;
    size_t orphan_reads = 0;
    size_t orphan_recovered_reads = 0;

    statistics_t& operator+=(const statistics_t& other)
    {
        processed_reads += other.processed_reads;
        aligned_reads += other.aligned_reads;
        orphan_reads += other.orphan_reads;
        orphan_recovered_reads += other.orphan_recovered_reads;
        return *this;
    }

    inline std::string to_string() const
    {
        return  "\n\t       Processed reads: " + std::to_string(processed_reads) +
                "\n\t         Aligned reads: " + std::to_string(aligned_reads) +
                "\n\t          Orphan reads: " + std::to_string(orphan_reads) +
                "\n\tOrphan recovered reads: " + std::to_string(orphan_recovered_reads);
    }

    inline void print() const
    {
        verbose("   Alignment statistics ");
        verbose("       Processed reads: ", processed_reads);
        verbose("         Aligned reads: ", aligned_reads);
        verbose("          Orphan reads: ", orphan_reads);
        verbose("Orphan recovered reads: ", orphan_recovered_reads);
    }

} statistics_t;



#endif /* end of include guard: _STATISTICS_HH */