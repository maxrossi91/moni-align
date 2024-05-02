/* csv - CSV format writer.
    Copyright (C) 2024 Rahul Varki
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
   \file csv.hpp
   \brief csv.hpp CSV format writer. Debugging/Logging purposes only.
   \author Rahul Varki
   \date 5/1/2024
*/

#ifndef _CSV_HH
#define _CSV_HH

#include <common.hpp>
#include <kpbseq.h>

typedef struct csv_t{
    const kseq_t *read = nullptr; // The read of the SAM entry. Contains: QNAME, SEQ, and QUAL
    size_t num_uniq_mems = 0; // The number of unique MEMs between the read and reference
    size_t total_mem_occ = 0; // The number of total occurrances of all the MEMs between the read and reference
    double max_mem_freq = 0;  // The maximum frequency of a particular MEM in the pangenome. 
    double min_mem_freq = 1; // The minimum frequency of a particular MEM in the pangenome. 
    size_t high_occ_mem = 0; // The occurance count of the highest occurance MEM on any one particular genome
    size_t low_occ_mem = 0;  // The occurance count of the lowest occurance MEM on any one particular genome
    size_t num_mems_filter = 0; // The number of MEMs filtered overall
    size_t num_chains_skipped = 0; // The number of chains skipped due to left MEM position similiarity and matching chain score

    // Constructor when passed read
    csv_t(const kseq_t* read_)
    {
        read = read_;
    }

    // Empty Constructor
    csv_t(){}

    // Destructor
    ~csv_t()
    {
        read = nullptr;
    }
} csv_t;

inline void write_csv(FILE *out, const csv_t c)
{
    fprintf(out, "%s,", c.read->name.s);
    fprintf(out, "%zu,", c.num_uniq_mems);
    fprintf(out, "%zu,", c.total_mem_occ);
    fprintf(out, "%f,", c.max_mem_freq);
    fprintf(out, "%f,", c.min_mem_freq);
    fprintf(out, "%zu,", c.high_occ_mem);
    fprintf(out, "%zu,", c.low_occ_mem);
    fprintf(out, "%zu,", c.num_mems_filter);
    fprintf(out, "%zu", c.num_chains_skipped);
    fprintf(out, "\n");
}


#endif /* end of include guard: _CSV_HH */