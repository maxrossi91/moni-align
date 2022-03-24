//
//  extract.cpp
//
//  Copyright 2020 Massimiliano Rossi. All rights reserved.
//

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


int main(int argc, char* argv[])
{
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <Index_base_name> <seq_name>\n", argv[0]);
        return 1;
    }

    std::string filename = argv[1];
    std::string chr = argv[2];

    plain_slp_t ra;
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();
    std::string filename_slp = filename + get_slp_file_extension<plain_slp_t>();
    // verbose("Loading random access file: " + filename_slp);

    if (not file_exists(filename_slp))
        error("File not found: ", filename_slp);

    std::ifstream fs(filename_slp);
    ra.load(fs);
    fs.close();

    size_t n = ra.getLen();

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    // verbose("Random access loading complete");
    // verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    liftidx idx;
    std::string filename_idx = filename + idx.get_file_extension();
    // verbose("Loading fasta index file: " + filename_idx);

    if (not file_exists(filename_idx))
      error("File not found: ", filename_idx);

    std::ifstream fs_idx(filename_idx);
    idx.load(fs_idx);
    fs_idx.close();

    t_insert_end = std::chrono::high_resolution_clock::now();

    // verbose("Fasta index loading complete");
    // verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    auto& names = idx.get_names();
    size_t i = 0;
    size_t j = 0;

    while (i < n and names[j].compare(chr))
    {
        i += idx.length(j);
        ++j;
    }

    if (i >= n)
        error("Sequence ", chr, " not found!");

    std::cout << ">" << chr << std::endl;

    size_t chr_len = idx.length(j);

    std::cerr << "Extracting " << chr << " of length " << chr_len << std::endl;

    size_t buff_size = 600;
    char *ref = (char *)malloc(buff_size + 1);

    size_t k = 0;
    while ( k < chr_len )
    {
        size_t len = std::min(chr_len - k, buff_size);
        ra.expandSubstr(i + k, len, ref);
        ref[len] = 0;
        j = 0;
        while ( j < len )
        {
            size_t l = std::min(len - j, (size_t)60);
            printf("%.*s\n", l, (ref + j) );
            j+= l;
        }
        k += j;
    }

    delete ref;
    
    return 0;
}