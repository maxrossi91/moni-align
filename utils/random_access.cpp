//
//  random_access.cpp
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
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <Index_base_name> <start> <length>\n", argv[0]);
        return 1;
    }

    std::string filename = argv[1];
    size_t start = std::stol(argv[2]);
    size_t length = std::stol(argv[3]);

    plain_slp_t ra;
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();
    std::string filename_slp = filename + get_slp_file_extension<plain_slp_t>();
    verbose("Loading random access file: " + filename_slp);

    if (not file_exists(filename_slp))
        error("File not found: ", filename_slp);

    std::ifstream fs(filename_slp);
    ra.load(fs);
    fs.close();

    size_t n = ra.getLen();

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Random access loading complete");
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    liftidx idx;
    std::string filename_idx = filename + idx.get_file_extension();
    verbose("Loading fasta index file: " + filename_idx);

    if (not file_exists(filename_idx))
      error("File not found: ", filename_idx);

    std::ifstream fs_idx(filename_idx);
    idx.load(fs_idx);
    fs_idx.close();

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Fasta index loading complete");
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    char *ref = (char *)malloc(length + 1);
    ra.expandSubstr(start, length, ref);
    ref[length] = 0;

    auto origin = idx.index(start);
    auto lift = idx.index(idx.lift(start));

    verbose("String: ", ref );
    verbose("Origin: ", origin.first, ":", origin.second );
    verbose("  Lift: ", lift.first, ":", lift.second );

    delete ref;
    
    return 0;
}

