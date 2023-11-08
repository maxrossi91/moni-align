//
//  blast-like-visualizer.cpp
//
//  Copyright 2020 Massimiliano Rossi. All rights reserved.
//

#include <iostream>
#include <fstream>
#define VERBOSE
#include <common.hpp>
#include <ksw2.h>
#include <ksw.h>


static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b);

int main(int argc, char* argv[])
{
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <in.seq> <out.seq>\n", argv[0]);
        return 1;
    }
    unsigned char seq_nt4_table[256] = {
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };

    unsigned char nt4_seq_table[4] = {'A', 'C', 'G', 'T'};

    std::string a = argv[1];
    std::string b = argv[2];
    verbose("Aligning : ", a );
    verbose("  Against: ", b );
    
    verbose("Reading sequences...");

    for(size_t i = 0; i < a.size(); ++i)
        a[i] = seq_nt4_table[(int)a[i]];
    for(size_t i = 0; i < b.size(); ++i)
        b[i] = seq_nt4_table[(int)b[i]];
    
    const int8_t smatch = 2;      // Match score default
    const int8_t smismatch = 4;   // Mismatch score default
    const int8_t gapo = 4;        // Gap open penalty
    const int8_t gapo2 = 13;      // Gap open penalty
    const int8_t gape = 2;        // Gap extension penalty
    const int8_t gape2 = 1;       // Gap extension penalty
    const int end_bonus = 400;    // Bonus to add at the extension score to declare the alignment
    
    const int w = -1;             // Band width
    const int zdrop = -1;

    void *km = 0;           // Kalloc

    const int m = 5;
    int8_t mat[25];

    int flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT;

    ksw_gen_simple_mat(m,mat,smatch,-smismatch);

    { // KSW
        std::cout << "ksw_align" << std::endl;
        const int xtra = KSW_XSTART;
        kswq_t *q = 0;
        kswr_t r = ksw_align(a.size(), (uint8_t*)a.data(), b.size(), (uint8_t*)b.data(), m, mat, gapo, gape, xtra, &q);
        
        // std::string blc = print_BLAST_like((uint8_t*)b.data(),(uint8_t*)a.data(),ez.cigar,ez.n_cigar);
        // for(size_t i = 0; i < blc.size(); ++i)
        // {   
        //     int d = (int)blc[i] - '0';
        //     if( d >= 0 and d <= 3)
        //         blc[i] = nt4_seq_table[d];
        // }
        
        // std::cout<<std::endl<<blc;
        std::cout << "tb: " << r.tb<< " te: " << r.te << std::endl<< std::endl;
        std::cout << "Score: " << r.score << std::endl<< std::endl;
    }
    

    { // KSW
        std::cout << "ksw_global" << std::endl;
        ksw_extz_t ez;
        memset(&ez, 0, sizeof(ksw_extz_t));

        ez.score = ksw_global(a.size(), (uint8_t*)a.data(), b.size(), (uint8_t*)b.data(), m, mat, gapo, gape, a.size(), &ez.n_cigar, &ez.cigar);
        
        std::string blc = print_BLAST_like((uint8_t*)b.data(),(uint8_t*)a.data(),ez.cigar,ez.n_cigar);
        for(size_t i = 0; i < blc.size(); ++i)
        {   
            int d = (int)blc[i] - '0';
            if( d >= 0 and d <= 3)
                blc[i] = nt4_seq_table[d];
        }
        
        std::cout<<std::endl<<blc;
        std::cout << "Score: " << ez.score << std::endl<< std::endl;
    }
    
    { // KSW2
        std::cout << "ksw_extz2_sse" << std::endl;
        ksw_extz_t ez;
        memset(&ez, 0, sizeof(ksw_extz_t));

        ksw_extz2_sse(km, a.size(), (uint8_t*)a.data(), b.size(), (uint8_t*)b.data(), m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez);

        std::string blc = print_BLAST_like((uint8_t*)b.data(),(uint8_t*)a.data(),ez.cigar,ez.n_cigar);
        for(size_t i = 0; i < blc.size(); ++i)
        {   
            int d = (int)blc[i] - '0';
            if( d >= 0 and d <= 3)
                blc[i] = nt4_seq_table[d];
        }
        
        std::cout<<std::endl<<blc;
        std::cout << "Score: " << ez.score << std::endl<< std::endl;
    }


    return 0;
}


static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
    int i, j;
    a = a < 0? -a : a;
    b = b > 0? -b : b;
    for (i = 0; i < m - 1; ++i) {
        for (j = 0; j < m - 1; ++j)
        mat[i * m + j] = i == j? a : b;
        mat[i * m + m - 1] = 0;
    }
    for (j = 0; j < m; ++j)
        mat[(m - 1) * m + j] = 0;
}