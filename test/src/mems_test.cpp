/* mems_test - Check number of MEMs generated
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
   \file mems_test.cpp
   \brief mems_test.cpp Check number of MEMs generated.
   \author Rahul Varki
   \date 13/6/2023
*/

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>

#include <iostream>

#define VERBOSE

#include <zlib.h>

#include <aligner_ksw2.hpp> //This file must be included first to avoid errors with common.hpp
#include <seed_finder.hpp>
#include <kseq.h>


//*********************** Global Variables *********************************
std::string test_dir = "../../../data/reads/";
std::string index_prefix = "../../../data/index/Chr21.10";
std::string filename_mate1 = "Chr21.10.HG002.R1.fastq.gz";
std::string filename_mate2 = "Chr21.10.HG002.R2.fastq.gz";
int unique_mems[10] = {37, 27, 24, 17, 18, 12, 36, 16, 36, 12};
//*********************** Testing chaining algorithm ************************

TEST_CASE("MEM Testing", "[mems]")
{   
    typedef struct paired_alignment_t
    {
        bool aligned = false;
        bool chained = false;
        bool best_score = false;
        bool second_best_score = false;

        kseq_t* mate1 = nullptr;
        kseq_t* mate2 = nullptr;

        kseq_t mate1_rev;
        kseq_t mate2_rev;

        sam_t sam_m1;
        sam_t sam_m2;

        // TODO: precompute the nt4 version of the mates

        int32_t min_score_m1 = 0;
        int32_t min_score_m2 = 0;
        int32_t min_score = 0;
        //paired_score_t score;
        int32_t score2 = 0;
        int32_t score2_m1 = 0;
        int32_t score2_m2 = 0;

        float frac_rep_m1 = 0.0;
        float frac_rep_m2 = 0.0;
        int sub_n = 0; // approximate number of suboptimal hits

        float mean = 0.0;
        float std_dev = 0.0;

        // Stats
        size_t n_seeds_dir1 = 0; // Number of seeds MATE1 MATE2_RC
        size_t n_seeds_dir2 = 0; // Number of seeds MATE2 MATE1_RC
        size_t n_mems_dir1 = 0; // Number of MEMs MATE1 MATE2_RC
        size_t n_mems_dir2 = 0; // Number of MEMs MATE2 MATE1_RC
        double avg_seed_length_dir1 = 0.0; // Average seed length MATE1 MATE2_RC
        double avg_seed_length_dir2 = 0.0;   // Average seed length  MATE2 MATE1_RC
        double avg_w_seed_length_dir1 = 0.0; // Average weighted seed length MATE1 MATE2_RC
        double avg_w_seed_length_dir2 = 0.0;   // Average weighted seed length  MATE2 MATE1_RC
        double armonic_avg_seed_length_dir1 = 0.0; // Average weighted seed length MATE1 MATE2_RC
        double armonic_avg_seed_length_dir2 = 0.0; // Average weighted seed length  MATE2 MATE1_RC

        std::vector<mem_t> mems;
        std::vector<std::pair<size_t, size_t>> anchors;
        std::vector<chain_t> chains;
        // (0) score (1) lift m1 (2) lift m2 (3) index in chains list (4) score m1 (5) score m2
        // std::vector<std::tuple<int32_t, size_t, size_t, size_t, size_t, size_t>> best_scores;
        //std::vector<paired_score_t> best_scores;
        
        paired_alignment_t(kseq_t *mate1_, kseq_t *mate2_, double mean_ = 0.0, double std_dev_ = 0.0) : 
        mate1(mate1_),
        mate2(mate2_),
        sam_m1(mate1),
        sam_m2(mate2),
        min_score_m1(20 + 8 * log(mate1->seq.l)),
        min_score_m2(20 + 8 * log(mate2->seq.l)),
        min_score(min_score_m1 + min_score_m2),
        mean(mean_),
        std_dev(std_dev_)
        {
        // TODO: Add parameter to decide whether the slash mate has to be removed
        remove_slash_mate(mate1);
        remove_slash_mate(mate2);

        rc_copy_kseq_t(&mate1_rev, mate1);
        rc_copy_kseq_t(&mate2_rev, mate2);

        // Initialize SAM infos
        // Fill sam fields RNEXT
        if(strcmp(mate1->name.s, mate2->name.s) == 0)
        {
            sam_m1.rnext = "=";
            sam_m2.rnext = "=";
        }
        else
        {
            sam_m1.rnext = std::string(mate2->name.s);
            sam_m2.rnext = std::string(mate1->name.s);
        }

        }

        paired_alignment_t() {}
        
        void init(kseq_t *mate1_, kseq_t *mate2_, double mean_ = 0.0, double std_dev_ = 0.0)
        {
            mate1 = mate1_;
            mate2 = mate2_;
            sam_m1 = sam_t(mate1);
            sam_m2 = sam_t(mate2);
            min_score_m1 = 20 + 8 * log(mate1->seq.l);
            min_score_m2 = 20 + 8 * log(mate2->seq.l);
            min_score = min_score_m1 + min_score_m2;
            mean = mean_;
            std_dev = std_dev_;
            // TODO: Add parameter to decide whether the slash mate has to be removed
            remove_slash_mate(mate1);
            remove_slash_mate(mate2);

            rc_copy_kseq_t(&mate1_rev, mate1);
            rc_copy_kseq_t(&mate2_rev, mate2);

            // Initialize SAM infos
            // Fill sam fields RNEXT
            if(strcmp(mate1->name.s, mate2->name.s) == 0)
            {
                sam_m1.rnext = "=";
                sam_m2.rnext = "=";
            }
            else
            {
                sam_m1.rnext = std::string(mate2->name.s);
                sam_m2.rnext = std::string(mate1->name.s);
            }

        }

        void write(FILE* out)
        {
            write_sam(out, sam_m1);
            write_sam(out, sam_m2);
        }

        void set_sam_not_aligned()
        {
            sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;
        }

        void release_memory()
        {
            free_kseq_t(&mate1_rev);
            free_kseq_t(&mate2_rev);
        }

        ~paired_alignment_t()
        {
            mate1 = nullptr;
            mate2 = nullptr;
            release_memory();
        }

    } paired_alignment_t;

    //*************************** Start Testing **************************************
    kseq_t *mate1 = nullptr;
    kseq_t *mate2 = nullptr;

    verbose("Attempting to open ", test_dir + filename_mate1);
    gzFile fp_mate1 = gzopen((test_dir + filename_mate1).c_str(), "r");
    if (!fp_mate1)
    {
        verbose("Failed to open ", test_dir + filename_mate1);
    }
    verbose("Attempting to open ", test_dir + filename_mate2);
    gzFile fp_mate2 = gzopen((test_dir + filename_mate2).c_str(), "r");
    if (!fp_mate2)
    {
        verbose("Failed to open ", test_dir + filename_mate1);
    }
    
    //Initialize the reads stored in the fastq files
    mate1 = kseq_init(fp_mate1);
    mate2 = kseq_init(fp_mate2);
    kpbseq_t *b = kpbseq_init();
    size_t b_size = 10; //Using 2 results in seg fault
    size_t l = 0;
    l = kpbseq_read(b, mate1, mate2, b_size);
    REQUIRE(l == b_size);
   
    std::vector<std::vector<paired_alignment_t>> alignments;
    std::vector<kpbseq_t *> memo;

    memo.push_back(kpbseq_init());
    copy_kpbseq_t(memo.back(), b);
    alignments.push_back(std::vector<paired_alignment_t>(l));
    kpbseq_t *batch = memo.back();

    verbose("Initializing the MEM finder object");
    seed_finder<plain_slp_t, ms_pointers<>> mem_finder = seed_finder<plain_slp_t, ms_pointers<>>(index_prefix, 25, false, 5000); //Have to provide the types for template class and functions

    int num_reads = batch->mate1->l;
    for (size_t i = 0; i < num_reads; ++i)
    {
        paired_alignment_t& alignment = alignments.back()[i];
        alignment.init(&batch->mate1->buf[i], &batch->mate2->buf[i]);
        verbose("Processing read: ", std::string(alignment.mate1->name.s));
        mem_finder.find_seeds(alignment.mate1, alignment.mems, 0, MATE_1 | MATE_F);
        mem_finder.find_seeds(&alignment.mate1_rev, alignment.mems, alignment.mate2->seq.l, MATE_1 | MATE_RC);
        mem_finder.find_seeds(alignment.mate2, alignment.mems, 0, MATE_2 | MATE_F);
        mem_finder.find_seeds(&alignment.mate2_rev, alignment.mems, alignment.mate1->seq.l, MATE_2 | MATE_RC);
        verbose("Found this many unique MEMs: ", alignment.mems.size());
        REQUIRE(alignment.mems.size() == unique_mems[i]);
    }

    //Call the destructors 
    kpbseq_destroy(memo.back());
    kseq_destroy(mate1);
    kseq_destroy(mate2);
    gzclose(fp_mate1);
    gzclose(fp_mate2);

    bool finish = true;
    REQUIRE(finish);
}

//************************* End Testing *************************************

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