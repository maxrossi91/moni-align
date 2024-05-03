/* aligner_ksw2 - Align the reads to the reference using the ksw2 library for SW
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
   \file aligner_ksw2.cpp
   \brief aligner_ksw2.cpp Align the reads to the reference using the ksw2 library for SW
   \author Massimiliano Rossi
   \date 13/07/2020
*/

#ifndef _ALIGNER_KSW2_HH
#define _ALIGNER_KSW2_HH

#define VERBOSE
#define MTIME

#include <vector>
#include <leviosam.hpp> // Included here to avoid conflicts with common.hpp
#include <common.hpp>

#include <sdsl/io.hpp>

#include <moni.hpp>
#include <sam.hpp>
#include <csv.hpp>
#include <mapq.hpp>

#include <malloc_count.h>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>
#include <PlainSlp.hpp>
#include <FixedBitLenCode.hpp>

#include <ksw.h>
#include <ksw2.h>

#include <libgen.h>
#include <kpbseq.h>
#include <liftidx.hpp>

MTIME_TSAFE_INIT(10);
MTIME_TSAFE_NAMES_INIT(
  "Seed extraction", 
  "Chaining", 
  "Chain alignment", 
  "Chaining - Anchor extraction", 
  "Chaining - Dynamic Programming", 
  "Chaining - Chain position", 
  "Chaining - Backtracking", 
  "Orphan recovery",
  "Seed extraction - Seeds finding",
  "Seed extraction - Seeds occurrences"
);
#include <slp_definitions.hpp>
#include <chain.hpp>
#include <statistics.hpp>

#include <htslib/sam.h>

#define _REALIGN

template <typename seed_finder_t>
class aligner
{
public:
    // using SelSd = SelectSdvec<>;
    // using DagcSd = DirectAccessibleGammaCode<SelSd>;
    // using SlpT = SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd>;
    using ll = long long int;


    typedef struct{

        size_t min_len = 25;    // Minimum MEM length
        size_t ext_len = 100;   // Extension length
        size_t top_k = 1;       // Report the top_k alignments
        size_t check_k = 5;     // Check the scores of the check_k chains with different scores alignments
        size_t region_dist = 10; // Maximum distance for two regions to be called "from the same region"
        // double mean = 5;         // The mean of the insert size distribution.
        // double std_dev = 10;     // The standard deviation of the insert size distribution.

        bool filter_dir = true; // Use MEMs average length to filter the orientation of the reads
        double dir_thr = 50.0; // Use MEMs average length distance to filter the orientation of the reads

        bool filter_seeds = true; // Filter seed if occurs more than threshold
        size_t n_seeds_thr = 5000;   // Filter seed if occurs more than threshold

        bool filter_freq = true;  // Filter seed if it occurs with frequency greater than threshold
        double freq_thr = 0.02;   // Filter seed if it occurs with frequency greater than threshold
        
        // ksw2 parameters
        int8_t smatch = 2;      // Match score default
        int8_t smismatch = 4;   // Mismatch score default
        int8_t gapo = 4;        // Gap open penalty
        int8_t gapo2 = 13;      // Gap open penalty
        int8_t gape = 2;        // Gap extension penalty
        int8_t gape2 = 1;       // Gap extension penalty
        int end_bonus = 400;    // Bonus to add at the extension score to declare the alignment

        int w = -1;             // Band width
        int zdrop = -1;         // Zdrop enable

        bool forward_only = true;      // Align only
        bool report_mems = false;      // Report MEMs instead of read alignment
        bool csv = false;             // Report MEM statistics in CSV file

        // Chaining parameters
        ll max_dist_x = 500;    // Max distance for two anchors to be chained
        ll max_dist_y = 100;    // Max distance for two anchors from the same read to be chained
        ll max_iter = 50;       // Max number of iterations of the chaining algorithhm
        ll max_pred = 50;       // Max number of predecessor to be considered
        ll min_chain_score = 40;// Minimum chain score
        ll min_chain_length = 1;// Minimum chain length
        bool secondary_chains = false; // Find secondary chains in paired-end setting

    } config_t;



    typedef struct score_t{
      int32_t score = 0;
      size_t pos = 0; // Position of the leftmost match of the read in the global reference
      size_t lft = 0; // Position of the lift in the reference coordinates
      bool unmapped_lft = false; // Set if the liftover results in unmapped read
    } score_t;

    /**
    * @brief Store the read, its reverse-complement, and the set of mems.
    * 
    */
    typedef struct alignment_t{
      bool aligned = false;
      bool chained = false;
      bool best_score = false;

      kseq_t* read;
      kseq_t read_rev;

      sam_t sam;
      csv_t csv;
      // bam1_t bam;

      const size_t min_score;
      score_t score;
      int32_t score2 = 0;

      float frac_rep = 0.0;
      int sub_n = 0; // approximate number of suboptimal hits

      std::vector<mem_t> mems;
      std::vector<std::pair<size_t, size_t>> anchors;
      std::vector<chain_t> chains;


      alignment_t(kseq_t *read_):
          read(read_),
          sam(read),
          csv(read),
          min_score(20 + 8 * log(read->seq.l))
      {
        rc_copy_kseq_t(&read_rev, read);
      }

      // void write(samFile* out, const sam_hdr_t *hdr)
      // {
      //   sam_write1(out, hdr, &bam);
      // }

      void write(FILE* out)
      {
        write_sam(out, sam);
      }

      void record_csv(FILE* out)
      {
        write_csv(out, csv);
      }

      void set_sam_not_aligned()
      {
        sam.flag = SAM_UNMAPPED;
      }

      void release_memory()
      {
        free_kseq_t(&read_rev);
      }
      
      ~alignment_t()
      {
        read = nullptr;
        release_memory();
      }

    } alignment_t;


    aligner(std::string filename, 
            config_t config = config_t()) : 
                min_len(config.min_len),        // Minimum MEM length
                ext_len(config.ext_len),        // Extension length
                top_k(config.top_k),            // Report the top_k alignments
                check_k(config.check_k),        // Check the scores of the check_k chains
                region_dist(config.region_dist),// Maximum distance for two regions to be called from the same region
                filter_dir(config.filter_dir),  // Use MEMs average length to filter the orientation of the reads
                dir_thr(config.dir_thr),        // Use MEMs average length distance to filter the orientation of the reads
                filter_seeds(config.filter_seeds),// Filter seed if occurs more than threshold
                n_seeds_thr(config.n_seeds_thr),// Filter seed if occurs more than threshold
                filter_freq(config.filter_freq),// Filter seed if it occurs with frequency greater than threshold
                freq_thr(config.freq_thr),      // Filter seed if it occurs with frequency greater than threshold
                max_iter(config.max_iter),      // Max number of iterations of the chaining algorithhm
                max_pred(config.max_pred),      // Max number of predecessor to be considered
                max_dist_x(config.max_dist_x),  // Max distance for two anchors to be chained
                max_dist_y(config.max_dist_y),  // Max distance for two anchors from the same read to be chained
                secondary_chains(config.secondary_chains), // Attempt to find secondary chains in paired-end setting
                smatch(config.smatch),          // Match score default
                smismatch(config.smismatch),    // Mismatch score default
                gapo(config.gapo),              // Gap open penalty
                gapo2(config.gapo2),            // Gap open penalty
                gape(config.gape),              // Gap extension penalty
                gape2(config.gape2),            // Gap extension penalty
                end_bonus(config.end_bonus),    // Bonus to add at the extension score to declare the alignment
                w(config.w),                    // Band width
                zdrop(config.zdrop),            // Zdrop enable
                max_penalty(std::max(smatch + smismatch, gapo + gape)), // Maximum penalty score
                forward_only(config.forward_only),
                report_mems(config.report_mems), // Report MEMs instead of read alignment
                csv(config.csv),                // Report MEM statistics in CSV file
                mem_finder(filename, config.min_len, config.filter_seeds, config.n_seeds_thr),
                ra(mem_finder.ra)
  {
    n = ra.getLen();




    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    std::string filename_idx = filename + idx.get_file_extension();
    verbose("Loading fasta index file: " + filename_idx);

    if (not file_exists(filename_idx))
      error("File not found: ", filename_idx);

    ifstream fs_idx(filename_idx);
    idx.load(fs_idx);
    fs_idx.close();

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Fasta index loading complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    verbose("Initialize the local aligner");
    t_insert_start = std::chrono::high_resolution_clock::now();

    ksw_gen_simple_mat(m,mat,smatch,-smismatch);

    chain_config.max_dist_x = config.max_dist_x;
    chain_config.max_dist_y = config.max_dist_y; 
    chain_config.max_iter = config.max_iter;
    chain_config.max_pred = config.max_pred;
    chain_config.min_chain_score = config.min_chain_score;
    chain_config.min_chain_length = config.min_chain_length;

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Local aligner initialization complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    verbose("Minimum MEM length: ", min_len);

    // memset(&ez_lc, 0, sizeof(ksw_extz_t));
    // memset(&ez_rc, 0, sizeof(ksw_extz_t));
    // memset(&ez, 0, sizeof(ksw_extz_t));
  }

  // Destructor
  ~aligner() 
  {
      // NtD
  }


  bool align(kseq_t *read, samFile* out, const sam_hdr_t *hdr)
  {
    // std::vector<mem_t> mems;
    alignment_t alignment(read);
    if(not align(alignment, hdr))
      alignment.set_sam_not_aligned();
    if (!report_mems)
      alignment.write(out, hdr);
    return alignment.aligned;
  }

  bool align(kseq_t *read, FILE* out, FILE* csv_out)
  {
    // std::vector<mem_t> mems;
    alignment_t alignment(read);
    if(not align(alignment, out))
      alignment.set_sam_not_aligned();
    if (!report_mems)
      alignment.write(out);
    if (csv)
      alignment.record_csv(csv_out);
    return alignment.aligned;
  }

  // Aligning unpaired sequences
  bool align(alignment_t &al, FILE* out = nullptr)
  {    
    MTIME_INIT(10);
    MTIME_START(0); // Timing helper
    MTIME_START(8);
    mem_finder.find_mems(al.read,al.mems, 0, MATE_1 | MATE_F);
    mem_finder.find_mems(&al.read_rev,al.mems, 0, MATE_1 | MATE_RC);
    MTIME_END(8);
    MTIME_START(9);
    mem_finder.populate_seeds(al.mems, report_mems);
    MTIME_END(9);

    if (csv)
      calculate_MEM_stats(al.mems, al.csv);
    if (filter_freq)
      seed_freq_filter(al.mems, freq_thr, al.csv);
    if (filter_seeds)
      seed_occ_filter(al.mems, n_seeds_thr, al.csv);

    // If reporting just the MEMs, at this point can directly write to the SAM file and skip the rest.
    if (report_mems)
    {
      MTIME_END(0); // Timing helper
      kseq_t read;
      nullptr_kseq_t(&read); // For FASTA reads to ensure nullptr for qual
      for (int i = 0; i < al.mems.size(); ++i){
        if (al.mems[i].mate & MATE_RC)
          copy_partial_kseq_t(&read, &al.read_rev, al.mems[i].idx, al.mems[i].len);
        else
          copy_partial_kseq_t(&read, al.read, al.mems[i].idx, al.mems[i].len);
        for (int j = 0; j < al.mems[i].occs.size(); ++j) {
            sam_t read_sam = sam_t(&read);
            read_sam.cigar = std::to_string(al.mems[i].len) + "M";
            const auto ref = idx.index(al.mems[i].occs[j]);
            read_sam.pos = ref.second + 1;
            read_sam.rname = ref.first;
            if (al.mems[i].mate & MATE_RC)
              read_sam.flag = SAM_SECONDARY_ALIGNMENT | SAM_REVERSED; 
            else
              read_sam.flag = SAM_SECONDARY_ALIGNMENT;
            write_sam(out, read_sam);
        }
        free_kseq_t(&read);
      }
      MTIME_TSAFE_MERGE;
      al.aligned = true;
      return al.aligned;
    }

    // Compute fraction of repetitive seeds
    al.frac_rep = compute_frac_rep(al.mems, al.read->seq.l, MATE_1);

    MTIME_END(0); //Timing helper
    MTIME_START(1); //Timing helper

    // Chain MEMs
    al.chained = find_chains(al.mems, al.anchors, al.chains, chain_config);

    MTIME_END(1);   //Timing helper
    MTIME_START(2); //Timing helper

    if (not al.chained)
    {
      MTIME_END(2); //Timing helper
      MTIME_TSAFE_MERGE;
      return false;
    }

    int32_t min_score = 20 + 8 * log(al.read->seq.l);

    // Compute the second best score
    std::vector<std::tuple<int32_t, size_t, size_t>> best_scores; // Score, position, chain index
    // get the occurrences of the top 4 best scores
    std::set<size_t> different_scores;
    size_t i = 0;

    std::vector<std::pair<size_t, size_t >> left_mem_vec;
    while (i < al.chains.size() and different_scores.size() < check_k)
    {
        different_scores.insert(al.chains[i].score);

        // Do pre-emptive check to see whether lifted left MEM pos of chain i match a previous chain, if so skip chain
        if (check_left_MEM(left_mem_vec, al, i)){
          ++i;
          al.csv.num_chains_skipped++;
          continue;
        }

        if (different_scores.size() < check_k)
        {
          // Align the chain
          auto chain = al.chains[i];
          // Reverse the chain order
          std::reverse(chain.anchors.begin(), chain.anchors.end());
          // Extract the anchors
          std::pair<ll, std::vector<size_t>> anchors_chain; 
          for(size_t i = 0; i < chain.anchors.size(); ++i)
            anchors_chain.second.push_back(chain.anchors[i]);
          // Compute the score of a chain.
          score_t score;
          if( (chain.mate & MATE_RC) )
            score = chain_score(anchors_chain, al.anchors, al.mems, min_score, &al.read_rev);
          else
            score = chain_score(anchors_chain, al.anchors, al.mems, min_score, al.read);
          
          score.lft = idx.lift(score.pos);
          // Check if there is another read scored in the same region
          bool replaced = false;
          for(size_t j = 0; j < best_scores.size(); ++j)
            if( (dist(std::get<1>(best_scores[j]),score.lft) < region_dist))
              if (score.score > std::get<0>(best_scores[j]) )
                if ( replaced )
                  best_scores[j] = std::make_tuple(0,0,i-1);
                else
                  best_scores[j] = std::make_tuple(score.score, score.lft, i++), replaced = true;
              else if (score.score <= std::get<0>(best_scores[j]) )
                j = best_scores.size(), replaced = true, i++;

          if( not replaced )
            best_scores.push_back(std::make_tuple(score.score, score.lft, i++));
        }
    }
    
    al.sub_n = best_scores.size() -1; 

    while (best_scores.size() < 2)
      best_scores.push_back(std::make_tuple(0, 0, al.chains.size()));

    std::sort(best_scores.begin(), best_scores.end(), std::greater<std::tuple<int32_t,size_t, size_t>>());

    assert(best_scores.size() > 1);

    if (std::get<0>(best_scores[0]) < min_score)
    {
      MTIME_END(2); //Timing helper
      MTIME_TSAFE_MERGE;
      return false;
    }

    al.best_score = true;

    al.score2 = std::get<0>(best_scores[1]);


    { 
      i = std::get<2>(best_scores[0]);
      // Align the chain
      auto chain = al.chains[i];
      // Reverse the chain order
      std::reverse(chain.anchors.begin(), chain.anchors.end());
      // Extract the anchors
      std::pair<ll, std::vector<size_t>> anchors_chain; 
      for(size_t i = 0; i < chain.anchors.size(); ++i)
        anchors_chain.second.push_back(chain.anchors[i]);
      // Compute the score of a chain.
      if( (chain.mate & MATE_RC) )
      {
        al.score = chain_score(anchors_chain, al.anchors, al.mems, min_score, &al.read_rev, false, al.score2, 1, nullptr, &al.sam);
        al.sam.read = &al.read_rev;
        al.sam.flag |= SAM_REVERSED;
      }
      else
        al.score = chain_score(anchors_chain, al.anchors, al.mems, min_score, al.read, false, al.score2, 0, nullptr, &al.sam);
    }

    al.aligned = (al.score.score >= min_score);

    MTIME_END(2); //Timing helper
    MTIME_TSAFE_MERGE;
    if (not al.aligned)
      return false;
    

    return al.aligned;
  
  }

  // Check whether the left most MEM lifted over REF position is within certain distance of previous left most MEM lifted over REF positions.
  // Pre-emptive check to see whether the chain we are going to fully align will liftover to the position of a previous chain. 
  // If this is the case, then do not spend time to do the full alignment if its likely going to be discarded anyways.
  bool check_left_MEM(
    std::vector<std::pair<size_t, size_t >>& left_mem_vec,
    alignment_t& al,
    size_t i)
  {
    // Get the chain
    auto& chain = al.chains[i];
    // Reverse the chain order [Have to undo this afterwards]
    chain.reverse();

    // Extract the ref position of the leftmost anchor of the current chain
    size_t left_mem_pos, left_mem_ref_pos;
    
    // Extract the leftmost MEM of the chain
    for(size_t j = 0; j < chain.anchors.size(); ++j)
    {
      size_t anchor_id = chain.anchors[j];
      left_mem_pos = al.mems[al.anchors[anchor_id].first].occs[al.anchors[anchor_id].second];
      const auto lift = idx.lift(left_mem_pos);
      const auto lft_ref = idx.index(lift);
      left_mem_ref_pos = lft_ref.second + 1;
      break;
    }

    // Check to see if the leftmost MEM of the chain has lifted over position corresponding to previously seen chain
    bool discovered = false;
    for (size_t j = 0; j < left_mem_vec.size(); ++j){
      if (dist(left_mem_vec[j].first, left_mem_ref_pos) < region_dist){
        if (left_mem_vec[j].second == al.chains[i].score){
          discovered = true;
        }
      }
    }

    // If discovered then likely current chain is the same as another previously seen chain
    if (discovered){
      chain.reset(); // Get original chain order
      return true;
    }
    else{
      chain.reset(); // Get the original chain order
      left_mem_vec.push_back(std::make_pair(left_mem_ref_pos, al.chains[i].score));
      return false;
    }
  }

  typedef struct orphan_paired_score_t{
    int32_t tot = 0;
    int64_t dist = 0;
    score_t m1;
    score_t m2; 
    size_t chain_i = 0; 
    std::pair<size_t,size_t> pos = std::make_pair(0,0); // Position of the mate

    bool operator<=(const orphan_paired_score_t &other)
    {
      return (this.tot < other.tot) or
             (this.tot == other.tot and this.m1.lft < this.m1.lft) or
             (this.tot == other.tot and this.m1.lft == this.m1.lft and this.m2.lft <= this.m2.lft);
    }

    friend bool operator>(const orphan_paired_score_t &lhs, const orphan_paired_score_t &rhs)
    {
      return (lhs.tot > rhs.tot) or
             (lhs.tot == rhs.tot and lhs.m1.lft > rhs.m1.lft) or
             (lhs.tot == rhs.tot and lhs.m1.lft == rhs.m1.lft and lhs.m2.lft > rhs.m2.lft);
    }

  } orphan_paired_score_t;

  typedef struct paired_score_t{
    int32_t tot = 0;    // Total score
    int64_t dist = 0;   // Pair distance
    score_t m1;         // Score of the first read in the pair
    score_t m2;         // Score of the second read in the pair
    size_t chain_i = 0; // Chain index
    bool paired = false;


    paired_score_t& operator=(orphan_paired_score_t const& copy)
    {
      tot = copy.tot;
      dist = copy.dist;
      m1 = copy.m1;
      m2 = copy.m2;
      chain_i = copy.chain_i;
      return *this;
    }

    bool operator<=(const paired_score_t &other)
    {
      return (this.tot < other.tot) or
             (this.tot == other.tot and this.m1.lft < this.m1.lft) or
             (this.tot == other.tot and this.m1.lft == this.m1.lft and this.m2.lft <= this.m2.lft);
    }

    friend bool operator>(const paired_score_t& lhs, const paired_score_t& rhs)
    {
      return (lhs.tot > rhs.tot) or
             (lhs.tot == rhs.tot and lhs.m1.lft > rhs.m1.lft) or
             (lhs.tot == rhs.tot and lhs.m1.lft == rhs.m1.lft and lhs.m2.lft > rhs.m2.lft);
    }
  } paired_score_t;

  typedef struct paired_alignment_t{
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

    csv_t csv_m1;
    csv_t csv_m2;

    // TODO: precompute the nt4 version of the mates

    int32_t min_score_m1 = 0;
    int32_t min_score_m2 = 0;
    int32_t min_score = 0;
    paired_score_t score;
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
    std::vector<paired_score_t> best_scores;
    
    paired_alignment_t(kseq_t *mate1_, kseq_t *mate2_, double mean_ = 0.0, double std_dev_ = 0.0) : 
      mate1(mate1_),
      mate2(mate2_),
      sam_m1(mate1),
      sam_m2(mate2),
      csv_m1(mate1),
      csv_m2(mate2),
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

    // MEMs within mate pair are stored together, therefore do not need csv_m2
    void record_csv(FILE* out)
    {
      write_csv(out, csv_m1);
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

  // Aligning pair-ended batched sequences
  // Return true if the fragment model has been learned.
  // Assumes alignment to be allocated with the exact number of elements
  bool learn_fragment_model(kpbseq_t *batch)
  {
    size_t n_aligned = 0;

    // Computing Mean and Variance using Welford's algorithm
    size_t count = 0;   // Number of samples
    double mean = 0.0;  // Accumulates the mean
    double m2 = 0.0;    // Accumulates the squared distance from the mean

    int l = batch->mate1->l;
    for (size_t i = 0; i < l; ++i)
    {
      // paired_alignment_t alignment(&batch->mate1->buf[i], &batch->mate2->buf[i]);
      paired_alignment_t alignment;
      alignment.init(&batch->mate1->buf[i], &batch->mate2->buf[i]);
      if(align(alignment, false) and ((not alignment.second_best_score) or ((alignment.best_scores[0].tot - alignment.best_scores[1].tot) > ins_learning_score_gap_threshold)))
      {
        // Get stats
        // mate_abs_distance.push_back((double)(alignment.sam_m1.tlen >= 0?alignment.sam_m1.tlen :-alignment.sam_m1.tlen));
        // double value = (double)(alignment.sam_m1.tlen >= 0?alignment.sam_m1.tlen :-alignment.sam_m1.tlen);
        double value = (double)(alignment.best_scores[0].dist);
        double delta = value - mean;
        mean += delta / (++count);
        m2 += delta * (value - mean);
      }
    }

    // Computes stats
    double variance = m2/count;
    double sampleVariance = m2/(count-1);
    double std_dev = sqrt(variance);

    __ins_mtx.lock();
    if (not ins_learning_complete)
    {
      if (ins_count > 0)
      {
        size_t t_count = ins_count + count;
        double delta = ins_mean - mean;
        ins_m2 += m2 + (delta * delta * ins_count * count) / t_count;
        ins_mean = (ins_count * ins_mean + count * mean)/ t_count;
        ins_count = t_count;
      }
      else
      {
        ins_mean = mean;
        ins_std_dev = std_dev;
        ins_m2 = m2;
        ins_count = count;
      }
      verbose("Number of high quality samples processed so far: ", ins_count);
      ins_learning_complete = ins_learning_complete or (ins_count >= ins_learning_n);
      if (ins_learning_complete)
      {
        ins_variance = ins_m2 / ins_count;
        ins_sample_variance = ins_m2 / (ins_count - 1);
        ins_std_dev = sqrt(ins_variance);

        verbose("Insertion size estimation complete!");
        verbose("Number of high quality samples: ", ins_count);
        verbose("                          Mean: ", ins_mean);
        verbose("                      Variance: ", ins_variance);
        verbose("               Sample Variance: ", ins_sample_variance);
        verbose("            Standard Deviation: ", ins_std_dev);
      }
    }
    __ins_mtx.unlock();

    return ins_learning_complete;
  }

  // Aligning pair-ended batched sequences
  statistics_t align(kpbseq_t *batch, FILE *out, FILE *csv_out)
  {
    statistics_t stats;

    int l = batch->mate1->l;
    for (size_t i = 0; i < l; ++i)
    {
      ++stats.processed_reads;
      // paired_alignment_t alignment(&batch->mate1->buf[i], &batch->mate2->buf[i]);
      paired_alignment_t alignment(&batch->mate1->buf[i], &batch->mate2->buf[i]);
      alignment.mean = ins_mean;
      alignment.std_dev = ins_std_dev;
      if(not align(alignment, true, out) and alignment.chained)
      {
        ++ stats.orphan_reads;
        orphan_recovery(alignment, ins_mean, ins_std_dev);
        if (alignment.aligned) ++stats.orphan_recovered_reads;
      }
      // Write alignment to file
      if (!report_mems){
        alignment.write(out);
      }
      // Write MEM stats to CSV file
      if (csv){
        alignment.record_csv(csv_out);
      }

      if(alignment.aligned) ++stats.aligned_reads;
    }
    
    return stats;
  }

  // // Aligning pair-ended batched sequences
  // size_t align(kpbseq_t *batch, FILE *out)
  // {
  //   size_t n_aligned = 0;

  //   std::vector<paired_alignment_t> alignments(batch->mate1->l);
  //   std::vector<double> mate_abs_distance;

  //   // Computing Mean and Variance using Welford's algorithm
  //   size_t count = 0;   // Number of samples
  //   double mean = 0.0;  // Accumulates the mean
  //   double m2 = 0.0;    // Accumulates the squared distance from the mean


  //   int l = batch->mate1->l;
  //   for (size_t i = 0; i < l; ++i)
  //   {
  //     // paired_alignment_t alignment(&batch->mate1->buf[i], &batch->mate2->buf[i]);
  //     paired_alignment_t& alignment = alignments[i];
  //     alignment.init(&batch->mate1->buf[i], &batch->mate2->buf[i]);
  //     if(align(alignment))
  //     {
  //       // Get stats
  //       // mate_abs_distance.push_back((double)(alignment.sam_m1.tlen >= 0?alignment.sam_m1.tlen :-alignment.sam_m1.tlen));
  //       double value = (double)(alignment.sam_m1.tlen >= 0?alignment.sam_m1.tlen :-alignment.sam_m1.tlen);
  //       double delta = value - mean;
  //       mean += delta / (++count);
  //       m2 += delta * (value - mean);
  //     }
  //     // Store the alignment informations
  //     // alignments.push_back(alignment);
  //   }

  //   // Computes stats
  //   double variance = m2/count;
  //   double sampleVariance = m2/(count-1);
  //   // double variance = 0.0;
  //   // double mean = mate_abs_distance[0];
  //   // for(size_t i = 0 ; i < l; ++i )
  //   // {
  //   //   mean += mate_abs_distance[i];
  //   //   double diff = ((i+1) * mate_abs_distance[i] ) - mean;
  //   //   variance += (diff * diff) / ((i + 1.0 ) * i);
  //   // }

  //   // variance /= (double)(l-1);
  //   // mean /= (double)l;
  //   double std_dev = sqrt(variance);

  //   // Perform local serach for unaligned mates.
  //   for(auto& alignment : alignments)
  //   {
  //     if(not alignment.aligned and alignment.chained)
  //       orphan_recovery(alignment, mean, std_dev);
  //     // Write alignment to file
  //     alignment.write(out);
  //     if(alignment.aligned) ++n_aligned;
  //   }
    
  //   return n_aligned;
  // }

  // Aligning pair-ended sequences
  bool align(kseq_t *mate1, kseq_t *mate2, FILE *out, FILE *csv_out)
  {
    paired_alignment_t alignment(mate1, mate2);
    if(not align(alignment, true, out))
      alignment.set_sam_not_aligned();
    if (!report_mems){
      alignment.write(out);
    }
    if (csv){
      alignment.record_csv(csv_out);
    }
    return alignment.aligned;
  }

  // Aligning pair-ended sequences
  bool align(paired_alignment_t &al, bool finalize = true, FILE* out = nullptr)
  { 
    MTIME_INIT(10);   
    MTIME_START(0); // Timing helper

    // Find MEMs
    if ( filter_dir )
    {
      MTIME_START(8); // Timing helper
      // Direction 1
      // find_seeds(al.mate1, al.mems, 0, MATE_1 | MATE_F);
      // find_seeds(&al.mate2_rev, al.mems, al.mate1->seq.l, MATE_2 | MATE_RC);
      mem_finder.find_mems(al.mate1, al.mems, 0, MATE_1 | MATE_F);
      mem_finder.find_mems(&al.mate2_rev, al.mems, al.mate1->seq.l, MATE_2 | MATE_RC);
      
      al.n_mems_dir1 = al.mems.size();
      al.n_seeds_dir1 = 0;
      for ( size_t i = 0; i < al.mems.size(); ++i )
      {
        al.n_seeds_dir1 += al.mems[i].occs.size();
        al.avg_seed_length_dir1 += al.mems[i].len;
        al.avg_w_seed_length_dir1 += al.mems[i].len * al.mems[i].occs.size();
        al.armonic_avg_seed_length_dir1 += (double)al.mate1->seq.l/(double)al.mems[i].len;
      }
      if (al.n_mems_dir1 > 0) 
      {
        al.avg_seed_length_dir1 /= al.n_mems_dir1;
        al.avg_w_seed_length_dir1 /= al.n_seeds_dir1;
        al.armonic_avg_seed_length_dir1 = (double)al.n_mems_dir1 / al.armonic_avg_seed_length_dir1;
      }
      // Direction 2
      // find_seeds(al.mate2, al.mems, 0, MATE_2 | MATE_F);
      // find_seeds(&al.mate1_rev, al.mems, al.mate2->seq.l, MATE_1 | MATE_RC );
      mem_finder.find_mems(al.mate2, al.mems, 0, MATE_2 | MATE_F);
      mem_finder.find_mems(&al.mate1_rev, al.mems, al.mate2->seq.l, MATE_1 | MATE_RC );

      al.n_mems_dir2 = al.mems.size() - al.n_mems_dir1;
      al.n_seeds_dir2 = 0;
      for ( size_t i = al.n_mems_dir1; i < al.mems.size(); ++i )
      {
        al.n_seeds_dir2 += al.mems[i].occs.size();
        al.avg_seed_length_dir2 += al.mems[i].len;
        al.avg_w_seed_length_dir2 += al.mems[i].len * al.mems[i].occs.size();
        al.armonic_avg_seed_length_dir2 += (double)al.mate1->seq.l / (double)al.mems[i].len;
      }
      if (al.n_mems_dir2 > 0) 
      {
        al.avg_seed_length_dir2 /= al.n_mems_dir2;
        al.avg_w_seed_length_dir2 /= al.n_seeds_dir2;
        al.armonic_avg_seed_length_dir2 = (double)al.n_mems_dir2 / al.armonic_avg_seed_length_dir2;
      }

      // Make a decision based on the MEMs
      if ((al.avg_seed_length_dir1 > al.avg_seed_length_dir2) and ((al.avg_seed_length_dir1 - al.avg_seed_length_dir2) > dir_thr))
        al.mems.erase(al.mems.begin() + al.n_mems_dir1, al.mems.end());
      if ((al.avg_seed_length_dir2 > al.avg_seed_length_dir1) and ((al.avg_seed_length_dir2 - al.avg_seed_length_dir1) > dir_thr))
        al.mems.erase(al.mems.begin(), al.mems.begin()  + al.n_mems_dir1);

      MTIME_END(8); //Timing helper
      MTIME_START(9); //Timing helper
      
      mem_finder.populate_seeds(al.mems, report_mems);
      if (csv)
        calculate_MEM_stats(al.mems, al.csv_m1);

      al.n_seeds_dir1 = 0;
      al.n_seeds_dir2 = 0;
      for (size_t i = 0; i < al.mems.size(); ++i)
      {
        al.n_seeds_dir1 += al.mems[i].occs.size();
        al.avg_w_seed_length_dir1 += al.mems[i].len * al.mems[i].occs.size();
      }
      for (size_t i = al.n_mems_dir1; i < al.n_mems_dir1 + al.n_mems_dir2; ++i)
      {
        al.n_seeds_dir2 += al.mems[i].occs.size();
        al.avg_w_seed_length_dir2 += al.mems[i].len * al.mems[i].occs.size();
      }
      if (al.n_mems_dir1 > 0) 
      {
        al.avg_w_seed_length_dir1 /= al.n_seeds_dir1;
      }
      if (al.n_mems_dir2 > 0) 
      {
        al.avg_w_seed_length_dir2 /= al.n_seeds_dir2;
      }
      MTIME_END(9); //Timing helper

      if (filter_freq)
        seed_freq_filter(al.mems, freq_thr, al.csv_m1);
      if (filter_seeds)
        seed_occ_filter(al.mems, n_seeds_thr, al.csv_m1);
    }
    else
    {
      MTIME_START(8);
      mem_finder.find_mems(al.mate1, al.mems, 0, MATE_1 | MATE_F);
      mem_finder.find_mems(&al.mate1_rev, al.mems, al.mate2->seq.l, MATE_1 | MATE_RC);
      mem_finder.find_mems(al.mate2, al.mems, 0, MATE_2 | MATE_F);
      mem_finder.find_mems(&al.mate2_rev, al.mems, al.mate1->seq.l, MATE_2 | MATE_RC);
      MTIME_END(8);
      MTIME_START(9);
      mem_finder.populate_seeds(al.mems, report_mems);
      MTIME_END(9);

      if (csv)
        calculate_MEM_stats(al.mems, al.csv_m1);
      if (filter_freq)
        seed_freq_filter(al.mems, freq_thr, al.csv_m1);
      if (filter_seeds)
        seed_occ_filter(al.mems, n_seeds_thr, al.csv_m1);
    }
    // find_mems(al.mate1, al.mems, 0, MATE_1 | MATE_F);
    // find_mems(&al.mate1_rev, al.mems, al.mate2->seq.l, MATE_1 | MATE_RC );
    // find_mems(al.mate2, al.mems, 0, MATE_2 | MATE_F);
    // find_mems(&al.mate2_rev, al.mems, al.mate1->seq.l, MATE_2 | MATE_RC);

    // If reporting just the MEMs, at this point can directly write to the SAM file and skip the rest.
    if (report_mems && out != nullptr)
    {
      MTIME_END(0); //Timing helper
      // Cannot modify mate 1 and mate 2 directly. Instead have to create new kseq_t variables to do this. 
      kseq_t mem_m1_read;
      kseq_t mem_m2_read;
      nullptr_kseq_t(&mem_m1_read); // For FASTA reads to ensure nullptr for qual
      nullptr_kseq_t(&mem_m2_read); // For FASTA reads to ensure nullptr for qual

      // Loop throught the MEMs of the paired-end read and write to the SAM file
      for (int i = 0; i < al.mems.size(); ++i){
        if (al.mems[i].mate == (MATE_1 | MATE_F) || al.mems[i].mate == (MATE_1 | MATE_RC)){
          if (al.mems[i].mate & MATE_RC)
            copy_partial_kseq_t(&mem_m1_read, &al.mate1_rev, al.mems[i].idx, al.mems[i].len);
          else
            copy_partial_kseq_t(&mem_m1_read, al.mate1, al.mems[i].idx, al.mems[i].len);
          for (int j = 0; j < al.mems[i].occs.size(); ++j){
            sam_t m1_sam = sam_t(&mem_m1_read);
            m1_sam.cigar = std::to_string(al.mems[i].len) + "M";
            const auto ref = idx.index(al.mems[i].occs[j]);
            m1_sam.pos = ref.second + 1;
            m1_sam.rname = ref.first;
            if (al.mems[i].mate & MATE_RC)
              m1_sam.flag = SAM_SECONDARY_ALIGNMENT | SAM_REVERSED;
            else
              m1_sam.flag = SAM_SECONDARY_ALIGNMENT;
            write_sam(out, m1_sam);
          }
          free_kseq_t(&mem_m1_read);
        }
        else{
          if (al.mems[i].mate & MATE_RC)
            copy_partial_kseq_t(&mem_m2_read, &al.mate2_rev, al.mems[i].idx, al.mems[i].len);
          else
            copy_partial_kseq_t(&mem_m2_read, al.mate2, al.mems[i].idx, al.mems[i].len);
          for (int j = 0; j < al.mems[i].occs.size(); ++j) {
            sam_t m2_sam = sam_t(&mem_m2_read);
            m2_sam.cigar = std::to_string(al.mems[i].len) + "M";
            const auto ref = idx.index(al.mems[i].occs[j]);
            m2_sam.pos = ref.second + 1;
            m2_sam.rname = ref.first;
            if (al.mems[i].mate & MATE_RC)
              m2_sam.flag = SAM_SECONDARY_ALIGNMENT | SAM_REVERSED;
            else
              m2_sam.flag = SAM_SECONDARY_ALIGNMENT;
            write_sam(out, m2_sam);
          }
          free_kseq_t(&mem_m2_read);
        }
      }

      MTIME_TSAFE_MERGE;
      al.aligned = true;
      return al.aligned;
    }

    // Compute fraction of repetitive seeds
    al.frac_rep_m1 = compute_frac_rep(al.mems, al.mate1->seq.l, MATE_1);
    al.frac_rep_m2 = compute_frac_rep(al.mems, al.mate2->seq.l, MATE_2);

    MTIME_END(0); //Timing helper
    MTIME_START(1); //Timing helper

    // Chain MEMs
    if (secondary_chains){
      al.chained = find_chains_secondary(al.mems, al.anchors, al.chains, chain_config);
    }
    else{
      al.chained = find_chains(al.mems, al.anchors, al.chains, chain_config);
    }

    MTIME_END(1);   //Timing helper
    MTIME_START(2); //Timing helper

    if (not al.chained)
    {
      MTIME_END(2); //Timing helper
      MTIME_TSAFE_MERGE;
      return false;
    }


    // int32_t min_score = al.min_score;

    // // // Compute the second best score
    // // (0) score (1) lift m1 (2) lift m2 (3) index in chains list (4) score m1 (5) score m2
    // std::vector<std::tuple<int32_t, size_t, size_t, size_t, size_t, size_t>> best_scores;
    // // get the occurrences of the top 4 best scores
    // std::set<size_t> different_scores;
    // size_t i = 0;
    // while (i < al.chains.size() and different_scores.size() < check_k)
    // {
    //     different_scores.insert(al.chains[i].score);
    //     if (different_scores.size() < check_k)
    //     {
    //       // Align the chain
    //       auto chain = al.chains[i];
    //       // Reverse the chain order
    //       std::reverse(chain.anchors.begin(), chain.anchors.end());
    //       // Compute the score of a chain.
    //       paired_score_t score;
    //       // if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
    //       //   score = paired_chain_score(chain, al.anchors, al.mems, min_score, al.min_score_m1, al.min_score_m2, al.mate1, &al.mate2_rev);
    //       // else
    //       //   score = paired_chain_score(chain, al.anchors, al.mems, min_score, al.min_score_m1, al.min_score_m2, &al.mate1_rev, al.mate2);
    //       score = paired_chain_score(chain, al);

    //       score.m1.lft = idx.lift(score.m1.pos);
    //       score.m2.lft = idx.lift(score.m2.pos);

    //       // Check if there is another read scored in the same region
    //       bool replaced = false;
    //       for(size_t j = 0; j < best_scores.size(); ++j)
    //       {
    //         auto d1 = dist(std::get<1>(best_scores[j]),score.m1.lft);
    //         auto d2 = dist(std::get<2>(best_scores[j]),score.m2.lft);
    //         if( (dist(std::get<1>(best_scores[j]),score.m1.lft) < region_dist) and 
    //             (dist(std::get<2>(best_scores[j]),score.m2.lft) < region_dist) )
    //           if ( score.tot > std::get<0>(best_scores[j]) )
    //             if ( replaced )
    //               best_scores[j] = std::make_tuple(0,0,0,i-1,0,0);
    //             else
    //               best_scores[j] = std::make_tuple(score.tot, score.m1.lft, score.m2.lft, i++, score.m1.score, score.m2.score), replaced = true;
    //           else if ( score.tot == std::get<0>(best_scores[j]) )
    //             j = best_scores.size(), replaced = true, i++;
    //       }
    //       if( not replaced )
    //         best_scores.push_back(std::make_tuple(score.tot, score.m1.lft, score.m2.lft, i++, score.m1.score, score.m2.score));
    //     }
    // }
    // if (best_scores.size() < 2)
    //   best_scores.push_back(std::make_tuple(0, 0, 0, al.chains.size(),0,0));



    // std::sort(best_scores.begin(), best_scores.end(), std::greater<std::tuple<int32_t, size_t, size_t, size_t, size_t, size_t>>());

    get_best_scores(al, check_k);

    auto& best_scores = al.best_scores;

    assert(best_scores.size() >= 1);

    if (best_scores[0].tot < al.min_score)
    {
      MTIME_END(2); //Timing helper
      MTIME_TSAFE_MERGE;
      return false;
    }

    uint8_t direction = 0;
    const size_t mate = al.chains[best_scores[0].chain_i].mate;
    if ( ((mate & MATE_1) and (mate & MATE_F)) or ((mate & MATE_2) and (mate & MATE_RC)))
      direction = 1;
    else
      direction = 2;

    // // Compute sub-optimal hits (https://github.com/bwa-mem2/bwa-mem2/blob/edc703f883e8aaed83067100d8e54e0e9e810ef5/src/bwamem.cpp#L1312)
    // // From BWA's manual: BWA will not search for suboptimal hits with a score lower than (bestScore-misMsc).
    // size_t k = 1;
    // al.sub_n = 0;
    // while( k < best_scores.size() and best_scores[k++].tot >= (best_scores[0].tot - max_penalty))
    //   ++al.sub_n;

    // al.best_score = true;

    // al.score2 = best_scores[1].tot;
    // al.score2_m1 = best_scores[1].m1.score;
    // al.score2_m2 = best_scores[1].m2.score;

    // al.second_best_score = (al.score2 >= al.min_score);

    if (finalize) { // Forward case
      size_t i = best_scores[0].chain_i;
      // // Align the chain
      // auto chain = al.chains[i];
      // // Reverse the chain order
      // std::reverse(chain.anchors.begin(), chain.anchors.end());
      // // Compute the score of a chain.
      // // if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
      // //   al.score = paired_chain_score(chain, al.anchors, al.mems, min_score, al.min_score_m1, al.min_score_m2, al.mate1, &al.mate2_rev, false, al.score2, 0, &al.sam_m1, &al.sam_m2);
      // // else
      // //   al.score = paired_chain_score(chain, al.anchors, al.mems, min_score, al.min_score_m1, al.min_score_m2, &al.mate1_rev, al.mate2, false, al.score2, 1, &al.sam_m1, &al.sam_m2);
      // al.score = paired_chain_score(chain, al, false);
      al.score = paired_chain_score(al, i, false);
      al.aligned = (al.score.tot >= al.min_score);
    } else
      al.aligned = (best_scores[0].tot >= al.min_score);

    MTIME_END(2); //Timing helper
    MTIME_TSAFE_MERGE;

    return al.aligned;
  }


  void get_best_scores(paired_alignment_t& al, size_t k)
  {
    // // Compute the second best score
    // get the occurrences of the top 4 best scores
    std::set<size_t> different_scores;
    size_t i = 0;
    std::vector<std::pair<size_t, size_t >> mate1_left_mem_vec;
    std::vector<std::pair<size_t, size_t >> mate2_left_mem_vec;

    while (i < al.chains.size() and different_scores.size() < k)
    {
        different_scores.insert(al.chains[i].score);

        // Do pre-emptive check to see whether lifted left MEM pos of mate1 and mate2 chain i match a previous chain, if so skip chain
        if (check_paired_left_MEM(mate1_left_mem_vec, mate2_left_mem_vec, al, i)){
          ++i;
          al.csv_m1.num_chains_skipped++;
          continue;
        }

        if (different_scores.size() < k)
        {
          // Align the chain
          paired_score_t score = paired_chain_score(al, i);

          // Check if there is another read scored in the same region
          if (score.tot >= al.min_score)
          {
            bool replaced = false;
            for(size_t j = 0; j < al.best_scores.size(); ++j)
            {
              paired_score_t zero;
              zero.chain_i = i;
              auto d1 = dist(al.best_scores[j].m1.lft,score.m1.lft);
              auto d2 = dist(al.best_scores[j].m2.lft,score.m2.lft);
              if ((dist(al.best_scores[j].m1.lft, score.m1.lft) < region_dist) and
                  (dist(al.best_scores[j].m2.lft, score.m2.lft) < region_dist))
                if ( score.tot > al.best_scores[j].tot )
                  if ( replaced )
                    al.best_scores[j] = zero;
                  else
                    al.best_scores[j] = score, replaced = true;
                else if ( score.tot <= al.best_scores[j].tot )
                  j = al.best_scores.size(), replaced = true;
            }
            if( not replaced )
              al.best_scores.push_back(score);
          }

          ++i;
        }
    }

    paired_score_t zero;
    zero.chain_i = al.chains.size();

    while (al.best_scores.size() < 2)
      al.best_scores.push_back(zero);

    // Sort the chains by score
    std::sort(al.best_scores.begin(), al.best_scores.end(), std::greater<paired_score_t>());

    // Compute sub-optimal hits (https://github.com/bwa-mem2/bwa-mem2/blob/edc703f883e8aaed83067100d8e54e0e9e810ef5/src/bwamem.cpp#L1312)
    // From BWA's manual: BWA will not search for suboptimal hits with a score lower than (bestScore-misMsc).
    size_t j = 1;
    al.sub_n = 0;
    while (j < al.best_scores.size() and al.best_scores[j++].tot >= (al.best_scores[0].tot - max_penalty))
      ++al.sub_n;

    al.best_score = true;

    al.score2 = al.best_scores[1].tot;
    al.score2_m1 = al.best_scores[1].m1.score;
    al.score2_m2 = al.best_scores[1].m2.score;

    al.second_best_score = (al.score2 >= al.min_score);
  }

  void update_best_scores(paired_alignment_t& al)
  {
    // // Compute the second best score
    for (size_t i = 0; i < al.best_scores.size(); ++i)
    {
      paired_score_t& score = al.best_scores[i];

      double ns = 0.0;
      if (al.std_dev > 0.0)
        ns = (score.dist - al.mean) / al.std_dev;
      score.tot = (int)(score.m1.score + score.m2.score + .721 * log(2. * erfc(fabs(ns) * M_SQRT1_2)) * smatch + .499);
      if (score.tot < 0)
        score.tot = 0;
    }

    // Sort the chains by score
    std::sort(al.best_scores.begin(), al.best_scores.end(), std::greater<paired_score_t>());

    // Compute sub-optimal hits (https://github.com/bwa-mem2/bwa-mem2/blob/edc703f883e8aaed83067100d8e54e0e9e810ef5/src/bwamem.cpp#L1312)
    // From BWA's manual: BWA will not search for suboptimal hits with a score lower than (bestScore-misMsc).
    size_t k = 1;
    al.sub_n = 0;
    while (k < al.best_scores.size() and al.best_scores[k++].tot >= (al.best_scores[0].tot - max_penalty))
      ++al.sub_n;

    al.best_score = true;

    al.score2 = al.best_scores[1].tot;
    al.score2_m1 = al.best_scores[1].m1.score;
    al.score2_m2 = al.best_scores[1].m2.score;

    al.second_best_score = (al.score2 >= al.min_score);
  }


  // Check whether the left most MEM lifted over REF position is within certain distance of previous left most MEM lifted over REF positions.
  // Pre-emptive check to see whether the chain we are going to fully align will liftover to the position of a previous chain. 
  // If this is the case, then do not spend time to do the full alignment if its likely going to be discarded anyways.
  bool check_paired_left_MEM(
    std::vector<std::pair<size_t, size_t >>& mate1_left_mem_vec,
    std::vector<std::pair<size_t, size_t >>& mate2_left_mem_vec,
    paired_alignment_t& al,
    size_t i)
  {
    // Get the chain
    auto& chain = al.chains[i];
    // Reverse the chain order [Have to undo this afterwards]
    chain.reverse();

    // Extract the ref position of the leftmost anchor of the current chain for each mate
    size_t mate1_left_mem_pos, mate1_left_mem_ref_pos, mate2_left_mem_pos, mate2_left_mem_ref_pos;
    
    // Extract the leftmost MEM of the chain from mate 1 
    for(size_t j = 0; j < chain.anchors.size(); ++j)
    {
      size_t anchor_id = chain.anchors[j];
      if ((al.mems[al.anchors[anchor_id].first].mate & MATE_2) == 0){
        mate1_left_mem_pos = al.mems[al.anchors[anchor_id].first].occs[al.anchors[anchor_id].second];
        const auto lift = idx.lift(mate1_left_mem_pos);
        const auto lft_ref = idx.index(lift);
        mate1_left_mem_ref_pos = lft_ref.second + 1;
        break;
      }
    }

    // Extract the leftmost MEM of the chain from mate 2
    for(size_t j = 0; j < chain.anchors.size(); ++j)
    {
      size_t anchor_id = chain.anchors[j];
      if ((al.mems[al.anchors[anchor_id].first].mate & MATE_2) != 0){
        mate2_left_mem_pos = al.mems[al.anchors[anchor_id].first].occs[al.anchors[anchor_id].second];
        const auto lift = idx.lift(mate2_left_mem_pos);
        const auto lft_ref = idx.index(lift);
        mate2_left_mem_ref_pos = lft_ref.second + 1;
        break;
      }
    }

    // Check to see if the leftmost MEM of each chain has lifted over position corresponding to previously seen chain
    bool discovered = false;
    for (size_t j = 0; j < mate1_left_mem_vec.size(); ++j){
      if ((dist(mate1_left_mem_vec[j].first, mate1_left_mem_ref_pos) < region_dist) 
      and (dist(mate2_left_mem_vec[j].first, mate2_left_mem_ref_pos) < region_dist)){
        if (mate1_left_mem_vec[j].second == al.chains[i].score){
          discovered = true;
        }
      }
    }

    // If discovered then likely current chain is the same as another previously seen chain
    if (discovered){
      chain.reset(); // Get original chain order
      return true;
    }
    else{
      chain.reset(); // Get the original chain order
      mate1_left_mem_vec.push_back(std::make_pair(mate1_left_mem_ref_pos, al.chains[i].score));
      mate2_left_mem_vec.push_back(std::make_pair(mate2_left_mem_ref_pos, al.chains[i].score));
      return false;
    }
  }


  bool orphan_recovery(
    paired_alignment_t& al, 
    const double mean, 
    const double std_dev)
  {
    MTIME_INIT(8);   
    MTIME_START(7); //Timing helper
    // We need to use the local alignment of ksw to perform local search
    
    // For all good chaining of both mates of length at least min_length, find the possible distance of the mate
    // and 

    // Compute the second best score
    // (0) score (1) pos (2) lift m1 (3) lift m2 (4) index in chains list (5) score m1 (6) score m2
    // std::vector<std::tuple<int32_t, std::pair<size_t,size_t>, size_t, size_t, size_t, size_t, size_t>> best_scores;
    std::vector<orphan_paired_score_t> best_scores;
    std::pair<size_t,size_t> pair_zero = {0,0};
    // std::vector<std::pair<std::pair<int32_t, std::pair<size_t,size_t> >, size_t>> best_scores;
    size_t i = 0;
    while (i < al.chains.size())
    {
      // // Align the chain
      // auto chain = al.chains[i];
      // // Reverse the chain order
      // std::reverse(chain.anchors.begin(), chain.anchors.end());
      // // Compute the score of a chain.
      // orphan_paired_score_t score;
      // // if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
      // //   score = paired_chain_orphan_score(chain, al.anchors, al.mems, al.min_score, al.min_score_m1, al.min_score_m2, al.mate1, &al.mate2_rev, mean, std_dev);
      // // else
      // //   score = paired_chain_orphan_score(chain, al.anchors, al.mems, al.min_score, al.min_score_m1, al.min_score_m2, &al.mate1_rev, al.mate2, mean, std_dev);
      // score = paired_chain_orphan_score(chain, al, mean, std_dev);
      
      // score.m1.lft = idx.lift(score.m1.pos);
      // score.m2.lft = idx.lift(score.m2.pos);

      orphan_paired_score_t score = paired_chain_orphan_score(al, i, mean, std_dev);

      // Check if there is another read scored in the same region
      if (score.tot >= al.min_score)
      {
        bool replaced = false;
        for(size_t j = 0; j < best_scores.size(); ++j)
        {
          orphan_paired_score_t zero;
          zero.chain_i = i;
          auto d1 = dist(best_scores[j].m1.lft,score.m1.lft);
          auto d2 = dist(best_scores[j].m2.lft,score.m2.lft);
          if ((dist(best_scores[j].m1.lft, score.m1.lft) < region_dist) and
              (dist(best_scores[j].m2.lft, score.m2.lft) < region_dist))
            if (score.tot > best_scores[j].tot)
              if (replaced)
                best_scores[j] = zero;
              else
                best_scores[j] = score, replaced = true;
            else if (score.tot <= best_scores[j].tot)
              j = best_scores.size(), replaced = true;
        }
        if (not replaced)
          best_scores.push_back(score);

      }

      ++i;
      // best_scores.push_back(std::make_pair(std::make_pair(score.tot,score.pos), i++));
    }

    orphan_paired_score_t zero;
    zero.chain_i = al.chains.size();

    while (best_scores.size() < 2)
      best_scores.push_back(zero);

    // Sort the chains by score
    std::sort(best_scores.begin(), best_scores.end(), std::greater<orphan_paired_score_t>());

    if (best_scores[0].tot < al.min_score)
    {
      MTIME_END(7); //Timing helper
      MTIME_TSAFE_MERGE;
      return false;
    }

    // Compute sub-optimal hits (https://github.com/bwa-mem2/bwa-mem2/blob/edc703f883e8aaed83067100d8e54e0e9e810ef5/src/bwamem.cpp#L1312)
    // From BWA's manual: BWA will not search for suboptimal hits with a score lower than (bestScore-misMsc).
    size_t j = 1;
    al.sub_n = 0;
    while (j < best_scores.size() and best_scores[j++].tot >= (best_scores[0].tot - max_penalty))
      ++al.sub_n;

    al.best_score = true;

    al.score2 = best_scores[1].tot;
    al.score2_m1 = best_scores[1].m1.score;
    al.score2_m2 = best_scores[1].m2.score;

    al.second_best_score = (al.score2 >= al.min_score);


    { // Forward case
      size_t i = best_scores[0].chain_i;
      ll start = best_scores[0].pos.first;
      ll end = best_scores[0].pos.second;
      // // Align the chain
      // auto chain = al.chains[i];
      // // Reverse the chain order
      // std::reverse(chain.anchors.begin(), chain.anchors.end());
      // // Compute the score of a chain.
      // // if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
      // //   al.score = paired_chain_orphan_score(chain, al.anchors, al.mems, al.min_score, al.min_score_m1, al.min_score_m2, al.mate1, &al.mate2_rev, mean, std_dev, false, al.score2, 0, start, end,  &al.sam_m1, &al.sam_m2);
      // // else
      // //   al.score = paired_chain_orphan_score(chain, al.anchors, al.mems, al.min_score, al.min_score_m1, al.min_score_m2, &al.mate1_rev, al.mate2, mean, std_dev, false, al.score2, 1, start, end, &al.sam_m1, &al.sam_m2);
      // al.score = paired_chain_orphan_score(chain, al, mean, std_dev, false, start, end);    
      al.score = paired_chain_orphan_score(al, i, mean, std_dev, false, start, end);    
    }

    al.aligned = (al.score.tot >= al.min_score);

    MTIME_END(7); //Timing helper
    MTIME_TSAFE_MERGE;

    return al.aligned;
  }

  // void find_mems(
  //   const kseq_t *read,
  //   std::vector<mem_t>& mems,
  //   size_t r_offset = 0,
  //   size_t mate = 0
  //   ) 
  // {
  //   auto pointers = ms.query(read->seq.s, read->seq.l);
  //   size_t l = 0;   // Current match length
  //   size_t pl = 0;  // Previous match length
  //   size_t n_Ns = 0;
  //   for (size_t i = 0; i < pointers.size(); ++i)
  //   {
  //     size_t pos = pointers[i];
  //     while ((i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
  //     {
  //       if(read->seq.s[i + l] == 'N') n_Ns++;
  //       else n_Ns = 0;
  //       ++l;
  //     }

  //     // Update MEMs
  //     if (l >= pl and n_Ns < l and l >= min_len)
  //     {
  //       size_t r = r_offset + (i + l - 1); // compatible with minimap2 chaining algorithm
  //       // size_t r = r_offset + ((reverse) ? (i) : (i + l - 1)); // compatible with minimap2 chaining algorithm
  //       // size_t r = i + l - 1;
  //       // r = (reverse)? (read->seq.l - (r + 1 - l) - 1) : (r); // compatible with minimap2 chaining algorithm
  //       mems.push_back(mem_t(pointers[i],l,i,mate,r));
  //       // find_MEM_occs(mems.back());
  //     }

  //     // Compute next match length
  //     pl = l;
  //     l = (l == 0 ? 0 : (l - 1));
  //   }

  // }

  double read_shredding_factor(
    std::vector<mem_t>& mems
    )
  {
    // for (size_t j = 0; j < n_MEMs; ++i)

  }

  // // Populate the seeds given a list of MEMs
  // void populate_seeds(
  //   std::vector<mem_t>& mems
  //   ) 
  // {
  //   size_t n_MEMs = mems.size();
  //   for (size_t j = 0; j < n_MEMs; ++j)
  //   {
  //     auto & mem = mems[j];
  //     size_t l = mem.len;
  //     size_t i = mem.idx;
  //     size_t mate = mem.mate;
  //     size_t pos = mem.pos;
  //     size_t r = mem.rpos;

  //     // size_t r = r_offset + (i + l - 1); // compatible with minimap2 chaining algorithm

  //     mem.occs.push_back(mem.pos);

  //     mem_finder.find_MEM_above(mem.pos, mem.len, mem.occs);
  //     size_t upper_suffix = mem.occs.back();
  //     mem_finder.find_MEM_below(mem.pos, mem.len, mem.occs);
  //     size_t lower_suffix = mem.occs.back();

  //     // Take two halves of the MEM
  //     if(l >= (min_len << 1))
  //     {
  //       size_t ll = l >> 1; 
  //       size_t rl = r - l + ll; // compatible with minimap2 chaining algorithm
  //       mems.push_back(mem_t(upper_suffix,ll,i,mate,rl));

  //       mem_t &mem = mems.back();
  //       if ((not mem_finder.find_MEM_above(upper_suffix, mem.len, mem.occs)) or
  //           (not mem_finder.find_MEM_below(lower_suffix, mem.len, mem.occs))){
  //           mems.pop_back(); continue;
  //           }

  //       size_t lr = l - ll;
  //       size_t rr = r; // compatible with minimap2 chaining algorithm
  //       mems.push_back(mem_t(pos + ll,lr,i + ll,mate,rr));
  //       mem_finder.find_MEM_occs(mems.back()); // TODO: Optimize this
  //       if ((not mem_finder.find_MEM_occs(mems.back())))
  //       {
  //         mems.pop_back();
  //         continue;
  //       }
  //     }
  //   }
  // }


  // void find_seeds(
  //   const kseq_t *read,
  //   std::vector<mem_t>& mems,
  //   size_t r_offset = 0,
  //   size_t mate = 0
  //   ) 
  // {
  //   auto pointers = ms.query(read->seq.s, read->seq.l);
  //   size_t l = 0;   // Current match length
  //   size_t pl = 0;  // Previous match length
  //   size_t n_Ns = 0;
  //   for (size_t i = 0; i < pointers.size(); ++i)
  //   {
  //     size_t pos = pointers[i];
  //     while ((i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
  //     {
  //       if(read->seq.s[i + l] == 'N') n_Ns++;
  //       else n_Ns = 0;
  //       ++l;
  //     }

  //     // Update MEMs
  //     if (l >= pl and n_Ns < l and l >= min_len)
  //     {
  //       size_t r = r_offset + (i + l - 1); // compatible with minimap2 chaining algorithm
  //       mems.push_back(mem_t(pointers[i],l,i,mate,r));

  //       mem_t& mem = mems.back();
  //       mem.occs.push_back(mem.pos);

  //       find_MEM_above(mem.pos, mem.len, mem.occs);
  //       size_t upper_suffix = mem.occs.back();
  //       find_MEM_below(mem.pos, mem.len, mem.occs);
  //       size_t lower_suffix = mem.occs.back();

  //       // Take two halves of the MEM
  //       if(l >= (min_len << 1))
  //       {
  //         size_t ll = l >> 1; 
  //         size_t rl = r_offset + (i + ll - 1); // compatible with minimap2 chaining algorithm
  //         mems.push_back(mem_t(upper_suffix,ll,i,mate,rl));

  //         mem_t &mem = mems.back();
  //         find_MEM_above(upper_suffix, mem.len, mem.occs);
  //         find_MEM_below(lower_suffix, mem.len, mem.occs);

  //         size_t lr = l - ll;
  //         size_t rr = r_offset + (i + l - 1); // compatible with minimap2 chaining algorithm
  //         mems.push_back(mem_t(pointers[i + ll],lr,i + ll,mate,rr));
  //         find_MEM_occs(mems.back()); // TODO: Optimize this
  //       }
  //     }

  //     // Compute next match length
  //     pl = l;
  //     l = (l == 0 ? 0 : (l - 1));
  //   }

  // }

  // void find_seeds_old(
  //   const kseq_t *read,
  //   std::vector<mem_t>& mems,
  //   size_t r_offset = 0,
  //   size_t mate = 0
  //   ) 
  // {
  //   auto pointers = ms.query(read->seq.s, read->seq.l);
  //   size_t l = 0;   // Current match length
  //   size_t pl = 0;  // Previous match length
  //   size_t n_Ns = 0;
  //   for (size_t i = 0; i < pointers.size(); ++i)
  //   {
  //     size_t pos = pointers[i];
  //     while ((i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
  //     {
  //       if(read->seq.s[i + l] == 'N') n_Ns++;
  //       else n_Ns = 0;
  //       ++l;
  //     }

  //     // Update MEMs
  //     if (l >= pl and n_Ns < l and l >= min_len)
  //     {
  //       size_t r = r_offset + (i + l - 1); // compatible with minimap2 chaining algorithm
  //       mems.push_back(mem_t(pointers[i],l,i,mate,r));
  //       find_MEM_occs(mems.back());
  //       // Take two halves of the MEM
  //       if(l >= (min_len << 1))
  //       {
  //         size_t ll = l >> 1; 
  //         size_t rl = r_offset + (i + ll - 1); // compatible with minimap2 chaining algorithm
  //         mems.push_back(mem_t(pointers[i],ll,i,mate,rl));
  //         find_MEM_occs(mems.back()); // TODO: Optimize this

  //         size_t lr = l - ll;
  //         size_t rr = r_offset + (i + l - 1); // compatible with minimap2 chaining algorithm
  //         mems.push_back(mem_t(pointers[i + ll],lr,i + ll,mate,rr));
  //         find_MEM_occs(mems.back()); // TODO: Optimize this
  //       }
  //     }

  //     // Compute next match length
  //     pl = l;
  //     l = (l == 0 ? 0 : (l - 1));
  //   }

  // }

  // Calculate fields for csv MEM stats file
  inline void calculate_MEM_stats(
    const std::vector<mem_t>& mems,
    csv_t& csv
  )
  {
    // Number of unique MEMs
    csv.num_uniq_mems = mems.size();
    // Total number of MEM occurances
    for ( size_t i = 0; i < mems.size(); ++i)
      csv.total_mem_occ += (mems[i].occs.size()+1); //+1 to account for the initial occ found in the ref with r-index
    // MEM frequency per pangenome and highest and lowest occurance of MEM on any one genome
    for ( size_t i = 0; i < mems.size(); ++i)
    {
      csv.max_mem_freq = (csv.max_mem_freq > (static_cast<double> (mems[i].occs.size() + 1) / csv.total_mem_occ) ? csv.max_mem_freq : (static_cast<double> (mems[i].occs.size() + 1) / csv.total_mem_occ));
      csv.min_mem_freq = (csv.min_mem_freq > (static_cast<double> (mems[i].occs.size() + 1) / csv.total_mem_occ) ? (static_cast<double> (mems[i].occs.size() + 1) / csv.total_mem_occ) : csv.min_mem_freq);
      std::map<std::string, size_t> count_dict;
      for (size_t j = 0; j < mems[i].occs.size(); ++j)
      {
        std::string ref = idx[mems[i].occs[j]];
        auto it = count_dict.find(ref);
        if (it != count_dict.end()) 
          count_dict[ref]++;
        else
          count_dict[ref] = 1;
      }
      // Iterate through key value pairs and find the high and lowest occurance of a MEM on any particular genome
      for (auto it = count_dict.begin(); it != count_dict.end(); ++it) 
      {
        if (csv.high_occ_mem == 0 && csv.low_occ_mem == 0)
        {
          csv.high_occ_mem = it->second;
          csv.low_occ_mem = it->second;
        }
        else
        {
          csv.high_occ_mem = (csv.high_occ_mem > it->second ? csv.high_occ_mem : it->second);
          csv.low_occ_mem = (csv.low_occ_mem > it->second ? it->second : csv.low_occ_mem);
        }
      }
    }
  }

  // Filter seeds by frequency of occurance
  inline void seed_freq_filter(
    std::vector<mem_t>& mems,
    const double freq,
    csv_t& csv
  )
  {
    size_t total_mem_occ = 1; // Set to 1 to avoid any potential divide by 0 operations
    std::vector<size_t> delete_ind;

    // Calculate the total number of occurance of MEMs
    for ( size_t i = 0; i < mems.size(); ++i)
      total_mem_occ += mems[i].occs.size() + 1; //+1 to account for the initial occ found in the ref with r-index

    // Find the indices of MEMs that have occurance freq greater than threshold freq
    for ( size_t i = 0; i < mems.size(); ++i)
    {
      if ( (static_cast<double> (mems[i].occs.size() + 1) / total_mem_occ) > freq)
      {
        delete_ind.push_back(i);
        csv.num_mems_filter += (mems[i].occs.size() + 1);
      }
    }

    // Reverse the delete indices to safely delete
    std::reverse(delete_ind.begin(), delete_ind.end());
    for (size_t idx: delete_ind)
      mems.erase(mems.begin() + idx);
  }

  // Filter seeds by number of occurances per ref
  inline void seed_occ_filter(
    std::vector<mem_t>& mems,
    const size_t n_seeds_thr,
    csv_t& csv
  )
  {
    for(size_t i = 0; i < mems.size(); ++i)
    {
      std::map<std::string, size_t> count_dict;
      std::vector<size_t> delete_ind;
      // Keep count of MEMs for each ref, keep track of indices of MEMs > threshold
      for (size_t j = 0; j < mems[i].occs.size(); ++j)
      {
        std::string ref = idx[mems[i].occs[j]];
        auto it = count_dict.find(ref);
        if (it != count_dict.end()) 
        {
          if (count_dict[ref] > n_seeds_thr)
          {
            delete_ind.push_back(j);
            csv.num_mems_filter += 1;
          }
          else
            count_dict[ref]++;
        }
        else
          count_dict[ref] = 1;
      }
      // Reverse the delete indices to safely delete
      std::reverse(delete_ind.begin(), delete_ind.end());
      for (size_t idx: delete_ind)
        mems[i].occs.erase(mems[i].occs.begin() + idx);
    }
  }

  // Compute the fraction of repetitive seeds
  // Inpired from https://github.com/lh3/bwa/blob/0747fcc09d96ff44ce555f0c258d0f9762c20611/bwamem.c#L291
  inline float compute_frac_rep(
    const std::vector<mem_t>& mems,
    const int len,
    const size_t mate = 0
    ) 
  {
    // TODO: find a right value for max_occ
    return 0.0; 


    size_t max_occ = -1;
    size_t begin = 0;
    size_t end = 0;
    size_t l_rep = 0;
    std::vector< std::tuple<size_t, size_t, size_t> > segments;
    segments.reserve(mems.size());
    // Collect the segments
    for ( size_t i = 0; i < mems.size(); ++i) 
    {
      const mem_t& mem = mems[i];
      if( mem.mate != mate )
        continue;
      segments.push_back( std::make_tuple(mem.idx, mem.idx + mem.len -1, mem.occs.size()) );
    }
    // Sort the segments by starting position
    std::sort( segments.begin(), segments.end() );
    // Compute the fraction of repetitive seeds
    for ( size_t i = 0; i < segments.size(); ++i )
    {    
      const size_t s_begin = std::get<0>(segments[i]);
      const size_t s_end = std::get<1>(segments[i]);
      const size_t occs = std::get<2>(segments[i]);
      if (occs <= max_occ) continue;
      if (s_begin > end)
      {
        l_rep += end - begin;
        begin = s_begin;
        end = s_end;
      } 
      else end = end > s_end? end : s_end;
    }
    l_rep += end - begin;
    return  (float)l_rep / len;
  }

  score_t chain_score(
      const std::pair<size_t, std::vector<size_t>> &chain,
      const std::vector<std::pair<size_t, size_t>> &anchors,
      const std::vector<mem_t> &mems,
      const int32_t min_score,
      const kseq_t *read,
      const bool score_only = true,
      const int32_t score2 = 0,
      const uint8_t strand = 0,
      FILE *out = nullptr,
      sam_t* sam = nullptr,           // The SAM information pointer
      const int32_t sub_n = 0,        // Number of sub-optimal alignments
      const double frac_rep = 0       // Length of the region covered by seeds (https://github.com/lh3/bwa/blob/0747fcc09d96ff44ce555f0c258d0f9762c20611/bwamem.c#L291)
      )     
  {
    bool aligned = false;
    // Extract the anchors
    std::vector<std::pair<size_t,size_t>> chain_anchors(chain.second.size());
    for(size_t i = 0; i < chain_anchors.size(); ++i)
      chain_anchors[i] = anchors[chain.second[i]];
    // Extracting left and right context of the read
    // lcs: left context sequence
    size_t lcs_len = mems[chain_anchors[0].first].idx;
    uint8_t* lcs = (uint8_t*)malloc(lcs_len);
    // verbose("lcs: " + std::string(read->seq.s).substr(0,lcs_len));
    // Convert A,C,G,T,N into 0,1,2,3,4
    // The left context is reversed
    for (size_t i = 0; i < lcs_len; ++i)
      lcs[lcs_len -i -1] = seq_nt4_table[(int)read->seq.s[i]];

    // rcs: right context sequence
    size_t rcs_occ = (mems[chain_anchors.back().first].idx + mems[chain_anchors.back().first].len); // The first character of the right context
    size_t rcs_len = read->seq.l - rcs_occ;
    uint8_t* rcs = (uint8_t*)malloc(rcs_len);
    // verbose("rcs: " + std::string(read->seq.s).substr(rcs_occ,rcs_len));
    // Convert A,C,G,T,N into 0,1,2,3,4
    for (size_t i = 0; i < rcs_len; ++i) 
      rcs[i] = seq_nt4_table[(int)read->seq.s[ rcs_occ + i]];

    // Fill between MEMs
    score_t score = fill_chain(
        mems,
        chain_anchors,
        lcs,   // Left context of the read
        lcs_len, // Left context of the read lngth
        rcs,   // Right context of the read
        rcs_len, // Right context of the read length
        read
      );

      if(!score_only and score.score >= min_score)
      {
        bool output = (sam == nullptr);

        if(output)
          sam = new sam_t(read);

        auto tmp_score = fill_chain(mems,chain_anchors,lcs,lcs_len,rcs,rcs_len,read,score_only,sam); 
        score.unmapped_lft = tmp_score.unmapped_lft;
        sam->flag = (strand?16:0);
        sam->zs = score2;
        // sam->mapq = compute_mapq(sam->as, sam->zs, min_score, sam->read->seq.l * smatch);
        // TODO: update the value of read->seq.l with the exact qlen
        sam->mapq = compute_mapq_se_bwa(sam->as, sam->zs, sam->rlen, read->seq.l, min_len, smatch, smismatch,
                                        mapq_coeff_len, mapq_coeff_fac, sub_n, 0, frac_rep);
        // sam->reverse = (strand != 0);
        // sam->rname = idx[sam->pos - 1];
        if(output && !score.unmapped_lft)
        {
          write_sam(out,*sam);
          free(sam);
        }
        
        // fill_chain(mems,chain_anchors,lcs,lcs_len,rcs,rcs_len,read,score_only,0,min_score,strand,out); 
      }

      free(rcs);
      free(lcs);

      return score;
  }

  // paired_score_t paired_chain_score(
  //     const chain_t &chain,
  //     const std::vector<std::pair<size_t, size_t>> &anchors,
  //     const std::vector<mem_t> &mems,
  //     const int32_t min_score,
  //     const int32_t min_score_m1,
  //     const int32_t min_score_m2,
  //     const kseq_t *mate1,
  //     const kseq_t *mate2,
  //     const bool score_only = true,
  //     const int32_t score2 = 0,
  //     const uint8_t strand = 0,
  //     sam_t *sam_m1_ = nullptr,
  //     sam_t *sam_m2_ = nullptr)
  // {
  paired_score_t paired_chain_score(
      // const chain_t &chain,
      paired_alignment_t& al,
      const size_t chain_i,
      const bool score_only = true)
  {
    const std::vector<std::pair<size_t, size_t>> &anchors = al.anchors;
    const std::vector<mem_t> &mems = al.mems;
    sam_t& sam_m1 = al.sam_m1;
    sam_t& sam_m2 = al.sam_m2;

    kseq_t *mate1;
    kseq_t *mate2;

    uint8_t strand = 0;

    // Get the chain
    auto& chain = al.chains[chain_i];
    // Reverse the chain order
    chain.reverse();


    if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
    {
      mate1 = al.mate1;
      mate2 = &al.mate2_rev;
    }
    else
    {
      mate1 = &al.mate1_rev;
      mate2 = al.mate2;
      strand = 1;
    }


    paired_score_t score;
    score.chain_i = chain_i;

    if(chain.paired){

      // Extract the anchors
      std::pair<ll, std::vector<size_t>> mate1_chain;
      std::pair<ll, std::vector<size_t>> mate2_chain;
      for(size_t i = 0; i < chain.anchors.size(); ++i)
      {
        size_t anchor_id = chain.anchors[i];
        if ((mems[anchors[anchor_id].first].mate & MATE_2) == 0)
          mate1_chain.second.push_back(anchor_id);
        else
          mate2_chain.second.push_back(anchor_id);
      }

      score.paired = chain.paired;

      if(score_only)
      {
          score.m1 = chain_score(mate1_chain, anchors, mems, al.min_score_m1, mate1);
          score.m2 = chain_score(mate2_chain, anchors, mems, al.min_score_m2, mate2);
          // score.dist = (score.m2.pos - (score.m1.pos + mate1->seq.l)) < 0 ? (score.m1.pos - (score.m2.pos + mate2->seq.l) ): (score.m2.pos - (score.m1.pos + mate1->seq.l));
          score.dist = dist(score.m2.pos,(score.m1.pos + mate1->seq.l));
          double ns = 0.0;
          if (al.std_dev > 0.0) ns = (score.dist - al.mean) / al.std_dev;
          score.tot = (int)(score.m1.score + score.m2.score + .721 * log(2. * erfc(fabs(ns) * M_SQRT1_2)) * smatch + .499);
          if (score.tot < 0)
            score.tot = 0;

          score.m1.lft = idx.lift(score.m1.pos);
          score.m2.lft = idx.lift(score.m2.pos);

      }
      else
      {
        score.m1 = chain_score(mate1_chain, anchors, mems, al.min_score_m1, mate1, false, al.score2_m1, strand, nullptr, &sam_m1, al.sub_n, al.frac_rep_m1);
        score.m2 = chain_score(mate2_chain, anchors, mems, al.min_score_m2, mate2, false, al.score2_m2, strand, nullptr, &sam_m2, al.sub_n, al.frac_rep_m2);
        // score.dist = (score.m2.pos - (score.m1.pos + mate1->seq.l)) < 0 ? (score.m1.pos - (score.m2.pos + mate2->seq.l)) : (score.m2.pos - (score.m1.pos + mate1->seq.l));
        score.dist = dist(score.m2.pos, (score.m1.pos + mate1->seq.l));
        double ns = 0.0;
        if (al.std_dev > 0.0)
          ns = (score.dist - al.mean) / al.std_dev;
        score.tot = (int)(score.m1.score + score.m2.score + .721 * log(2. * erfc(fabs(ns) * M_SQRT1_2)) * smatch + .499);
        if (score.tot < 0) score.tot = 0;

        score.m1.lft = idx.lift(score.m1.pos);
        score.m2.lft = idx.lift(score.m2.pos);
        
        sam_m1.read = mate1;
        sam_m2.read = mate2;
        // Fill sam fields RNEXT, PNEXT and TLEN
        // sam_m1.rnext = std::string(mate2->name.s);
        // sam_m2.rnext = std::string(mate1->name.s);

        if(score.m1.score >= al.min_score_m1 && !score.m1.unmapped_lft && score.m2.score >= al.min_score_m2 && !score.m2.unmapped_lft)
        {

          sam_m1.pnext = sam_m2.pos;
          sam_m2.pnext = sam_m1.pos;

          ll tlen;
          if (sam_m2.pos > sam_m1.pos)
          {
            tlen = (sam_m2.pos + mate2->seq.l) - sam_m1.pos;
            sam_m1.tlen = tlen;
            sam_m2.tlen = -tlen;
          }
          else
          {
            tlen = (sam_m1.pos + mate1->seq.l) - sam_m2.pos;
            sam_m1.tlen = -tlen;
            sam_m2.tlen = tlen;
          }
           

          //sam_m1.rname = idx[sam_m1.pos - 1]; // Check if necessary [This is causing the discrepency between OA Ref and Ref. Anyways REF already calculated by fill_chain function]
          //sam_m2.rname = idx[sam_m2.pos - 1];

          // size_t mapq = compute_mapq(score.tot, score2, al.min_score, (sam_m1.read->seq.l + sam_m2.read->seq.l) * smatch);
          // sam_m1.mapq = mapq;
          // sam_m2.mapq = mapq;

          int32_t score_un = 0; // score.m1 + score.m2 - unpaired_penalty;
          size_t mapq = compute_mapq_pe_bwa(score.tot, al.score2, score_un, smatch, al.sub_n, al.frac_rep_m1, al.frac_rep_m2, score.m1.score, score.m2.score, al.score2_m1, al.score2_m2, sam_m1.mapq, sam_m2.mapq);

          sam_m1.as = score.tot;
          sam_m2.as = score.tot;

          sam_m1.zs = al.score2;
          sam_m2.zs = al.score2;

          sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_MAPPED_PAIRED;
          if(strand)
          {
            sam_m1.flag |= SAM_REVERSED | SAM_FIRST_IN_PAIR;
            sam_m2.flag |= SAM_MATE_REVERSED | SAM_SECOND_IN_PAIR;
          }
          else
          {
            sam_m1.flag |= SAM_MATE_REVERSED | SAM_FIRST_IN_PAIR;
            sam_m2.flag |= SAM_REVERSED | SAM_SECOND_IN_PAIR;
          }
        }else if(score.m1.score >= al.min_score_m1 && !score.m1.unmapped_lft) {
          //sam_m1.rname = idx[sam_m1.pos - 1]; //Do not need since calculated in fill_chain function
          sam_m1.zs = al.score2_m1;
          
          sam_m1.flag = SAM_PAIRED | SAM_MATE_UNMAPPED | SAM_FIRST_IN_PAIR;
          sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_SECOND_IN_PAIR;
          if(strand)
            sam_m1.flag |= SAM_REVERSED;

          // Using GATK unmapped mate convention: https://gatk.broadinstitute.org/hc/en-us/articles/360035891031-Mate-unmapped-records
          sam_m2.rname = sam_m1.rname;
          sam_m2.pos = sam_m1.pos;
          sam_m2.mapq = sam_m1.mapq;
          sam_m2.cigar = "*";
          sam_m2.pnext = sam_m1.pnext = sam_m1.pos;
          sam_m2.tlen = sam_m1.tlen = 0;
        }else if(score.m2.score >= al.min_score_m2 && !score.m2.unmapped_lft) {
          //sam_m2.rname = idx[sam_m2.pos - 1]; //Do not need since calculated in fill_chain function
          sam_m1.zs = al.score2_m2;

          sam_m1.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_FIRST_IN_PAIR;
          sam_m2.flag = SAM_PAIRED | SAM_MATE_UNMAPPED | SAM_SECOND_IN_PAIR;
          if(not strand)
            sam_m2.flag |= SAM_REVERSED;

          // Using GATK unmapped mate convention: https://gatk.broadinstitute.org/hc/en-us/articles/360035891031-Mate-unmapped-records
          sam_m1.rname = sam_m2.rname;
          sam_m1.pos = sam_m2.pos;
          sam_m1.mapq = sam_m2.mapq;
          sam_m1.cigar = "*";
          sam_m1.pnext = sam_m2.pnext = sam_m2.pos;
          sam_m1.tlen = sam_m2.tlen = 0;
        }else {
          sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;
        }

        // write_sam(out,sam_m1);
        // write_sam(out,sam_m2);
      }
    }
    else
    {
      if(chain.mate == 0)
      {
        // score.m1 = chain_score(mate1_chain, anchors, mems, min_score, mate1, false, score2, strand, nullptr, &sam_m1);
      }
      else
      {
        // score.m2 = chain_score(mate2_chain, anchors, mems, min_score, mate2, false, score2, strand, nullptr, &sam_m2);

      }
    }

    return score;
  }



// orphan_paired_score_t paired_chain_orphan_score(
//       const chain_t &chain,
//       const std::vector<std::pair<size_t, size_t>> &anchors,
//       const std::vector<mem_t> &mems,
//       const int32_t min_score,
//       const int32_t min_score_m1,
//       const int32_t min_score_m2,
//       const kseq_t *mate1,
//       const kseq_t *mate2,
//       const double mean,
//       const double std_dev,
//       const bool score_only = true,
//       const int32_t score2 = 0,
//       const uint8_t strand = 0,
//       long long int start = 0,  // Start position of the reference
//       long long int end = 0,    //End position of the reference
//       sam_t *sam_m1_ = nullptr,
//       sam_t *sam_m2_ = nullptr)
orphan_paired_score_t paired_chain_orphan_score(
      // const chain_t &chain,
      paired_alignment_t& al,
      const size_t chain_i,
      const double mean,
      const double std_dev,
      const bool score_only = true,
      long long int start = 0,  // Start position of the reference
      long long int end = 0    //End position of the reference
    )
  {
    const std::vector<std::pair<size_t, size_t>> &anchors = al.anchors;
    const std::vector<mem_t> &mems = al.mems;
    sam_t& sam_m1 = al.sam_m1;
    sam_t& sam_m2 = al.sam_m2;

    kseq_t *mate1;
    kseq_t *mate2;

    uint8_t strand = 0;

    // Get the chain
    auto& chain = al.chains[chain_i];
    // Reverse the chain order
    chain.reverse();

    if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
    {
      mate1 = al.mate1;
      mate2 = &al.mate2_rev;
    }
    else
    {
      mate1 = &al.mate1_rev;
      mate2 = al.mate2;
      strand = 1;
    }

    orphan_paired_score_t score;
    score.chain_i = chain_i;

    // Extract the anchors
    std::pair<ll, std::vector<size_t>> mate1_chain;
    std::pair<ll, std::vector<size_t>> mate2_chain;
    size_t lm_pos = -1; // Leftmost position
    size_t rm_pos = 0;  // Rightmost position
    bool forward = true;
    for(size_t i = 0; i < chain.anchors.size(); ++i)
    {
      size_t anchor_id = chain.anchors[i];
      auto& mem = mems[anchors[anchor_id].first];
      maxl(rm_pos, mem.occs[anchors[anchor_id].second] + mem.len);
      minl(lm_pos, mem.occs[anchors[anchor_id].second]);
      if ((mem.mate & MATE_2) == 0)
        mate1_chain.second.push_back(anchor_id);
      else
        mate2_chain.second.push_back(anchor_id);
    }

    if(score_only)
    {
      if(mate1_chain.second.size() > 0)
      {
        // Chain 
        score.m1 = chain_score(mate1_chain, anchors, mems, al.min_score_m1, mate1);
        // Perform local search
        start = rm_pos + (ll)std::floor(mean - 4*std_dev);
        end = rm_pos + (ll)std::ceil(mean + 4*std_dev);
        maxl(start,(ll)0);
        minl(end,(ll)n);
        
        if (start < end)
          score.m2 = fill_orphan(start,end,mate2);  
        score.pos = std::pair(start,end);
      }else
      {
        score.m2 = chain_score(mate2_chain, anchors, mems, al.min_score_m2, mate2);
        // Perform local search
        start = lm_pos + (ll)std::floor(-mean - 4*std_dev);
        end = lm_pos + (ll)std::ceil(-mean + 4*std_dev);
        maxl(start,(ll)0);
        minl(end,(ll)n);

        if (start < end)
          score.m1 = fill_orphan(start,end,mate1);  
        score.pos = std::pair(start,end);
      }
      // score.dist = (score.m2.pos - (score.m1.pos + mate1->seq.l)) < 0 ? (score.m1.pos - (score.m2.pos + mate2->seq.l)) : (score.m2.pos - (score.m1.pos + mate1->seq.l));
      score.dist = dist(score.m2.pos, (score.m1.pos + mate1->seq.l));
      double ns = 0.0;
      if (al.std_dev > 0.0)
        ns = (score.dist - al.mean) / al.std_dev;
      score.tot = (int)(score.m1.score + score.m2.score + .721 * log(2. * erfc(fabs(ns) * M_SQRT1_2)) * smatch + .499);
      if (score.tot < 0)
        score.tot = 0;

      score.m1.lft = idx.lift(score.m1.pos);
      score.m2.lft = idx.lift(score.m2.pos);
    }
    else
    {
      // sam_t& sam_m1 = *sam_m1_;
      // sam_t& sam_m2 = *sam_m2_;

      if(mate1_chain.second.size() > 0)
      {
        score.m1 = chain_score(mate1_chain, anchors, mems, al.min_score_m1, mate1, false, al.score2_m1, strand, nullptr, &sam_m1, al.sub_n, al.frac_rep_m1);
        if (start < end)
          score.m2 = fill_orphan(start,end,mate2,false,&sam_m2);

        sam_m2.mapq = compute_mapq_se_bwa(sam_m2.as, al.score2_m2, sam_m2.rlen, mate2->seq.l, min_len, smatch, smismatch,
                                mapq_coeff_len, mapq_coeff_fac, al.sub_n, 0, al.frac_rep_m2);
      }
      else
      {
        if (start < end)
          score.m1 = fill_orphan(start,end,mate1,false,&sam_m1);
        score.m2 = chain_score(mate2_chain, anchors, mems, al.min_score_m2, mate2, false, al.score2_m2, strand, nullptr, &sam_m2, al.sub_n, al.frac_rep_m2);

        sam_m1.mapq = compute_mapq_se_bwa(sam_m1.as, al.score2_m1, sam_m1.rlen, mate1->seq.l, min_len, smatch, smismatch,
                        mapq_coeff_len, mapq_coeff_fac, al.sub_n, 0, al.frac_rep_m1);
      }
      // score.dist = (score.m2.pos - (score.m1.pos + mate1->seq.l)) < 0 ? (score.m1.pos - (score.m2.pos + mate2->seq.l)) : (score.m2.pos - (score.m1.pos + mate1->seq.l));
      score.dist = dist(score.m2.pos, (score.m1.pos + mate1->seq.l));
      double ns = 0.0;
      if (al.std_dev > 0.0)
        ns = (score.dist - al.mean) / al.std_dev;
      score.tot = (int)(score.m1.score + score.m2.score + .721 * log(2. * erfc(fabs(ns) * M_SQRT1_2)) * smatch + .499);
      if (score.tot < 0) score.tot = 0;

      score.m1.lft = idx.lift(score.m1.pos);
      score.m2.lft = idx.lift(score.m2.pos);

      sam_m1.read = mate1;
      sam_m2.read = mate2;
      // Fill sam fields RNEXT, PNEXT and TLEN
      // sam_m1.rnext = std::string(mate2->name.s);
      // sam_m2.rnext = std::string(mate1->name.s);

      if(score.m1.score >= al.min_score_m1 && !score.m1.unmapped_lft and score.m2.score >= al.min_score_m2 && !score.m2.unmapped_lft)
      {

        sam_m1.pnext = sam_m2.pos;
        sam_m2.pnext = sam_m1.pos;

        ll tlen;
        if (sam_m2.pos > sam_m1.pos)
        {
          tlen = (sam_m2.pos + mate2->seq.l) - sam_m1.pos;
          sam_m1.tlen = tlen;
          sam_m2.tlen = -tlen;
        }
        else
        {
          tlen = (sam_m1.pos + mate1->seq.l) - sam_m2.pos;
          sam_m1.tlen = -tlen;
          sam_m2.tlen = tlen;
        }


        // sam_m1.rname = idx[sam_m1.pos - 1];
        // sam_m2.rname = idx[sam_m2.pos - 1];
        // size_t mapq = compute_mapq(score.tot, score2, min_score, (sam_m1.read->seq.l + sam_m2.read->seq.l) * smatch);
        
        // sam_m1.mapq = mapq;
        // sam_m2.mapq = mapq;

        int32_t score_un = 0; // score.m1 + score.m2 - unpaired_penalty;
        size_t mapq = compute_mapq_pe_bwa(score.tot, al.score2, score_un, smatch, al.sub_n, al.frac_rep_m1, al.frac_rep_m2, score.m1.score, score.m2.score, al.score2_m1, al.score2_m2, sam_m1.mapq, sam_m2.mapq);

        sam_m1.as = score.tot;
        sam_m2.as = score.tot;

        sam_m1.zs = al.score2;
        sam_m2.zs = al.score2;

        sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_MAPPED_PAIRED;
        if(strand)
        {
          sam_m1.flag |= SAM_REVERSED | SAM_FIRST_IN_PAIR;
          sam_m2.flag |= SAM_MATE_REVERSED | SAM_SECOND_IN_PAIR;
        }
        else
        {
          sam_m1.flag |= SAM_MATE_REVERSED | SAM_FIRST_IN_PAIR;
          sam_m2.flag |= SAM_REVERSED | SAM_SECOND_IN_PAIR;
        }
      }else if(score.m1.score >= al.min_score_m1 && !score.m1.unmapped_lft) {
        // sam_m1.rname = idx[sam_m1.pos - 1];
        sam_m1.zs = al.score2_m1;

        sam_m1.flag = SAM_PAIRED | SAM_MATE_UNMAPPED | SAM_FIRST_IN_PAIR;
        sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_SECOND_IN_PAIR;
        if(strand)
          sam_m1.flag |= SAM_REVERSED;

        // Using GATK unmapped mate convention: https://gatk.broadinstitute.org/hc/en-us/articles/360035891031-Mate-unmapped-records
        sam_m2.rname   = sam_m1.rname;
        sam_m2.pos   = sam_m1.pos;
        sam_m2.mapq  = sam_m1.mapq;
        sam_m2.cigar = "*";
        sam_m2.pnext = sam_m1.pnext = sam_m1.pos;
        sam_m2.tlen  = sam_m1.tlen = 0;

      }else if(score.m2.score >= al.min_score_m2 && !score.m2.unmapped_lft) {
        // sam_m2.rname = idx[sam_m2.pos - 1];
        sam_m1.zs = al.score2_m2;

        sam_m1.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_FIRST_IN_PAIR;
        sam_m2.flag = SAM_PAIRED | SAM_MATE_UNMAPPED | SAM_SECOND_IN_PAIR;
        if(not strand)
          sam_m2.flag |= SAM_REVERSED;

        // Using GATK unmapped mate convention: https://gatk.broadinstitute.org/hc/en-us/articles/360035891031-Mate-unmapped-records
        sam_m1.rname = sam_m2.rname;
        sam_m1.pos = sam_m2.pos;
        sam_m1.mapq = sam_m2.mapq;
        sam_m1.cigar = "*";
        sam_m1.pnext = sam_m2.pnext = sam_m2.pos;
        sam_m1.tlen = sam_m2.tlen = 0;
      }else {
        sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;
      }
    }



    return score;
  }

  // If score_only is true we compute the score of the alignment. 
  // If score_only is false, we extend again the read and we write the result
  // in the SAM file, so we need to give the second best score. 
  // Return the start and end of the alignment found.
  score_t fill_orphan(
    long long int& start, // Start position of the reference
    long long int& end, //End position of the reference
    const kseq_t *paired_mate, // The read that has been aligned
    const bool score_only = true, // Report only the score
    sam_t* sam = nullptr,     // The SAM information pointer
    bool realign = false   // Realign globally the read
  )
  {
    score_t score = {0,0};
    // Extract reference
    size_t ref_occ = start;
    size_t ref_len = end - start + 1;
    char *ref = (char *)malloc(ref_len);
    if(ref_len <= 0)
      assert(ref_len > 0);
    ra.expandSubstr(ref_occ, ref_len, ref);
    // Convert A,C,G,T,N into 0,1,2,3,4
    for (size_t i = 0; i < ref_len; ++i)
      ref[i] = seq_nt4_table[(int)ref[i]];
    
    // Convert paired_mate in nt4 format
    size_t seq_len = paired_mate->seq.l;
    uint8_t *seq = (uint8_t *)malloc(paired_mate->seq.l);
    // Convert A,C,G,T,N into 0,1,2,3,4
    for (size_t i = 0; i < (int)paired_mate->seq.l; ++i)
        seq[i] = seq_nt4_table[(int)paired_mate->seq.s[i]];

    if(score_only or realign)
    { // TODO: Fix the double alignment
      kswq_t *q = 0;
      kswr_t r;
      const int xtra = KSW_XSTART;

      r = ksw_align(seq_len, (uint8_t *)seq, ref_len, (uint8_t *)ref, 5, mat, gapo, gape, xtra, &q);

      // Update start and end
      end = start + r.te;
      start += r.tb;

      int flag = KSW_EZ_SCORE_ONLY;
      size_t ref_len_ = r.te - r.tb + 1;
      char *ref_ = ref + r.tb;
      ksw_extz_t ez;
      memset(&ez, 0, sizeof(ksw_extz_t));

      ksw_reset_extz(&ez);
      // ksw_extd2_sse(km, seq_len, (uint8_t *)seq, ref_len, (uint8_t*)ref, m, mat, gapo, gape, gapo2, gape2, w, zdrop, end_bonus, flag, &ez);
      ksw_extz2_sse(km, seq_len, (uint8_t *)seq, ref_len_, (uint8_t*)ref_, m, mat, gapo, gape,  w, zdrop, end_bonus, flag, &ez);

      score.score= ez.score;
      score.pos = start;

      bool is_valid = idx.valid(start, end-start+1);
      if (not is_valid)
        score.score = std::numeric_limits<int32_t>::min();
      free(q);
    }


    if(not score_only)
    {
        // Realign the whole sequence globally
        int flag = KSW_EZ_RIGHT;

        ksw_extz_t ez;
        memset(&ez, 0, sizeof(ksw_extz_t));

        // char* tmp = (char*)calloc(max(ref_len,seq_len),1);

        ksw_reset_extz(&ez);
        // ksw_extd2_sse(km, seq_len, (uint8_t *)seq, ref_len, (uint8_t*)ref, m, mat, gapo, gape, gapo2, gape2, w, zdrop, end_bonus, flag, &ez);
        ksw_extz2_sse(km, seq_len, (uint8_t *)seq, ref_len, (uint8_t*)ref, m, mat, gapo, gape,  w, zdrop, end_bonus, flag, &ez);

        // ez.score = ksw_global(seq_len, (uint8_t *)seq, ref_len, (uint8_t *)ref, 5, mat, gapo, gape, seq_len, &ez.n_cigar, &ez.cigar);
        // Concatenate the CIGAR strings

        // std::string cigar_s;
        sam->lift_cigar = "";
        for (size_t i = 0; i < ez.n_cigar; ++i)
          sam->lift_cigar += std::to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];

        // sam->n_cigar = ez.n_cigar;
        // sam->cigar_b = (uint32_t*) malloc(ez.n_cigar * sizeof(uint32_t));
        // std::memcopy(sam->cigar_b, ez.cigar, ez.n_cigar * sizeof(uint32_t));

        // Compute the MD:Z field and the number of mismatches
        sam->lift_nm = write_MD_core((uint8_t *)ref, seq, ez.cigar, ez.n_cigar, nullptr, 0, sam->lift_md);

        const auto ref = idx.index(ref_occ);
        sam->as = ez.score;
        sam->lift_pos = ref.second + 1; // ref_occ + 1; // ref_occ is 1 based
        sam->lift_rname = ref.first;    // idx[ref_occ];
        sam->lift_rlen = ref_len;

        bam1_t* bam = bam_init1();
        bam_set1(bam, paired_mate->name.l, paired_mate->name.s, sam->flag, 0, ref.second, 0, ez.n_cigar, ez.cigar, 0, 0, 0, 0, NULL, NULL, 0);
        idx.lift_cigar(bam, ref_occ);


        const auto lift = idx.lift(ref_occ);
        const auto lft_ref = idx.index(lift);
        sam->pos = lft_ref.second + 1; // ref_occ + 1; // ref_occ is 1 based
        sam->rname = lft_ref.first;    // idx[ref_occ];

        uint32_t *cigar = bam_get_cigar(bam);
        sam->cigar = "";
        for (size_t i = 0; i < bam->core.n_cigar; ++i)
          sam->cigar += std::to_string(cigar[i] >> 4) + "MID"[cigar[i] & 0xf];

        ref_occ = lift; // This is correct since it is the position in the concatenation.
        ref_len = bam_cigar2rlen(bam->core.n_cigar, cigar);
        if(ref_len > 0)
        {
          char * l_ref = (char *)malloc( ref_len );
          assert(ref_len > 0);
          ra.expandSubstr(ref_occ, ref_len, l_ref);

          // Convert A,C,G,T,N into 0,1,2,3,4
          for (size_t i = 0; i < ref_len; ++i)
            l_ref[i] = seq_nt4_table[(int)l_ref[i]];
          // Compute the MD:Z field and the number of mismatches
          sam->nm = write_MD_core((uint8_t *)l_ref, seq, cigar, bam->core.n_cigar, nullptr, 0, sam->md);
          sam->rlen = ref_len;

          score.score = ez.score;
          score.pos = ref_occ;
          // score.pos = start;
          free(l_ref);
        }  else { // Read is unmapped because it align on an insertion of length > readlength
          sam->pos = 0;
          sam->rname = "*";
          sam->cigar = "*";
          sam->rlen = 0;
          sam->unmapped_lft = true;
          score.unmapped_lft = true;
        }
        // free(tmp);
        bam_destroy1(bam);
    }

    free(ref);
    free(seq);
    return score;  
  }

  size_t get_aligned_reads()
  {
    return aligned_reads;
  }

  typedef struct fill_chain_state_t {
    ksw_extz_t ez_lc;
    ksw_extz_t ez_rc;
    ksw_extz_t ez;
    std::vector<ksw_extz_t> ez_cc;

    fill_chain_state_t(size_t cc_size):
      ez_cc(cc_size)
    {
      memset(&ez_lc, 0, sizeof(ksw_extz_t));
      memset(&ez_rc, 0, sizeof(ksw_extz_t));
      memset(&ez, 0, sizeof(ksw_extz_t));
      for(size_t i = 0; i < cc_size; ++i)
        memset(&ez_cc[i], 0, sizeof(ksw_extz_t));
    }

    void del(){
      if (ez_lc.m_cigar > 0)
        free(ez_lc.cigar);
      if (ez_rc.m_cigar > 0)
        free(ez_rc.cigar);
      if (ez.m_cigar > 0)
        free(ez.cigar);
      for( size_t i = 0; i < ez_cc.size(); ++i )
        if(ez_cc[i].n_cigar > 0)
          free(ez_cc[i].cigar);
    }


  } fill_chain_state_t;

  // If score_only is true we compute the score of the alignment. 
  // If score_only is false, we extend again the read and we write the result
  // in the SAM file, so we need to give the second best score. 
  score_t fill_chain(
    const std::vector<mem_t>& mems, // pairs of pos/lengths of the mems
    const std::vector<std::pair<size_t,size_t>>& anchors, // pairs of mem index/occ index
    const uint8_t* lcs,   // Left context of the read
    const size_t lcs_len, // Left context of the read lngth
    const uint8_t* rcs,   // Right context of the read
    const size_t rcs_len, // Right context of the read length
    const kseq_t *read, // The read that has been aligned
    const bool score_only = true, // Report only the score
    sam_t* sam = nullptr,     // The SAM information pointer
    // const int32_t score2 = 0,    // The score of the second best alignment
    // const int32_t min_score = 0, // The minimum score to call an alignment
    // const int8_t strand = 0,    // 0: forward aligned ; 1: reverse complement aligned
    // const FILE *out = nullptr,   // The SAM file pointer
    bool realign = false   // Realign globally the read
    // const bool realign = false   // Realign globally the read
  )
  {
    int flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT;

    if(score_only) 
      flag = KSW_EZ_SCORE_ONLY;
    
    score_t score;

    int score_lc = 0;
    int score_rc = 0;

    ksw_extz_t ez_lc;
    ksw_extz_t ez_rc;
    ksw_extz_t ez;
    memset(&ez_lc, 0, sizeof(ksw_extz_t));
    memset(&ez_rc, 0, sizeof(ksw_extz_t));
    memset(&ez, 0, sizeof(ksw_extz_t));
    // TODO: Update end_bonus according to the MEM contribution to the score
    
    // Extract the context from the reference
    // lc: left context of the first mem
    ksw_reset_extz(&ez_lc);
    if(lcs_len > 0)
    {
      size_t mem_pos = mems[anchors[0].first].occs[anchors[0].second];

      size_t lc_occ = (mem_pos > ext_len ? mem_pos - ext_len : 0);
      size_t lc_len = (mem_pos > ext_len ? ext_len : ext_len - mem_pos);
      char *tmp_lc = (char *)malloc(ext_len);
      assert(lc_len > 0);
      ra.expandSubstr(lc_occ, lc_len, tmp_lc);
      // verbose("lc: " + std::string(lc));
      // Convert A,C,G,T,N into 0,1,2,3,4
      // The left context is reversed
      uint8_t *lc = (uint8_t *)malloc(ext_len);
      for (size_t i = 0; i < lc_len; ++i)
        lc[lc_len -i -1] = seq_nt4_table[(int)tmp_lc[i]];
      
      free(tmp_lc);

      // Query: lcs
      // Target: lc
      // verbose("aligning lc and lcs");
      ksw_extz2_sse(km, lcs_len, (uint8_t*)lcs, lc_len, (uint8_t*)lc, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez_lc);
      score_lc =  ez_lc.mqe;
      // verbose("lc score: " + std::to_string(score_lc));
      // Check if the extension reached the end of the query
      assert(score_only or ez_lc.reach_end);

      // std::string blc = print_BLAST_like((uint8_t*)lc,(uint8_t*)lcs,ez_lc.cigar,ez_lc.n_cigar);
      // std::cout<<blc;

      free(lc);
    }

    // rc: right context of the last mem
    ksw_reset_extz(&ez_rc);
    if(rcs_len > 0)
    {
      size_t mem_pos =  mems[anchors.back().first].occs[anchors.back().second];
      size_t mem_len =  mems[anchors.back().first].len;

      size_t rc_occ = mem_pos + mem_len;
      size_t rc_len = (rc_occ < n-ext_len ? ext_len : n - rc_occ);
      char *rc = (char *)malloc(ext_len);
      assert(rc_len > 0);
      ra.expandSubstr(rc_occ, rc_len, rc);
      // verbose("rc: " + std::string(rc));
      // Convert A,C,G,T,N into 0,1,2,3,4
      for (size_t i = 0; i < rc_len; ++i)
        rc[i] = seq_nt4_table[(int)rc[i]];

      // Query: rcs
      // Target: rc
      // verbose("aligning rc and rcs");
      ksw_extz2_sse(km, rcs_len, (uint8_t*)rcs, rc_len, (uint8_t*)rc, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez_rc);
      score_rc = ez_rc.mqe;
      // verbose("rc score: " + std::to_string(score_rc));
      // Check if the extension reached the end of the query
      assert(score_only or ez_rc.reach_end);

      // std::string brc = print_BLAST_like((uint8_t*)rc,(uint8_t*)rcs,ez_rc.cigar,ez_rc.n_cigar);
      // std::cout<<brc;
      free(rc);
    }


    // Compute the partial score score
    score.score = score_lc + score_rc;
    // int32_t score = mem_len * smatch + score_lc + score_rc;

    // Compute starting position in reference
    size_t mem_pos = mems[anchors[0].first].occs[anchors[0].second];
    size_t mem_len = mems[anchors.back().first].occs[anchors.back().second] + mems[anchors.back().first].len - mem_pos; // from the strart of the first MEM to the end of the last MEM.
    size_t ref_pos;
    // If the condidition is true, it results in the read going unmapped from my testing. 
    // Acceptable, but ideally clip the left end of the read.
    if ((lcs_len > 0 ? ez_lc.mqe_t + 1 : 0) > mem_pos)
      ref_pos = 0;
    else
      ref_pos = mem_pos - (lcs_len > 0 ? ez_lc.mqe_t + 1 : 0);
    size_t ref_len = (lcs_len > 0 ? ez_lc.mqe_t + 1 : 0) + mem_len + (rcs_len > 0 ? ez_rc.mqe_t + 1: 0);
    char *ref = (char *)malloc(ref_len);
    assert(ref_len > 0);
    ra.expandSubstr(ref_pos, ref_len, ref);

    // Convert A,C,G,T,N into 0,1,2,3,4
    for (size_t i = 0; i < ref_len; ++i)
      ref[i] = seq_nt4_table[(int)ref[i]];

    // Convert the read
    size_t seq_len = read->seq.l;
    uint8_t* seq = (uint8_t*) malloc(seq_len);
    for (size_t i = 0; i < seq_len; ++i)
      seq[i] = seq_nt4_table[(int)read->seq.s[i]];

    // Store the starting position of the allignment
    score.pos = ref_pos;

    bool mems_overlap = false;
    size_t last_ref = mem_pos + mems[anchors[0].first].len;
    size_t last_seq = mems[anchors[0].first].idx + mems[anchors[0].first].len;
    for(size_t i = 1; i < anchors.size() and not mems_overlap; ++i)
    {
      const size_t& ref_occ  = mems[anchors[i].first].occs[anchors[i].second];
      const size_t& seq_occ = mems[anchors[i].first].idx;
      const size_t& mem_len = mems[anchors[i].first].len;
      // Check with > because last is the first position not covered
      if (last_ref > ref_occ or last_seq > seq_occ) mems_overlap = true;
      last_ref = ref_occ + mem_len;
      last_seq = seq_occ + mem_len;
    }

    // TODO: fix the gap fill when MEMs overlap. The quick fix is a full alignment between the first and last MEM
    // Fill the gaps between each mem
    std::vector<ksw_extz_t> ez_cc(anchors.size()-1);
    if (not mems_overlap and not realign)
    {
      // Check if we need to extend between anchors or not.
      size_t last_ref = mem_pos + mems[anchors[0].first].len;
      size_t last_seq = mems[anchors[0].first].idx + mems[anchors[0].first].len;

      for(size_t i = 1; i < anchors.size(); ++i)
      {

        const size_t& ref_occ  = mems[anchors[i].first].occs[anchors[i].second];
        const size_t& seq_occ = mems[anchors[i].first].idx;
        const size_t& mem_len = mems[anchors[i].first].len;
        if(last_ref == ref_occ)
        {
          if(last_seq < seq_occ)
          {
            size_t l = (seq_occ - last_seq);
            ez_cc[i-1].score = -std::min(gapo + l*gape, gapo2 + l*gape2);
            ez_cc[i-1].m_cigar = 1;
            ez_cc[i-1].n_cigar = 1;
            ez_cc[i-1].cigar = (uint32_t*) malloc(ez_cc[i-1].m_cigar * sizeof(uint32_t));
            ez_cc[i-1].cigar[0] = (l << 4) |  1; // Insertion to ref
          }
          else
          {
            ez_cc[i-1].score = 0;
            ez_cc[i-1].m_cigar = 0;
            ez_cc[i-1].n_cigar = 0;
          }
        }
        else // last_ref < seq_occ
        {
          if(last_seq == seq_occ)
          {
            size_t l = (seq_occ - last_seq);
            ez_cc[i-1].score = -std::min(gapo + l*gape, gapo2 + l*gape2);
            ez_cc[i-1].m_cigar = 1;
            ez_cc[i-1].n_cigar = 1;
            ez_cc[i-1].cigar = (uint32_t*) malloc(ez_cc[i-1].m_cigar * sizeof(uint32_t));
            ez_cc[i-1].cigar[0] = (l << 4) |  2; // Deletion from ref
          }
          else
          {
            flag = KSW_EZ_RIGHT;
            ksw_reset_extz(&ez_cc[i-1]);

            size_t cc_occ = mems[anchors[i-1].first].occs[anchors[i-1].second] + mems[anchors[i-1].first].len;
            size_t cc_len = mems[anchors[i].first].occs[anchors[i].second] - cc_occ; 
            cc_occ -= ref_pos;
            // char *cc = (char *)malloc(cc_len);
            // // Convert A,C,G,T,N into 0,1,2,3,4
            // for (size_t i = 0; i < cc_len; ++i)
            //   cc[i] = ref[cc_occ + i];


            size_t ccs_pos = mems[anchors[i-1].first].idx + mems[anchors[i-1].first].len;
            size_t ccs_len = mems[anchors[i].first].idx - ccs_pos;

            // Query: rcs
            // Target: rc
            ksw_extz2_sse(km, ccs_len, (uint8_t*)(seq + ccs_pos), cc_len, (uint8_t*)(ref + cc_occ), m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez_cc[i-1]);

          }
        }
        last_ref = ref_occ + mem_len;
        last_seq = seq_occ + mem_len;

        // Update score
        // Here we use the score because we need to take into account eventual deletions
        // at the end of the reference.
        score.score +=  mems[anchors[i-1].first].len * smatch + ez_cc[i-1].score;

      }

      score.score +=  mems[anchors.back().first].len * smatch;
    }
    else
    {
  //******************************************************************************
  // BEGIN RAPID PROTOTYPING HACK
  //******************************************************************************
      // Realign the whole sequence globally
      ksw_reset_extz(&ez);
      ksw_extz2_sse(km, seq_len, (uint8_t *)seq, ref_len, (uint8_t *)ref, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez);

      score.score = ez.score;
      // TODO: Check if we have to update also the position
      realign = true;

  //******************************************************************************
  // END RAPID PROTOTYPING HACK
  //******************************************************************************
    }
    bool is_valid = idx.valid(ref_pos, ref_len);
    if (not is_valid) score.score = std::numeric_limits<int32_t>::min();

    if(is_valid and not score_only)
    {
      // Compute starting position in reference

      // char* tmp = (char*)calloc(max(ref_len,seq_len),1);

      size_t n_cigar = 0;
      uint32_t *cigar = nullptr;

      if(realign)
      {
        // Realign the whole sequence globally
        flag = KSW_EZ_RIGHT;
        ksw_reset_extz(&ez);
        ksw_extz2_sse(km, seq_len, (uint8_t*)seq, ref_len, (uint8_t*)ref, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez);


        // std::string bfull = print_BLAST_like((uint8_t*)ref,seq,ez.cigar,ez.n_cigar);
        // std::cout << bfull;

        // Example were ez.score is lrger than score:
        // Left context alignment
        // 22333022233022233302223302223
        // ||||  ||||||||| ||||||*|||||*
        // 2233  222330222 3302220302220
        // Right context alignment
        // 33022233022233022233022233022233      0222334
        // *||||||||||||||||*||||||||||||||      ||||||*
        // 130222330222330220330222330222332222330222330
        // [INFO] 16:26:16 - Message: old score:  130  new score:  140
        // Global alignment
        // 2233    3022233022233302223  30222330222330222330222330222  330222  330222330222330222330222330222330222334
        // ||||    |||||||||||*| ||||*  |||||||||||||||||||||||||||||  *|||||  |||||||||||*||||||||||||||*|||||||||||*
        // 223322233022233022203 02220  30222330222330222330222330222  130222  330222330220330222330222332222330222330
        // The original occurrence of the MEM has been shifted to the left by 6 positions, 
        // reducing the gap in the right context, and moving in to the left context.

        assert(ez.score >= score.score);

        // Concatenate the CIGAR strings
        n_cigar = ez.n_cigar;
        cigar = ez.cigar;
        score.score = ez.score;
        ez.m_cigar = 0; // To avoid double free
      }
      else
      {
        // Concatenate the CIGAR strings
        n_cigar = ez_lc.n_cigar + ez_rc.n_cigar + 1;
        for(size_t i = 0; i < anchors.size()-1; ++i)
          n_cigar += ez_cc[i].n_cigar + 1;

        cigar = (uint32_t*)calloc(n_cigar,sizeof(uint32_t));
        size_t i = 0;

        for(size_t j = 0; j < ez_lc.n_cigar; ++j)
          cigar[i++] = ez_lc.cigar[ez_lc.n_cigar -j -1];

        // CIGARs between MEMs
        for(size_t j = 0; j < anchors.size(); ++j)
        {
          // CIGAR for the MEM
          const size_t& mem_len = mems[anchors[j].first].len;
          if(i > 0 and ((cigar[i-1]& 0xf) == 0))
          { // If the previous operation is also an M then merge the two operations
            cigar[i-1] += (((uint32_t)mem_len) << 4);
            --n_cigar;
          }
          else
            cigar[i++] = (((uint32_t)mem_len) << 4);

          // CIGAR of the gap between the j-th MEM and the j+1-th MEM
          if(j < anchors.size() - 1)
          {
            if(ez_cc[j].n_cigar > 0)
            {
              // Check the first element
              if((ez_cc[j].cigar[0]& 0xf) == 0)
              { // If the next operation is also an M then merge the two operations
                cigar[i-1] += ez_cc[j].cigar[0];
                --n_cigar;
              }
              else
                cigar[i++] = ez_cc[j].cigar[0];
            } 
                      
              // Copy all the other elements
              for(size_t k = 1; k < ez_cc[j].n_cigar; ++k)
                cigar[i++] = ez_cc[j].cigar[k];
          }
        }

        // CIGAR of the right context
        if(ez_rc.n_cigar > 0)
        {
          if((ez_rc.cigar[0]& 0xf) == 0)
          { // If the next operation is also an M then merge the two operations
            cigar[i-1] += ez_rc.cigar[0];
            --n_cigar;
          }
          else
            cigar[i++] = ez_rc.cigar[0];
        }

        for(size_t j = 1; j < ez_rc.n_cigar; ++j)
          cigar[i++] = ez_rc.cigar[j];
        
        assert(i <= n_cigar);

      }




      // report output
      sam->lift_cigar = "";
      for(size_t i = 0; i < n_cigar; ++i)
        sam->lift_cigar += std::to_string(cigar[i] >> 4) + "MID"[cigar[i] & 0xf];

      // sam->n_cigar = n_cigar;
      // sam->cigar_b = (uint32_t*) malloc(n_cigar * sizeof(uint32_t));
      // std::memcopy(sam->cigar_b, cigar, n_cigar * sizeof(uint32_t));

      // Compute the MD:Z field and the number of mismatches
      sam->lift_nm = write_MD_core((uint8_t*)ref,seq,cigar,n_cigar,nullptr,0,sam->lift_md);

      const auto ref = idx.index(ref_pos);
      sam->as = score.score;
      sam->lift_pos = ref.second + 1; //ref_pos + 1; // ref_pos is 1 based
      sam->lift_rname = ref.first; // idx[ref_pos];
      sam->lift_rlen = ref_len;

      bam1_t* bam = bam_init1();
      bam_set1(bam, read->name.l, read->name.s, sam->flag, 0, ref.second, 0, n_cigar, cigar, 0, 0, 0, 0, NULL, NULL, 0);
      idx.lift_cigar(bam, ref_pos);


      const auto lift = idx.lift(ref_pos);
      const auto lft_ref = idx.index(lift);
      sam->pos = lft_ref.second + 1; //ref_occ + 1; // ref_occ is 1 based
      sam->rname = lft_ref.first; // idx[ref_occ];

      uint32_t *lft_cigar = bam_get_cigar(bam); // This has not to be deleted later
      sam->cigar = "";
      for (size_t i = 0; i < bam->core.n_cigar; ++i)
        sam->cigar += std::to_string(lft_cigar[i] >> 4) + "MID"[lft_cigar[i] & 0xf];

      ref_pos = lift; // This is correct since it is the position in the concatenation.
      ref_len = bam_cigar2rlen(bam->core.n_cigar, lft_cigar);
      if (ref_len > 0)
      {
        char* l_ref = (char *)malloc(ref_len);
        ra.expandSubstr(ref_pos, ref_len, l_ref);

        // Convert A,C,G,T,N into 0,1,2,3,4
        for (size_t i = 0; i < ref_len; ++i)
          l_ref[i] = seq_nt4_table[(int)l_ref[i]];
        // Compute the MD:Z field and the number of mismatches
        sam->nm = write_MD_core((uint8_t *)l_ref, seq, lft_cigar, bam->core.n_cigar, nullptr, 0, sam->md);
        sam->rlen = ref_len;

        score.score = ez.score;
        score.pos = ref_pos;

        free(l_ref);
      } else { // Read is unmapped because it align on an insertion of length > readlength
        sam->pos = 0;
        sam->rname = "*";
        sam->cigar = "*";
        sam->rlen = 0;
        sam->unmapped_lft = true;
        score.unmapped_lft = true;
      }
      bam_destroy1(bam);
    



      free(cigar);
      // free(tmp);
    }
    free(ref);
    free(seq);

    if (ez_lc.m_cigar > 0)
      free(ez_lc.cigar);
    if (ez_rc.m_cigar > 0)
      free(ez_rc.cigar);
    if (ez.m_cigar > 0)
      free(ez.cigar);
    for( size_t i = 0; i < ez_cc.size(); ++i )
      if(ez_cc[i].n_cigar > 0)
        free(ez_cc[i].cigar);

    return score;
  }

  // From https://github.com/lh3/ksw2/blob/master/cli.c
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

  std::string to_sam()
  {
      std::string res = "@HD\tVN:1.6\tSO:unknown\n";
      res += idx.to_sam();
      res += "@PG\tID:moni\tPN:moni\tVN:0.1.0\n";
      return res; 
  }

  void to_sam_header(sam_hdr_t *h)
  {
    sam_hdr_add_line(h, "HD", "VN", "1.6", "SO", "unknown", NULL);
    const std::vector<std::string>& names = idx.get_names();
    for (size_t i = 0; i < names.size(); ++i)
        sam_hdr_add_line(h, "SQ", "SN", names[i].c_str(), "LN", std::to_string(idx.length(i)).c_str(), NULL);
    sam_hdr_add_line(h, "PG", "ID", "moni", "VN", "0.1.0", NULL);
  }

  std::string to_csv()
  {
    std::string header = "Read_Name,Num_Unique_MEMs,Total_Num_MEMs,Max_Freq_MEM,Min_Freq_MEM,Highest_Occ_MEM,Lowest_Occ_MEM,Num_Filtered_MEM,Num_Chains_Skipped\n";
    return header;
  }

protected:
    seed_finder_t mem_finder;
    typename seed_finder_t::slp_type& ra;
    liftidx idx;
    // seqidx idx;
    // SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd> ra;

    const size_t min_len = 0;
    const size_t ext_len = 100;   // Extension length
    size_t aligned_reads = 0;
    size_t n = 0;
    size_t top_k = 1; // report the top_k alignments
    size_t check_k = 5; // report the top_k alignments
    size_t region_dist = 10; // report the top_k alignments
    float mapq_coeff_len = 50.0; // report the top_k alignments
    int32_t mapq_coeff_fac = log(mapq_coeff_len); // report the top_k alignments
    
    // Insertion size variables
    double ins_mean = 0.0;        // The mean of the insert size distribution.
    double ins_std_dev = 0.0;     // The standard deviation of the insert size distribution.
    double ins_variance = 0.0;     // The variance of the insert size distribution.
    double ins_sample_variance = 0.0;     // The variance of the insert size distribution.
    double ins_m2 = 0.0;     // The standard deviation of the insert size distribution.
    size_t ins_count = 0;     // The standard deviation of the insert size distribution.
    bool ins_learning_complete = false;
    size_t ins_learning_n = 1000; // Number of unique alignments required to learn the insert size distribution
    size_t ins_learning_score_gap_threshold = 0;
    std::mutex __ins_mtx;


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
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

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

    const int8_t max_penalty = 17; // Largest combined penalty for sub_n computation
    // int8_t max_rseq = 0;

    const int m = 5;
    int8_t mat[25];
    // int minsc = 0, xtra = KSW_XSTART;
    // uint8_t *rseq = 0;

    bool forward_only;

    bool filter_dir = true;
    double dir_thr = 50.0;

    bool filter_seeds = true;
    size_t n_seeds_thr = 5000;

    bool filter_freq = true;  // Filter seed if it occurs with frequency greater than threshold
    double freq_thr = 0.02;   // Filter seed if it occurs with frequency greater than threshold

    ll max_iter = 50;       // Max number of iterations of the chaining algorithhm
    ll max_pred = 50;       // Max number of predecessor to be considered
    ll max_dist_x = 500;    // Max distance for two anchors to be chained
    ll max_dist_y = 100;    // Max distance for two anchors from the same read to be chained
    bool secondary_chains = false; // Find secondary chains in paired-end setting
    bool report_mems = false; // Report MEMs instead of read alignment
    bool csv = false;         // Report MEM statistics in CSV file

    chain_config_t chain_config;
};

#endif /* end of include guard: _ALIGNER_KSW2_HH */
