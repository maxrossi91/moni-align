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

MTIME_TSAFE_INIT(4);
#include <slp_definitions.hpp>
#include <chain.hpp>

#include <htslib/sam.h>

#define _REALIGN


// KSEQ_INIT(gzFile, gzread);

// //*************** Borrowed from minimap2 ***************************************
// static const char LogTable256[256] = {
// #define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
// 	-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
// 	LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
// 	LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
// };

// static inline int ilog2_32(uint32_t v)
// {
// 	uint32_t t, tt;
// 	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
// 	return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
// }
// //******************************************************************************
// ////////////////////////////////////////////////////////////////////////////////
// /// kseq extra
// ////////////////////////////////////////////////////////////////////////////////

// static inline size_t ks_tell(kseq_t *seq)
// {
//   return gztell(seq->f->f) - seq->f->end + seq->f->begin;
// }

// void copy_kstring_t(kstring_t &l, kstring_t &r)
// {
//   l.l = r.l;
//   l.m = r.m;
//   l.s = (char *)malloc(l.m);
//   for (size_t i = 0; i < r.m; ++i)
//     l.s[i] = r.s[i];
// }
// void copy_kseq_t(kseq_t *l, kseq_t *r)
// {
//   copy_kstring_t(l->name, r->name);
//   copy_kstring_t(l->comment, r->comment);
//   copy_kstring_t(l->seq, r->seq);
//   copy_kstring_t(l->qual, r->qual);
//   l->last_char = r->last_char;
// }
// ////////////////////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////////////
// /// helper functions
// ////////////////////////////////////////////////////////////////////////////////

// char complement(char n)
// {
//   switch (n)
//   {
//   case 'A':
//     return 'T';
//   case 'T':
//     return 'A';
//   case 'G':
//     return 'C';
//   case 'C':
//     return 'G';
//   default:
//     return n;
//   }
// }
// ////////////////////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////////////
// /// SLP definitions
// ////////////////////////////////////////////////////////////////////////////////

// using SelSd = SelectSdvec<>;
// using DagcSd = DirectAccessibleGammaCode<SelSd>;
// using Fblc = FixedBitLenCode<>;

// using shaped_slp_t = SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd>;
// using plain_slp_t = PlainSlp<uint32_t, Fblc, Fblc>;

// template< typename slp_t>
// std::string get_slp_file_extension()
// {
//   return std::string(".slp");
// }

// template <>
// std::string get_slp_file_extension<shaped_slp_t>()
// {
//   return std::string(".slp");
// }

// template <>
// std::string get_slp_file_extension<plain_slp_t>()
// {
//   return std::string(".plain.slp");
// }
// ////////////////////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////////////
// /// SAM flags
// ////////////////////////////////////////////////////////////////////////////////

// #define SAM_PAIRED 1                      // template having multiple segments in sequencing
// #define SAM_MAPPED_PAIRED 2               // each segment properly aligned according to the aligner
// #define SAM_UNMAPPED 4                    // segment unmapped
// #define SAM_MATE_UNMAPPED 8               // next segment in the template unmapped
// #define SAM_REVERSED 16                   // SEQ being reverse complemented
// #define SAM_MATE_REVERSED 32              // SEQ of the next segment in the template being reverse complemented
// #define SAM_FIRST_IN_PAIR 64              // the first segment in the template
// #define SAM_SECOND_IN_PAIR 128            // the last segment in the template
// #define SAM_SECONDARY_ALIGNMENT 256       // secondary alignment
// #define SAM_FAILS_CHECKS 512              // not passing filters, such as platform/vendor quality controls
// #define SAM_DUPLICATE 1024                // PCR or optical duplicate
// #define SAM_SUPPLEMENTARY_ALIGNMENT 2048  // supplementary alignment

// ////////////////////////////////////////////////////////////////////////////////

template <typename slp_t,
          typename ms_t>
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
        size_t region_dist = 10; // Maximum distance for two regions to be called from the same region

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

    } config_t;



    typedef struct score_t{
      int32_t score = 0;
      size_t pos = 0; // Position of the leftmost match of the read in the global reference
      size_t lft = 0; // Position of the lift in the reference coordinates
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
      bam1_t bam;

      const size_t min_score;
      score_t score;
      int32_t score2 = 0;

      std::vector<mem_t> mems;
      std::vector<std::pair<size_t, size_t>> anchors;
      std::vector<chain_t> chains;


      alignment_t(kseq_t *read_):
          read(read_),
          sam(read),
          min_score(20 + 8 * log(read->seq.l))
      {
        rc_copy_kseq_t(&read_rev, read);
      }

      void write(samFile* out, const sam_hdr_t *hdr)
      {
        sam_write1(out, hdr, &bam);
      }

      void write(FILE* out)
      {
        write_sam(out, sam);
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
                smatch(config.smatch),          // Match score default
                smismatch(config.smismatch),    // Mismatch score default
                gapo(config.gapo),              // Gap open penalty
                gapo2(config.gapo2),            // Gap open penalty
                gape(config.gape),              // Gap extension penalty
                gape2(config.gape2),            // Gap extension penalty
                end_bonus(config.end_bonus),    // Bonus to add at the extension score to declare the alignment
                w(config.w),                    // Band width
                zdrop(config.zdrop),            // Zdrop enable
                forward_only(config.forward_only)
  {
    verbose("Loading the matching statistics index");
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    std::string filename_ms = filename + ms.get_file_extension();

    ifstream fs_ms(filename_ms);
    ms.load(fs_ms);
    fs_ms.close();

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Matching statistics index loading complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    std::string filename_slp = filename + get_slp_file_extension<slp_t>();
    verbose("Loading random access file: " + filename_slp);
    t_insert_start = std::chrono::high_resolution_clock::now();


    ifstream fs(filename_slp);
    ra.load(fs);
    fs.close();

    n = ra.getLen();

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Random access loading complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    std::string filename_idx = filename + idx.get_file_extension();
    verbose("Loading fasta index file: " + filename_idx);
    t_insert_start = std::chrono::high_resolution_clock::now();


    ifstream fs_idx(filename_idx);
    idx.load(fs_idx);
    fs_idx.close();

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Fasta index loading complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    verbose("Initialize the local aligner");
    t_insert_start = std::chrono::high_resolution_clock::now();

    ksw_gen_simple_mat(m,mat,smatch,-smismatch);

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

  // Fill the vector of occurrences of the mem_t data structure
  // TODO: remove the checks for the lengths once Dominik fix the lce queries
  void find_MEM_occs(mem_t &mem)
  {
    mem.occs.push_back(mem.pos);

    const auto sa_first = ms.get_first_run_sample();
    const auto sa_last = ms.get_last_run_sample();

    // Phi direction
    size_t curr = mem.pos;
    if(curr != sa_first)
    {
      size_t next = ms.Phi(curr);
      if((n-curr) >= mem.len and (n-next) >= mem.len)
      {
        size_t lcp =  lceToRBounded(ra,curr,next,mem.len);
        while(lcp >= mem.len)
        {
          mem.occs.push_back(next);
          
          curr = next;
          if(curr != sa_first)
          {
            next = ms.Phi(curr);
            if((n-curr) >= mem.len and (n-next) >= mem.len)
              lcp  = lceToRBounded(ra,curr,next,mem.len);
            else
              lcp = 0;
          }else lcp = 0;
        }
      }
    }

    // Phi_inv direction
    curr = mem.pos;
    if(curr != sa_last)
    {
      size_t next = ms.Phi_inv(curr);
      if((n-curr) >= mem.len and (n-next) >= mem.len)
      {
        size_t lcp =  lceToRBounded(ra,curr,next,mem.len);
        while(lcp >= mem.len)
        {
          mem.occs.push_back(next);
          
          curr = next;
          if(curr != sa_last)
          {
            next = ms.Phi_inv(curr);
            if((n-curr) >= mem.len and (n-next) >= mem.len)
              lcp  = lceToRBounded(ra,curr,next,mem.len);
            else
              lcp = 0;
          }else lcp = 0;
        }
      }
    }

  }



  bool align(kseq_t *read, samFile* out, const sam_hdr_t *hdr)
  {
    // std::vector<mem_t> mems;
    alignment_t alignment(read);
    if(not align(alignment))
      alignment.set_sam_not_aligned();
    alignment.write(out, hdr);
    return alignment.aligned;
  }

  bool align(kseq_t *read, FILE* out)
  {
    // std::vector<mem_t> mems;
    alignment_t alignment(read);
    if(not align(alignment))
      alignment.set_sam_not_aligned();
    alignment.write(out);
    return alignment.aligned;
  }

  // Aligning unpaired sequences
  bool align(alignment_t &al)
  {    
    MTIME_INIT(3);
    MTIME_START(0); // Timing helper

    find_seeds(al.read,al.mems, 0, MATE_1 | MATE_F);
    find_seeds(&al.read_rev,al.mems, 0, MATE_1 | MATE_RC);

    MTIME_END(0); //Timing helper
    MTIME_START(1); //Timing helper

    // Chain MEMs
    al.chained = find_chains(al.mems, al.anchors, al.chains);
 
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

    while (i < al.chains.size() and different_scores.size() < check_k)
    {
        different_scores.insert(al.chains[i].score);
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
              else if (score.score == std::get<0>(best_scores[j]) )
                j = best_scores.size(), replaced = true, i++;

          if( not replaced )
            best_scores.push_back(std::make_tuple(score.score, score.lft, i++));
        }
    }
    
    if (best_scores.size() < 2)
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


  typedef struct orphan_paired_score_t{
    int32_t tot = 0;
    score_t m1;
    score_t m2; 
    std::pair<size_t,size_t> pos = std::make_pair(0,0); // Position of the mate
  } orphan_paired_score_t;

  typedef struct paired_score_t{
    int32_t tot = 0;
    score_t m1;
    score_t m2; 
    bool paired = false;


    paired_score_t& operator=(orphan_paired_score_t const& copy)
    {
      tot = copy.tot;
      m1 = copy.m1;
      m2 = copy.m2;
      return *this;
    }
  } paired_score_t;

  typedef struct paired_alignment_t{
    bool aligned = false;
    bool chained = false;
    bool best_score = false;

    kseq_t* mate1 = nullptr;
    kseq_t* mate2 = nullptr;

    kseq_t mate1_rev;
    kseq_t mate2_rev;

    sam_t sam_m1;
    sam_t sam_m2;

    // TODO: precompute the nt4 version of the mates

    size_t min_score = 0;
    paired_score_t score;
    int32_t score2 = 0;

    std::vector<mem_t> mems;
    std::vector<std::pair<size_t, size_t>> anchors;
    std::vector<chain_t> chains;

    paired_alignment_t(kseq_t *mate1_, kseq_t *mate2_):
      mate1(mate1_),
      mate2(mate2_),
      sam_m1(mate1),
      sam_m2(mate2),
      min_score(20 + 8 * log(mate1->seq.l))
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
    
    void init(kseq_t *mate1_, kseq_t *mate2_)
    {
      mate1 = mate1_;
      mate2 = mate2_;
      sam_m1 = sam_t(mate1);
      sam_m2 = sam_t(mate2);
      min_score = 20 + 8 * log(mate1->seq.l);
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
      release_memory();
    }

  } paired_alignment_t;

  // Aligning pair-ended batched sequences
  size_t align(kpbseq_t* batch, FILE *out)
  {
    size_t n_aligned = 0;

    std::vector<paired_alignment_t> alignments(batch->mate1->l);
    std::vector<double> mate_abs_distance;

    // Computing Mean and Variance using Welford's algorithm
    size_t count = 0;   // Number of samples
    double mean = 0.0;  // Accumulates the mean
    double m2 = 0.0;    // Accumulates the squared distance from the mean


    int l = batch->mate1->l;
    for (size_t i = 0; i < l; ++i)
    {
      // paired_alignment_t alignment(&batch->mate1->buf[i], &batch->mate2->buf[i]);
      paired_alignment_t& alignment = alignments[i];
      alignment.init(&batch->mate1->buf[i], &batch->mate2->buf[i]);
      if(align(alignment))
      {
        // Get stats
        // mate_abs_distance.push_back((double)(alignment.sam_m1.tlen >= 0?alignment.sam_m1.tlen :-alignment.sam_m1.tlen));
        double value = (double)(alignment.sam_m1.tlen >= 0?alignment.sam_m1.tlen :-alignment.sam_m1.tlen);
        double delta = value - mean;
        mean += delta / (++count);
        m2 += delta * (value - mean);
      }
      // Store the alignment informations
      // alignments.push_back(alignment);
    }

    // Computes stats
    double variance = m2/count;
    double sampleVariance = m2/(count-1);
    // double variance = 0.0;
    // double mean = mate_abs_distance[0];
    // for(size_t i = 0 ; i < l; ++i )
    // {
    //   mean += mate_abs_distance[i];
    //   double diff = ((i+1) * mate_abs_distance[i] ) - mean;
    //   variance += (diff * diff) / ((i + 1.0 ) * i);
    // }

    // variance /= (double)(l-1);
    // mean /= (double)l;
    double std_dev = sqrt(variance);

    // Perform local serach for unaligned mates.
    for(auto& alignment : alignments)
    {
      if(not alignment.aligned and alignment.chained)
        orphan_recovery(alignment, mean, std_dev);
      // Write alignment to file
      alignment.write(out);
      if(alignment.aligned) ++n_aligned;
    }
    
    return n_aligned;
  }

  // Aligning pair-ended sequences
  bool align(kseq_t *mate1, kseq_t *mate2, FILE *out)
  {
    paired_alignment_t alignment(mate1, mate2);
    if(not align(alignment))
      alignment.set_sam_not_aligned();
    alignment.write(out);
    return alignment.aligned;
  }

  // Aligning pair-ended sequences
  bool align(paired_alignment_t &al)
  { 
    MTIME_INIT(3);   
    MTIME_START(0); // Timing helper

    // Find MEMs
    find_seeds(al.mate1, al.mems, 0, MATE_1 | MATE_F);
    find_seeds(&al.mate1_rev, al.mems, al.mate2->seq.l, MATE_1 | MATE_RC );
    find_seeds(al.mate2, al.mems, 0, MATE_2 | MATE_F);
    find_seeds(&al.mate2_rev, al.mems, al.mate1->seq.l, MATE_2 | MATE_RC);
    // find_mems(al.mate1, al.mems, 0, MATE_1 | MATE_F);
    // find_mems(&al.mate1_rev, al.mems, al.mate2->seq.l, MATE_1 | MATE_RC );
    // find_mems(al.mate2, al.mems, 0, MATE_2 | MATE_F);
    // find_mems(&al.mate2_rev, al.mems, al.mate1->seq.l, MATE_2 | MATE_RC);

    MTIME_END(0); //Timing helper
    MTIME_START(1); //Timing helper

    // Chain MEMs
    al.chained = find_chains(al.mems, al.anchors, al.chains);

    MTIME_END(1);   //Timing helper
    MTIME_START(2); //Timing helper

    if (not al.chained)
    {
      MTIME_END(2); //Timing helper
      MTIME_TSAFE_MERGE;
      return false;
    }


    int32_t min_score = 20 + 8 * log(al.mate1->seq.l);

    // // Compute the second best score
    std::vector<std::tuple<int32_t, size_t, size_t, size_t>> best_scores;
    // get the occurrences of the top 4 best scores
    std::set<size_t> different_scores;
    size_t i = 0;
    while (i < al.chains.size() and different_scores.size() < check_k)
    {
        different_scores.insert(al.chains[i].score);
        if (different_scores.size() < check_k)
        {
          // Align the chain
          auto chain = al.chains[i];
          // Reverse the chain order
          std::reverse(chain.anchors.begin(), chain.anchors.end());
          // Compute the score of a chain.
          paired_score_t score;
          if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
            score = paired_chain_score(chain, al.anchors, al.mems, min_score, al.mate1, &al.mate2_rev);
          else
            score = paired_chain_score(chain, al.anchors, al.mems, min_score, &al.mate1_rev, al.mate2);
          

          score.m1.lft = idx.lift(score.m1.pos);
          score.m2.lft = idx.lift(score.m2.pos);

          // Check if there is another read scored in the same region
          bool replaced = false;
          for(size_t j = 0; j < best_scores.size(); ++j)
          {
            auto d1 = dist(std::get<1>(best_scores[j]),score.m1.lft);
            auto d2 = dist(std::get<2>(best_scores[j]),score.m2.lft);
            if( (dist(std::get<1>(best_scores[j]),score.m1.lft) < region_dist) and 
                (dist(std::get<2>(best_scores[j]),score.m2.lft) < region_dist) )
              if ( score.tot > std::get<0>(best_scores[j]) )
                if ( replaced )
                  best_scores[j] = std::make_tuple(0,0,0,i-1);
                else
                  best_scores[j] = std::make_tuple(score.tot, score.m1.lft, score.m2.lft, i++), replaced = true;
              else if ( score.tot == std::get<0>(best_scores[j]) )
                j = best_scores.size(), replaced = true, i++;
          }
          if( not replaced )
            best_scores.push_back(std::make_tuple(score.tot, score.m1.lft, score.m2.lft, i++));
        }
    }
    
    if (best_scores.size() < 2)
      best_scores.push_back(std::make_tuple(0, 0, 0, al.chains.size()));

    std::sort(best_scores.begin(), best_scores.end(), std::greater<std::tuple<int32_t, size_t, size_t, size_t>>());

    assert(best_scores.size() > 1);

    if (std::get<0>(best_scores[0]) < min_score)
    {
      MTIME_END(2); //Timing helper
      MTIME_TSAFE_MERGE;
      return false;
    }

    al.best_score = true;

    al.score2 = std::get<0>(best_scores[1]);


    { // Forward case
      i = std::get<3>(best_scores[0]);
      // Align the chain
      auto chain = al.chains[i];
      // Reverse the chain order
      std::reverse(chain.anchors.begin(), chain.anchors.end());
      // Compute the score of a chain.
      if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
        al.score = paired_chain_score(chain, al.anchors, al.mems, min_score, al.mate1, &al.mate2_rev, false, al.score2, 0, &al.sam_m1, &al.sam_m2);
      else
        al.score = paired_chain_score(chain, al.anchors, al.mems, min_score, &al.mate1_rev, al.mate2, false, al.score2, 1, &al.sam_m1, &al.sam_m2);
    }

    al.aligned = (al.score.tot >= min_score);

    MTIME_END(2); //Timing helper
    MTIME_TSAFE_MERGE;
    if (not al.aligned)
      return false;
    

    return al.aligned;
  }

  bool orphan_recovery(
    paired_alignment_t& alignment, 
    const double mean, 
    const double std_dev)
  {
    MTIME_INIT(3);   
    MTIME_START(2); //Timing helper
    // We need to use the local alignment of ksw to perform local search
    
    // For all good chaining of both mates of length at least min_length, find the possible distance of the mate
    // and 

    // Compute the second best score
    std::vector<std::pair<std::pair<int32_t, std::pair<size_t,size_t> >, size_t>> best_scores;
    size_t i = 0;
    while (i < alignment.chains.size())
    {
      // Align the chain
      auto chain = alignment.chains[i];
      // Reverse the chain order
      std::reverse(chain.anchors.begin(), chain.anchors.end());
      // Compute the score of a chain.
      orphan_paired_score_t score;
      if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
        score = paired_chain_orphan_score(chain, alignment.anchors, alignment.mems, alignment.min_score, alignment.mate1, &alignment.mate2_rev, mean, std_dev);
      else
        score = paired_chain_orphan_score(chain, alignment.anchors, alignment.mems, alignment.min_score, &alignment.mate1_rev, alignment.mate2, mean, std_dev);
      
      best_scores.push_back(std::make_pair(std::make_pair(score.tot,score.pos), i++));
    }

    orphan_paired_score_t zero;
    if (best_scores.size() < 2)
      best_scores.push_back(std::make_pair(std::make_pair(0,std::make_pair(0,0)), alignment.chains.size()));

    std::sort(best_scores.begin(), best_scores.end(), std::greater<std::pair<std::pair<int32_t, std::pair<size_t,size_t> >, size_t>>());

    assert(best_scores.size() > 1);

    if (best_scores[0].first.first < alignment.min_score)
    {
      MTIME_END(2); //Timing helper
      MTIME_TSAFE_MERGE;
      return false;
    }
    alignment.best_score = true;

    alignment.score2 = best_scores[1].first.first;


    { // Forward case
      i = best_scores[0].second;
      ll start = best_scores[0].first.second.first;
      ll end = best_scores[0].first.second.second;
      // Align the chain
      auto chain = alignment.chains[i];
      // Reverse the chain order
      std::reverse(chain.anchors.begin(), chain.anchors.end());
      // Compute the score of a chain.
      if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
        alignment.score = paired_chain_orphan_score(chain, alignment.anchors, alignment.mems, alignment.min_score, alignment.mate1, &alignment.mate2_rev, mean, std_dev, false, alignment.score2, 0, start, end,  &alignment.sam_m1, &alignment.sam_m2);
      else
        alignment.score = paired_chain_orphan_score(chain, alignment.anchors, alignment.mems, alignment.min_score, &alignment.mate1_rev, alignment.mate2, mean, std_dev, false, alignment.score2, 1, start, end, &alignment.sam_m1, &alignment.sam_m2);
    }

    alignment.aligned = (alignment.score.tot >= alignment.min_score);

    MTIME_END(2); //Timing helper
    MTIME_TSAFE_MERGE;
    if (not alignment.aligned)
      return false;
    

    return alignment.aligned;
  }

  void find_mems(
    const kseq_t *read,
    std::vector<mem_t>& mems,
    size_t r_offset = 0,
    size_t mate = 0
    ) 
  {
    auto pointers = ms.query(read->seq.s, read->seq.l);
    size_t l = 0;   // Current match length
    size_t pl = 0;  // Previous match length
    size_t n_Ns = 0;
    for (size_t i = 0; i < pointers.size(); ++i)
    {
      size_t pos = pointers[i];
      while ((i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
      {
        if(read->seq.s[i + l] == 'N') n_Ns++;
        else n_Ns = 0;
        ++l;
      }

      // Update MEMs
      if (l >= pl and n_Ns < l and l >= min_len)
      {
        size_t r = r_offset + (i + l - 1); // compatible with minimap2 chaining algorithm
        // size_t r = r_offset + ((reverse) ? (i) : (i + l - 1)); // compatible with minimap2 chaining algorithm
        // size_t r = i + l - 1;
        // r = (reverse)? (read->seq.l - (r + 1 - l) - 1) : (r); // compatible with minimap2 chaining algorithm
        mems.push_back(mem_t(pointers[i],l,i,mate,r));
        find_MEM_occs(mems.back());
      }

      // Compute next match length
      pl = l;
      l = (l == 0 ? 0 : (l - 1));
    }

  }

  void find_seeds(
    const kseq_t *read,
    std::vector<mem_t>& mems,
    size_t r_offset = 0,
    size_t mate = 0
    ) 
  {
    auto pointers = ms.query(read->seq.s, read->seq.l);
    size_t l = 0;   // Current match length
    size_t pl = 0;  // Previous match length
    size_t n_Ns = 0;
    for (size_t i = 0; i < pointers.size(); ++i)
    {
      size_t pos = pointers[i];
      while ((i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
      {
        if(read->seq.s[i + l] == 'N') n_Ns++;
        else n_Ns = 0;
        ++l;
      }

      // Update MEMs
      if (l >= pl and n_Ns < l and l >= min_len)
      {
        size_t r = r_offset + (i + l - 1); // compatible with minimap2 chaining algorithm
        mems.push_back(mem_t(pointers[i],l,i,mate,r));
        find_MEM_occs(mems.back());
        // Take two halves of the MEM
        if(l >= (min_len << 1))
        {
          size_t ll = l >> 1; 
          size_t rl = r_offset + (i + ll - 1); // compatible with minimap2 chaining algorithm
          mems.push_back(mem_t(pointers[i],ll,i,mate,rl));
          find_MEM_occs(mems.back()); // TODO: Optimize this

          size_t lr = l - ll;
          size_t rr = r_offset + (i + l - 1); // compatible with minimap2 chaining algorithm
          mems.push_back(mem_t(pointers[i + ll],lr,i + ll,mate,rr));
          find_MEM_occs(mems.back()); // TODO: Optimize this
        }
      }

      // Compute next match length
      pl = l;
      l = (l == 0 ? 0 : (l - 1));
    }

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
      sam_t* sam = nullptr)     // The SAM information pointer)
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

        fill_chain(mems,chain_anchors,lcs,lcs_len,rcs,rcs_len,read,score_only,sam); 
        
        sam->flag = (strand?16:0);
        sam->zs = score2;
        sam->mapq = compute_mapq(sam->as, sam->zs, min_score, sam->read->seq.l * smatch);
        // sam->reverse = (strand != 0);
        // sam->rname = idx[sam->pos - 1];
        if(output)
        {
          write_sam(out,*sam);
          delete sam;
        }
        
        // fill_chain(mems,chain_anchors,lcs,lcs_len,rcs,rcs_len,read,score_only,0,min_score,strand,out); 
      }

      delete rcs;
      delete lcs;

      return score;
  }

  paired_score_t paired_chain_score(
      const chain_t &chain,
      const std::vector<std::pair<size_t, size_t>> &anchors,
      const std::vector<mem_t> &mems,
      const int32_t min_score,
      const kseq_t *mate1,
      const kseq_t *mate2,
      const bool score_only = true,
      const int32_t score2 = 0,
      const uint8_t strand = 0,
      sam_t *sam_m1_ = nullptr,
      sam_t *sam_m2_ = nullptr)
  {
    paired_score_t score;
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
          score.m1 = chain_score(mate1_chain, anchors, mems, min_score, mate1);
          score.m2 = chain_score(mate2_chain, anchors, mems, min_score, mate2);
          score.tot = score.m1.score + score.m2.score;
      }
      else
      {
        sam_t& sam_m1 = *sam_m1_;
        sam_t& sam_m2 = *sam_m2_;


        score.m1 = chain_score(mate1_chain, anchors, mems, min_score, mate1, false, score2, strand, nullptr, &sam_m1);
        score.m2 = chain_score(mate2_chain, anchors, mems, min_score, mate2, false, score2, strand, nullptr, &sam_m2);
        score.tot = score.m1.score + score.m2.score;

        sam_m1.read = mate1;
        sam_m2.read = mate2;
        // Fill sam fields RNEXT, PNEXT and TLEN
        // sam_m1.rnext = std::string(mate2->name.s);
        // sam_m2.rnext = std::string(mate1->name.s);

        if(score.m1.score >= min_score and score.m2.score >= min_score)
        {

          sam_m1.pnext = sam_m2.pos;
          sam_m2.pnext = sam_m1.pos;

          ll tlen = (sam_m2.pos + mate2->seq.l) - sam_m1.pos;

          sam_m1.tlen = tlen;
          sam_m2.tlen = -tlen;

          sam_m1.rname = idx[sam_m1.pos - 1]; // Check if necessary
          sam_m2.rname = idx[sam_m2.pos - 1];

          size_t mapq = compute_mapq(score.tot, score2, min_score, (sam_m1.read->seq.l + sam_m2.read->seq.l) * smatch);
          
          sam_m1.mapq = mapq;
          sam_m2.mapq = mapq;

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
        }else if(score.m1.score >= min_score) {
          sam_m1.rname = idx[sam_m1.pos - 1];

          sam_m1.flag = SAM_PAIRED | SAM_MATE_UNMAPPED | SAM_FIRST_IN_PAIR;
          sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_SECOND_IN_PAIR;
          if(strand)
            sam_m1.flag |= SAM_REVERSED;
        }else if(score.m2.score >= min_score) {
          sam_m2.rname = idx[sam_m2.pos - 1];

          sam_m1.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_FIRST_IN_PAIR;
          sam_m2.flag = SAM_PAIRED | SAM_MATE_UNMAPPED | SAM_SECOND_IN_PAIR;
          if(not strand)
            sam_m2.flag |= SAM_REVERSED;
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



orphan_paired_score_t paired_chain_orphan_score(
      const chain_t &chain,
      const std::vector<std::pair<size_t, size_t>> &anchors,
      const std::vector<mem_t> &mems,
      const int32_t min_score,
      const kseq_t *mate1,
      const kseq_t *mate2,
      const double mean,
      const double std_dev,
      const bool score_only = true,
      const int32_t score2 = 0,
      const uint8_t strand = 0,
      long long int start = 0,  // Start position of the reference
      long long int end = 0,    //End position of the reference
      sam_t *sam_m1_ = nullptr,
      sam_t *sam_m2_ = nullptr)
  {
    orphan_paired_score_t score;

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
        score.m1 = chain_score(mate1_chain, anchors, mems, min_score, mate1);
        // Perform local search
        // TODO: Check whether we do not cross one boundary
        start = rm_pos + (ll)std::floor(mean - 4*std_dev);
        end = rm_pos + (ll)std::ceil(mean + 4*std_dev);
        maxl(start,(ll)0);
        minl(end,(ll)n);

        score.m2 = fill_orphan(start,end,mate2);  
        score.pos = std::pair(start,end);

      }else
      {
        score.m2 = chain_score(mate2_chain, anchors, mems, min_score, mate2);
        // Perform local search
        // TODO: Check whether we do not cross one boundary
        start = lm_pos + (ll)std::floor(-mean - 4*std_dev);
        end = lm_pos + (ll)std::ceil(-mean + 4*std_dev);
        maxl(start,(ll)0);
        minl(end,(ll)n);

        score.m1 = fill_orphan(start,end,mate1);  
        score.pos = std::pair(start,end);
      }
      score.tot = score.m1.score + score.m2.score;
    }
    else
    {
      sam_t& sam_m1 = *sam_m1_;
      sam_t& sam_m2 = *sam_m2_;

      if(mate1_chain.second.size() > 0)
      {
        score.m1 = chain_score(mate1_chain, anchors, mems, min_score, mate1, false, score2, strand, nullptr, &sam_m1);
        score.m2 = fill_orphan(start,end,mate2,false,&sam_m2);
      }
      else
      {
        score.m1 = fill_orphan(start,end,mate1,false,&sam_m1);
        score.m2 = chain_score(mate2_chain, anchors, mems, min_score, mate2, false, score2, strand, nullptr, &sam_m2);
      }
      score.tot = score.m1.score + score.m2.score;

      sam_m1.read = mate1;
      sam_m2.read = mate2;
      // Fill sam fields RNEXT, PNEXT and TLEN
      // sam_m1.rnext = std::string(mate2->name.s);
      // sam_m2.rnext = std::string(mate1->name.s);

      if(score.m1.score >= min_score and score.m2.score >= min_score)
      {

        sam_m1.pnext = sam_m2.pos;
        sam_m2.pnext = sam_m1.pos;

        ll tlen = (sam_m2.pos + mate2->seq.l) - sam_m1.pos;

        sam_m1.tlen = tlen;
        sam_m2.tlen = -tlen;

        // sam_m1.rname = idx[sam_m1.pos - 1];
        // sam_m2.rname = idx[sam_m2.pos - 1];
        size_t mapq = compute_mapq(score.tot, score2, min_score, (sam_m1.read->seq.l + sam_m2.read->seq.l) * smatch);
        
        sam_m1.mapq = mapq;
        sam_m2.mapq = mapq;


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
      }else if(score.m1.score >= min_score) {
        // sam_m1.rname = idx[sam_m1.pos - 1];

        sam_m1.flag = SAM_PAIRED | SAM_MATE_UNMAPPED | SAM_FIRST_IN_PAIR;
        sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_SECOND_IN_PAIR;
        if(strand)
          sam_m1.flag |= SAM_REVERSED;
      }else if(score.m2.score >= min_score) {
        // sam_m2.rname = idx[sam_m2.pos - 1];

        sam_m1.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_FIRST_IN_PAIR;
        sam_m2.flag = SAM_PAIRED | SAM_MATE_UNMAPPED | SAM_SECOND_IN_PAIR;
        if(not strand)
          sam_m2.flag |= SAM_REVERSED;
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
    {
      kswq_t *q = 0;
      kswr_t r;
      const int xtra = KSW_XSTART;

      r = ksw_align(seq_len, (uint8_t *)seq, ref_len, (uint8_t *)ref, 5, mat, gapo, gape, xtra, &q);

      // Update start and end
      end = start + r.te;
      start += r.tb;
      score.score= r.score;
      score.pos = start;
    }

    if(not score_only)
    {
        // Realign the whole sequence globally
        int flag = KSW_EZ_RIGHT;

        ksw_extz_t ez;
        memset(&ez, 0, sizeof(ksw_extz_t));

        char* tmp = (char*)calloc(max(ref_len,seq_len),1);

        ksw_reset_extz(&ez);
        ksw_extd2_sse(km, seq_len, (uint8_t *)seq, ref_len, (uint8_t*)ref, m, mat, gapo, gape, gapo2, gape2, w, zdrop, end_bonus, flag, &ez);

        // Concatenate the CIGAR strings

        // std::string cigar_s;
        sam->lift_cigar = "";
        for (size_t i = 0; i < ez.n_cigar; ++i)
          sam->lift_cigar += std::to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];

        // sam->n_cigar = ez.n_cigar;
        // sam->cigar_b = (uint32_t*) malloc(ez.n_cigar * sizeof(uint32_t));
        // std::memcopy(sam->cigar_b, ez.cigar, ez.n_cigar * sizeof(uint32_t));

        // Compute the MD:Z field and the number of mismatches
        sam->lift_nm = write_MD_core((uint8_t *)ref, seq, ez.cigar, ez.n_cigar, tmp, 0, sam->lift_md);

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
        char * l_ref = (char *)malloc( ref_len );
        ra.expandSubstr(ref_occ, ref_len, l_ref);

        // Convert A,C,G,T,N into 0,1,2,3,4
        for (size_t i = 0; i < ref_len; ++i)
          l_ref[i] = seq_nt4_table[(int)l_ref[i]];
        // Compute the MD:Z field and the number of mismatches
        sam->nm = write_MD_core((uint8_t *)l_ref, seq, cigar, bam->core.n_cigar, tmp, 0, sam->md);
        sam->rlen = ref_len;

        score.score = ez.score;
        score.pos = ref_occ;
        delete tmp;
        delete l_ref;
        bam_destroy1(bam);
    }

    delete ref;
    delete seq;
    return score;  
  }

  size_t get_aligned_reads()
  {
    return aligned_reads;
  }


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
      ra.expandSubstr(lc_occ, lc_len, tmp_lc);
      // verbose("lc: " + std::string(lc));
      // Convert A,C,G,T,N into 0,1,2,3,4
      // The left context is reversed
      uint8_t *lc = (uint8_t *)malloc(ext_len);
      for (size_t i = 0; i < lc_len; ++i)
        lc[lc_len -i -1] = seq_nt4_table[(int)tmp_lc[i]];
      
      delete tmp_lc;

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

      delete lc;
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
      delete rc;
    }


    // Compute the partial score score
    score.score = score_lc + score_rc;
    // int32_t score = mem_len * smatch + score_lc + score_rc;

    // Compute starting position in reference
    size_t mem_pos = mems[anchors[0].first].occs[anchors[0].second];
    size_t mem_len = mems[anchors.back().first].occs[anchors.back().second] + mems[anchors.back().first].len - mem_pos; // from the strart of the first MEM to the end of the last MEM.

    size_t ref_pos = mem_pos - (lcs_len > 0 ? ez_lc.mqe_t + 1 : 0);
    size_t ref_len = (lcs_len > 0 ? ez_lc.mqe_t + 1 : 0) + mem_len + (rcs_len > 0 ? ez_rc.mqe_t + 1: 0);
    char *ref = (char *)malloc(ref_len);
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

    if(not score_only)
    {
      // Compute starting position in reference

      char* tmp = (char*)calloc(max(ref_len,seq_len),1);

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
      sam->lift_nm = write_MD_core((uint8_t*)ref,seq,cigar,n_cigar,tmp,0,sam->lift_md);

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
      char* l_ref = (char *)malloc(ref_len);
      ra.expandSubstr(ref_pos, ref_len, l_ref);

      // Convert A,C,G,T,N into 0,1,2,3,4
      for (size_t i = 0; i < ref_len; ++i)
        l_ref[i] = seq_nt4_table[(int)l_ref[i]];
      // Compute the MD:Z field and the number of mismatches
      sam->nm = write_MD_core((uint8_t *)l_ref, seq, lft_cigar, bam->core.n_cigar, tmp, 0, sam->md);
      sam->rlen = ref_len;

      delete l_ref;
      bam_destroy1(bam);
    



      delete cigar;
      delete tmp;
    }
    delete ref;
    delete seq;

    if (ez_lc.m_cigar > 0)
      delete ez_lc.cigar;
    if (ez_rc.m_cigar > 0)
      delete ez_rc.cigar;
    if (ez.m_cigar > 0)
      delete ez.cigar;
    for( size_t i = 0; i < ez_cc.size(); ++i )
      if(ez_cc[i].n_cigar > 0)
        delete ez_cc[i].cigar;

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

protected:
    ms_t ms;
    slp_t ra;
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

    // int8_t max_rseq = 0;

    const int m = 5;
    int8_t mat[25];
    // int minsc = 0, xtra = KSW_XSTART;
    // uint8_t *rseq = 0;

    bool forward_only;
};

#endif /* end of include guard: _ALIGNER_KSW2_HH */
