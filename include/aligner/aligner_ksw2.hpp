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

#include <ksw2.h>

#include <libgen.h>
#include <kpbseq.h>

#include <seqidx.hpp>
#include <slp_definitions.hpp>
#include <chain.hpp>

#define _REALIGN

MTIME_INIT(3);

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


    aligner(std::string filename, 
            config_t config = config_t()) : 
                min_len(config.min_len),        // Minimum MEM length
                ext_len(config.ext_len),        // Extension length
                top_k(config.top_k),            // Report the top_k alignments
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


  // bool align_chains(kseq_t *read, FILE* out, uint8_t strand)
  bool align(kseq_t *read, FILE* out)
  {
    // std::vector<mem_t> mems;


    bool aligned = false;

    // Find MEMs in the forward direction
    std::vector<mem_t> mems;

    find_mems(read,mems);

    // for(size_t i = 0; i < mems.size(); ++i){
    //   std::cout << "MEM[" << i <<"]: \n";
    //   std::cout << "    len:" << mems[i].len <<"\n";
    //   std::cout << "    pos:" << mems[i].pos <<"\n";
    //   std::cout << "    idx:" << mems[i].idx <<"\n";
    //   std::cout << "   occs:" << mems[i].occs.size() <<"\n";
    // }

    std::vector< std::pair< size_t, size_t > > anchors;
    std::vector<chain_t> chains;

    // if(!find_chains(mems,anchors,chains))
    //   return false;
    const bool fwd_chains = find_chains(mems, anchors, chains);

    // Find MEMs in reverse direction
    std::vector<mem_t> mems_rev;

    //copy seq
    kseq_t read_rev;
    copy_kseq_t(&read_rev, read);

    for (size_t i = 0; i < read->seq.l; ++i)
      read_rev.seq.s[i] = complement(read->seq.s[read->seq.l - i - 1]);

    if (read_rev.seq.m > read_rev.seq.l)
      read_rev.seq.s[read_rev.seq.l] = 0;


    find_mems(&read_rev,mems_rev,1);

    // for (size_t i = 0; i < mems_rev.size(); ++i)
    // {
    //   std::cout << "MEM[" << i << "]: \n";
    //   std::cout << "    len:" << mems_rev[i].len << "\n";
    //   std::cout << "    pos:" << mems_rev[i].pos << "\n";
    //   std::cout << "    idx:" << mems_rev[i].idx << "\n";
    //   std::cout << "   occs:" << mems_rev[i].occs.size() << "\n";
    // }

    std::vector< std::pair< size_t, size_t > > anchors_rev;
    std::vector<chain_t> chains_rev;

    const bool rev_chains = find_chains(mems_rev, anchors_rev, chains_rev);
    
    if ((not fwd_chains) and (not rev_chains))
    {
      std::string dummy = "";
      write_sam(0, 0, 0, 0, "*", read, 0, out, dummy, dummy, 0, "*", 0, 0);
      return false;
    } 

    int32_t min_score = 20 + 8 * log(read->seq.l);

    // Compute the second best score
    std::vector<std::pair<int32_t, size_t>> best_scores;
    // get the occurrences of the top 4 best scores
    std::set<size_t> different_scores;
    size_t i = 0;
    size_t j = 0;
    size_t off = chains.size();
    while (i < chains.size() and j < chains_rev.size() and different_scores.size() < 3)
    {
      if(chains[i].score > chains_rev[j].score)
      {
        different_scores.insert(chains[i].score);
        if (different_scores.size() < 3)
        {
          // Align the chain
          auto chain = std::make_pair(chains[i].score, chains[i].anchors);
          // Reverse the chain order
          std::reverse(chain.second.begin(), chain.second.end());
          // Compute the score of a chain.
          int32_t score = chain_score(chain, anchors, mems, min_score, read);
          best_scores.push_back(std::make_pair(score, i++));
        }
      }
      else
      {
        different_scores.insert(chains_rev[j].score);
        if (different_scores.size() < 3)
        {
          // Align the chain
          auto chain = std::make_pair(chains_rev[j].score, chains_rev[j].anchors);
          // Reverse the chain order
          std::reverse(chain.second.begin(), chain.second.end());
          // Compute the score of a chain.
          int32_t score = chain_score(chain, anchors_rev, mems_rev, min_score, &read_rev);
          best_scores.push_back(std::make_pair(score, off + j++));
        }
      }
    }
    while (different_scores.size() < 3 and i < chains.size())
    {
      different_scores.insert(chains[i].score);
      if (different_scores.size() < 3)
      {
        // Align the chain
        auto chain = std::make_pair(chains[i].score, chains[i].anchors);
        // Reverse the chain order
        std::reverse(chain.second.begin(), chain.second.end());
        // Compute the score of a chain.
        int32_t score = chain_score(chain, anchors, mems, min_score, read);
        best_scores.push_back(std::make_pair(score, i++));
      }
    }
    while (different_scores.size() < 3 and j < chains_rev.size())
    {
      different_scores.insert(chains_rev[j].score);
      if(different_scores.size() < 3)
      {
        // Align the chain
        auto chain = std::make_pair(chains_rev[j].score, chains_rev[j].anchors);
        // Reverse the chain order
        std::reverse(chain.second.begin(), chain.second.end());
        // Compute the score of a chain.
        int32_t score = chain_score(chain, anchors_rev, mems_rev, min_score, &read_rev);
        best_scores.push_back(std::make_pair(score, off + j++));
      }
    }

    if(best_scores.size() < 2)
      best_scores.push_back(std::make_pair(0,chains.size() + chains_rev.size()));

    std::sort(best_scores.begin(), best_scores.end(), std::greater<std::pair<int32_t, size_t>>());

    assert(best_scores.size() > 1);

    if(best_scores[0].first < min_score)
    {
      std::string dummy = "";
      write_sam(0, 0, 0, 0, "*", read, 0, out, dummy, dummy, 0, "*", 0, 0);
      return false;
    }

    int32_t score2 = best_scores[1].first;
    int32_t score = 0;

    if(best_scores[0].second >= off)
    { // Reverse case
      j = best_scores[0].second - off;
      // Align the chain
      auto chain = std::make_pair(chains_rev[j].score, chains_rev[j].anchors);
      // Reverse the chain order
      std::reverse(chain.second.begin(), chain.second.end());
      // Compute the score of a chain.
      score = chain_score(chain, anchors_rev, mems_rev, min_score, &read_rev, false, score2, 1, out);
    }
    else
    { // Forward case
      i = best_scores[0].second;
      // Align the chain
      auto chain = std::make_pair(chains[i].score, chains[i].anchors);
      // Reverse the chain order
      std::reverse(chain.second.begin(), chain.second.end());
      // Compute the score of a chain.
      score = chain_score(chain, anchors, mems, min_score, read, false, score2, 0, out);
    }

    if (score >= min_score) aligned = true;
    // TODO: Implement the topk retrival

    // std::vector<std::pair<size_t, size_t>> top4;
    // for(size_t i = 0; i < min(chains.size(),(size_t)2); ++i)
    //   top4.push_back(std::make_pair(chains[i].first,i));
    // for(size_t i = 0; i < min(chains_rev.size(),(size_t)2); ++i)
    //   top4.push_back(std::make_pair(chains_rev[i].first,2 + i));

    // std::sort(top4.begin(), top4.end(),std::greater<std::pair<size_t,size_t>>());

    // int32_t score2 = 0;
    // int32_t score = 0;

    // //TODO: Rewrite this part that is so ugly
    // if(top4.size() > 1)
    // {
    //   if(top4[1].second > 1)
    //   {
    //     auto chain = chains_rev[top4[1].second-2];
    //     // Reverse the chain order
    //     std::reverse(chain.second.begin(), chain.second.end());
    //     // Compute the score of a chain.
    //     score2 = chain_score(chain,anchors_rev, mems_rev, min_score,(not (top_k > 1)),-1,&read_rev,16,out);
    //   }
    //   else
    //   {
    //     auto chain = chains[top4[1].second];
    //     // Reverse the chain order
    //     std::reverse(chain.second.begin(), chain.second.end());
    //     // Compute the score of a chain.
    //     score2 = chain_score(chain,anchors, mems, min_score,(not (top_k > 1)),-1,read,0,out);
    //   }
    // }

    // // Report the high-score chain
    // if(top4[0].second > 1)
    // {
    //   auto chain = chains_rev[top4[0].second-2];
    //   // Reverse the chain order
    //   std::reverse(chain.second.begin(), chain.second.end());
    //   // Compute the score of a chain.
    //   score = chain_score(chain,anchors_rev, mems_rev, min_score,false,score2,&read_rev,16,out);
    // }else{
    //   auto chain = chains[top4[0].second];
    //   // Reverse the chain order
    //   std::reverse(chain.second.begin(), chain.second.end());
    //   // Compute the score of a chain.
    //   score = chain_score(chain,anchors, mems, min_score,false,score2,read,0,out);
      
    // }
    // if(score > min_score) aligned = true;
    



    // // Compute the alignment score for the top_k chains
    // // NOTE: for seconadry alignments the MAPQ score is set to 255 (as in Bowtie2)
    // // QUESTION: Should I compute the scores for all the chains?
    // size_t j = 0;
    // size_t j_rev = 0;
    // for(size_t i = 0; i < min(chains.size() + chains_rev.size(),top_k); ++i)
    // {
    //   // QUESTION: Should I compute the score for the top_k distinct alignments?

    //   if(j < chains.size() and j_rev < chains_rev.size())
    //   {
    //     if(chains[j].first > chains_rev[j_rev].first)
    //     {
    //       if(i > 1){
    //         auto chain = chains[j];
    //         // Reverse the chain order
    //         std::reverse(chain.second.begin(), chain.second.end());
    //         // Compute the score of a chain.
    //         chain_score(chain,anchors, mems, min_score,false,-1,read,0,out);
    //       }
    //       j++;
    //     }
    //     else
    //     {
    //       if(i > 1){
    //         auto chain = chains_rev[j_rev];
    //         // Reverse the chain order
    //         std::reverse(chain.second.begin(), chain.second.end());
    //         // Compute the score of a chain.
    //         score2 = chain_score(chain,anchors_rev, mems_rev, min_score,false,-1,&read_rev,16,out);
    //       }
    //       j_rev++;
    //     }
    //   }
    //   else if (j < chains.size())
    //   {
    //     if(i > 1)
    //     {
    //         auto chain = chains[j];
    //         // Reverse the chain order
    //         std::reverse(chain.second.begin(), chain.second.end());
    //         // Compute the score of a chain.
    //         chain_score(chain,anchors, mems, min_score,false,-1,read,0,out);
    //     }
    //     j++;
    //   }
    //   else if (j_rev < chains_rev.size())
    //   {
    //     if(i > 1)
    //     {
    //       auto chain = chains_rev[j_rev];
    //       // Reverse the chain order
    //       std::reverse(chain.second.begin(), chain.second.end());
    //       // Compute the score of a chain.
    //       chain_score(chain,anchors_rev, mems_rev, min_score,false,-1,&read_rev,16,out);
          
    //     }
    //     j_rev++;
    //   }

    //   // // Extract the anchors
    //   // std::vector<std::pair<size_t,size_t>> chain_anchors(chain.second.size());
    //   // for(size_t i = 0; i < chain_anchors.size(); ++i)
    //   //   chain_anchors[i] = anchors[chain.second[i]];
    //   // // Extracting left and right context of the read
    //   // // lcs: left context sequence
    //   // size_t lcs_len = mems[chain_anchors[0].first].idx;
    //   // uint8_t* lcs = (uint8_t*)malloc(lcs_len);
    //   // // verbose("lcs: " + std::string(read->seq.s).substr(0,lcs_len));
    //   // // Convert A,C,G,T,N into 0,1,2,3,4
    //   // // The left context is reversed
    //   // for (size_t i = 0; i < lcs_len; ++i)
    //   //   lcs[lcs_len -i -1] = seq_nt4_table[(int)read->seq.s[i]];

    //   // // rcs: right context sequence
    //   // size_t rcs_occ = (mems[chain_anchors.back().first].idx + mems[chain_anchors.back().first].len); // The first character of the right context
    //   // size_t rcs_len = read->seq.l - rcs_occ;
    //   // uint8_t* rcs = (uint8_t*)malloc(rcs_len);
    //   // // verbose("rcs: " + std::string(read->seq.s).substr(rcs_occ,rcs_len));
    //   // // Convert A,C,G,T,N into 0,1,2,3,4
    //   // for (size_t i = 0; i < rcs_len; ++i) 
    //   //   rcs[i] = seq_nt4_table[(int)read->seq.s[ rcs_occ + i]];

    //   // int32_t min_score = 20 + 8 * log(read->seq.l);
    //   // // Fill between MEMs
    //   // int32_t score = fill_chain(
    //   //     mems,
    //   //     chain_anchors,
    //   //     lcs,   // Left context of the read
    //   //     lcs_len, // Left context of the read lngth
    //   //     rcs,   // Right context of the read
    //   //     rcs_len, // Right context of the read length
    //   //     read
    //   //   );

    //   //   if(score > min_score)
    //   //   {
    //   //     fill_chain(mems,chain_anchors,lcs,lcs_len,rcs,rcs_len,read,false,0,strand,out); 
    //   //     aligned = true;
    //   //   }

    //   //   delete rcs;
    //   //   delete lcs;
      
    // }

    /****************************************************************************/
    /**** TODO: Tailored dynamic programming
    // Sort MEM occurrences
    for(auto mem: mems)
      std::sort(mem.occs.begin(), mem.occs.end());

    // Find non-overlapping MEM chains
    std::vector<std::vector<size_t>> chains;
    // first: chain idx
    // second: mems idx
    std::queue<std::pair<size_t,size_t>> next;
    for(size_t i = 0; i < mems.size(); ++i)
    {
      chains.push_back(std::vector<size_t>());
      next.push(make_pair(i,i)); 
    }

    while(not next.empty())
    {
      auto e = next.front(); next.pop();
      chains[e.first].push_back(e.second);
      // Find next non-overlapping MEM
      size_t mem_end = mems[e.second].idx + mems[e.second].len;
      size_t i = e.second + 1;
      while( i < mems.size())
      {
        if(mems[i].idx >= mem_end)
        {
          chains.push_back(std::vector<size_t>(chains[e.first]));
          next.push(make_pair(chains.size()-1,i));
        }
        ++i;
      }
    }
    // Debug
    verbose("New chains: ",chains.size());
    for(auto chain: chains)
    {
      for(auto mem: chain)
        verbose(mems[mem].idx, mems[mem].len, mems[mem].occs.size());
      verbose("new chain");
    }

    // Reverse chains order
    std::reverse(chains.begin(),chains.end());

    // Find compatible occurrences of the MEMs in the chain
    size_t max_gap = 10;

    for(size_t i = 0; i< chains.size(); ++i)
    {
      if(chains[i].size() < 2)
        continue;

      auto& chain = chains[i];

      std::vector<size_t> j(chain.size(),0); // Index of the occurrences

      std:vector<std::vector<size_t>> valid();
      // Lambda helpers
      // Compute the distance between the i-th and i+1-th MEM in the chain
      auto mem_distance = [] (size_t i) -> size_t {
        return mems[chain[i+1]].idx - mems[chain[i]].idx;
      };
      // Compute the distance between the i-th and i+1-th MEM occurrences in j
      auto occ_distance = [] (size_t i) -> size_t {
        return mems[chain[i+1]].occs[j[i+1]] - mems[chain[i]].occs[j[i]];
      };
      // Increment j indexes
      auto inc_j = [&] () -> bool {
        
      };

      

    }
    ********/




    // The distance between the MEMs in the read
    

    // // Align the read
    // if (mem_len >= min_len)
    // {
    //   // // Extract all the occurrences of the MEM
    //   // std::vector<size_t> occs;
    //   // occs.push_back(mem_pos);

    //   // // Phi direction
    //   // size_t curr = mem_pos;
    //   // size_t next = ms.Phi(curr);
    //   // size_t lcp =  lceToRBounded(ra,curr,next,mem_len);
    //   // while(lcp >= mem_len)
    //   // {
    //   //   occs.push_back(next);
        
    //   //   curr = next;
    //   //   next = ms.Phi(curr);
    //   //   lcp  = lceToRBounded(ra,curr,next,mem_len);
    //   //   // verbose("Phi: " + std::to_string(lcp));
    //   //   // if(occs.size() > 100)
    //   //   //   error("More than 100 occs Phi" + std::string(read->seq.s));
    //   // }

    //   // // Phi_inv direction
    //   // curr = mem_pos;
    //   // next = ms.Phi_inv(curr);
    //   // lcp =  lceToRBounded(ra,curr,next,mem_len);
    //   // while(lcp >= mem_len)
    //   // {
    //   //   occs.push_back(next);
        
    //   //   curr = next;
    //   //   next = ms.Phi_inv(curr);
    //   //   lcp  = lceToRBounded(ra,curr,next,mem_len);
    //   //   // verbose("Phi_inv: " + std::to_string(next));
    //   //   // if(occs.size() > 100)
    //   //   //   error("More than 100 occs Phi_inv" + std::string(read->seq.s));
    //   // }

    //   // Extractin left and right context of the read
    //   // lcs: left context sequence
    //   size_t lcs_len = mem_idx;
    //   uint8_t* lcs = (uint8_t*)malloc(lcs_len);
    //   // verbose("lcs: " + std::string(read->seq.s).substr(0,lcs_len));
    //   // Convert A,C,G,T,N into 0,1,2,3,4
    //   // The left context is reversed
    //   for (size_t i = 0; i < lcs_len; ++i)
    //     lcs[lcs_len -i -1] = seq_nt4_table[(int)read->seq.s[i]];

    //   // rcs: right context sequence
    //   size_t rcs_occ = (mem_idx + mem_len); // The first character of the right context
    //   size_t rcs_len = read->seq.l - rcs_occ;
    //   uint8_t* rcs = (uint8_t*)malloc(rcs_len);
    //   // verbose("rcs: " + std::string(read->seq.s).substr(rcs_occ,rcs_len));
    //   // Convert A,C,G,T,N into 0,1,2,3,4
    //   for (size_t i = 0; i < rcs_len; ++i) 
    //     rcs[i] = seq_nt4_table[(int)read->seq.s[ rcs_occ + i]];

    //   int32_t min_score = 20 + 8 * log(read->seq.l);
    //   // verbose("Number of occurrences: " + std::to_string(occs.size()));
    //   // For all the occurrences align
    //   for(auto curr_mem_pos: occs)
    //   {
    //     int32_t score = extend(
    //       curr_mem_pos,
    //       mem_len,
    //       lcs,   // Left context of the read
    //       lcs_len, // Left context of the read lngth
    //       rcs,   // Right context of the read
    //       rcs_len // Right context of the read length
    //     );

    //     if(score > min_score)
    //     {
    //       extend(curr_mem_pos,mem_len,lcs,lcs_len,rcs,rcs_len,false,0,min_score,read,strand,out); 
    //       aligned = true;
    //     }

    //   }
    //   delete lcs;
    //   delete rcs;
    // }

    if(not aligned)
    {
      std::string dummy = "";
      write_sam(0, 0, 0, 0, "*", read, 0, out, dummy, dummy, 0, "*", 0, 0);
    }

    return aligned;
  }

  typedef struct{
    bool aligned = false;
    enum type_t {Unpaired, Paired} type;

    kseq_t* read;
    kseq_t* mate1;
    kseq_t* mate2;




  } alignment_t;

  // Aligning pair-ended sequences
  bool align(kseq_t *mate1, kseq_t *mate2, FILE *out)
  {

    bool aligned = false;
    
    MTIME_START(0); // Timing helper

    // Generate rc reads
    kseq_t mate1_rev, mate2_rev;
    rc_copy_kseq_t(&mate1_rev, mate1);
    rc_copy_kseq_t(&mate2_rev, mate2);

    // Find MEMs
    std::vector<mem_t> mems;

    find_mems(mate1, mems, 0, MATE_1 | MATE_F);
    find_mems(&mate1_rev, mems, mate2->seq.l, MATE_1 | MATE_RC );
    find_mems(mate2, mems, 0, MATE_2 | MATE_F);
    find_mems(&mate2_rev, mems, mate1->seq.l, MATE_2 | MATE_RC);

    MTIME_END(0); //Timing helper
    MTIME_START(1); //Timing helper

    std::vector<std::pair<size_t, size_t>> anchors;
    std::vector<chain_t> chains;

    // Chain MEMs
    const bool chained = find_chains(mems, anchors, chains);

    MTIME_END(1);   //Timing helper
    MTIME_START(2); //Timing helper

    if (not chained)
    {
      sam_t sam_m1(mate1);
      sam_t sam_m2(mate2);

      // Fill sam fields RNEXT, PNEXT and TLEN
      sam_m1.rnext = std::string(mate2->name.s);
      sam_m2.rnext = std::string(mate1->name.s);

      sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;
      
      write_sam(out, sam_m1);
      write_sam(out, sam_m2);
      return false;
    }

    int32_t min_score = 20 + 8 * log(mate1->seq.l);

    // // Compute the second best score
    std::vector<std::pair<int32_t, size_t>> best_scores;
    // get the occurrences of the top 4 best scores
    std::set<size_t> different_scores;
    size_t i = 0;
    while (i < chains.size() and different_scores.size() < 3)
    {
        different_scores.insert(chains[i].score);
        if (different_scores.size() < 3)
        {
          // Align the chain
          auto chain = chains[i];
          // Reverse the chain order
          std::reverse(chain.anchors.begin(), chain.anchors.end());
          // Compute the score of a chain.
          paired_score_t score;
          if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
            score = paired_chain_score(chain, anchors, mems, min_score, mate1, &mate2_rev);
          else
            score = paired_chain_score(chain, anchors, mems, min_score, &mate1_rev, mate2);
          
          best_scores.push_back(std::make_pair(score.tot, i++));
        }
    }
    
    if (best_scores.size() < 2)
      best_scores.push_back(std::make_pair(0, chains.size()));

    std::sort(best_scores.begin(), best_scores.end(), std::greater<std::pair<int32_t, size_t>>());

    assert(best_scores.size() > 1);

    if (best_scores[0].first < min_score)
    {
      sam_t sam_m1(mate1);
      sam_t sam_m2(mate2);

      // Fill sam fields RNEXT, PNEXT and TLEN
      sam_m1.rnext = std::string(mate2->name.s);
      sam_m2.rnext = std::string(mate1->name.s);

      sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;

      write_sam(out, sam_m1);
      write_sam(out, sam_m2);
      return false;
    }

    int32_t score2 = best_scores[1].first;
    paired_score_t score;


    { // Forward case
      i = best_scores[0].second;
      // Align the chain
      auto chain = chains[i];
      // Reverse the chain order
      std::reverse(chain.anchors.begin(), chain.anchors.end());
      // Compute the score of a chain.
      if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
        score = paired_chain_score(chain, anchors, mems, min_score, mate1, &mate2_rev, false, score2, 0, out);
      else
        score = paired_chain_score(chain, anchors, mems, min_score, &mate1_rev, mate2, false, score2, 1, out);
    }

    if (score.tot >= min_score)
      aligned = true;

    if (not aligned)
    {
      sam_t sam_m1(mate1);
      sam_t sam_m2(mate2);

      // Fill sam fields RNEXT, PNEXT and TLEN
      sam_m1.rnext = std::string(mate2->name.s);
      sam_m2.rnext = std::string(mate1->name.s);

      sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;

      write_sam(out, sam_m1);
      write_sam(out, sam_m2);
      return false;
    }

    MTIME_END(2); //Timing helper

    return aligned;
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

  int32_t chain_score(
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
    int32_t score = fill_chain(
        mems,
        chain_anchors,
        lcs,   // Left context of the read
        lcs_len, // Left context of the read lngth
        rcs,   // Right context of the read
        rcs_len, // Right context of the read length
        read
      );

      if(!score_only and score >= min_score)
      {
        bool output = (sam == nullptr);

        if(output)
          sam = new sam_t(read);

        fill_chain(mems,chain_anchors,lcs,lcs_len,rcs,rcs_len,read,score_only,sam); 
        
        sam->flag = (strand?16:0);
        sam->zs = score2;
        sam->mapq = compute_mapq(sam->as, sam->zs, min_score, sam->read->seq.l * smatch);
        // sam->reverse = (strand != 0);
        sam->rname = idx[sam->pos - 1];
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

  typedef struct{
    int32_t tot = 0;
    int32_t m1 = 0;
    int32_t m2 = 0; 
    bool paired = false;
  } paired_score_t;

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
      FILE *out = nullptr)
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
      }
      else
      {
        sam_t sam_m1(mate1);
        sam_t sam_m2(mate2);


        score.m1 = chain_score(mate1_chain, anchors, mems, min_score, mate1, false, score2, strand, nullptr, &sam_m1);
        score.m2 = chain_score(mate2_chain, anchors, mems, min_score, mate2, false, score2, strand, nullptr, &sam_m2);

        // Fill sam fields RNEXT, PNEXT and TLEN
        sam_m1.rnext = std::string(mate2->name.s);
        sam_m2.rnext = std::string(mate1->name.s);

        if(score.m1 >= min_score and score.m2 >= min_score)
        {

          sam_m1.pnext = sam_m2.pos;
          sam_m2.pnext = sam_m1.pos;

          ll tlen = (sam_m2.pos + mate2->seq.l) - sam_m1.pos;

          sam_m1.tlen = tlen;
          sam_m2.tlen = -tlen;

          sam_m1.rname = idx[sam_m1.pos - 1];
          sam_m2.rname = idx[sam_m2.pos - 1];

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
        }else if(score.m1 >= min_score) {
          sam_m1.rname = idx[sam_m1.pos - 1];

          sam_m1.flag = SAM_PAIRED | SAM_MATE_UNMAPPED | SAM_FIRST_IN_PAIR;
          sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_SECOND_IN_PAIR;
          if(strand)
            sam_m1.flag |= SAM_REVERSED;
        }else if(score.m2 >= min_score) {
          sam_m2.rname = idx[sam_m2.pos - 1];

          sam_m1.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_FIRST_IN_PAIR;
          sam_m2.flag = SAM_PAIRED | SAM_MATE_UNMAPPED | SAM_SECOND_IN_PAIR;
          if(not strand)
            sam_m2.flag |= SAM_REVERSED;
        }else {
          sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;
        }

        write_sam(out,sam_m1);
        write_sam(out,sam_m2);
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


    score.tot = score.m1 + score.m2;
    return score;
  }

  size_t get_aligned_reads()
  {
    return aligned_reads;
  }

  // If score_only is true we compute the score of the alignment. 
  // If score_only is false, we extend again the read and we write the result
  // in the SAM file, so we need to give the second best score.
  int32_t extend(
      const size_t mem_pos,
      const size_t mem_len,
      const uint8_t *lcs,           // Left context of the read
      const size_t lcs_len,         // Left context of the read lngth
      const uint8_t *rcs,           // Right context of the read
      const size_t rcs_len,         // Right context of the read length
      const bool score_only = true, // Report only the score
      const int32_t score2 = 0,     // The score of the second best alignment
      const int32_t min_score = 0,  // The minimum score to call an alignment
      const kseq_t *read = nullptr, // The read that has been aligned
      int8_t strand = 0,            // 0: forward aligned ; 1: reverse complement aligned
      FILE *out = nullptr,          // The SAM file pointer
      const bool realign = false    // Realign globally the read
  )
  {
    int flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT;

    if(score_only) 
      flag = KSW_EZ_SCORE_ONLY;

    int score_lc = 0;
    int score_rc = 0;

    // TODO: Update end_bonus according to the MEM contribution to the score
    ksw_extz_t ez_lc;
    ksw_extz_t ez_rc;
    ksw_extz_t ez;
    memset(&ez_lc, 0, sizeof(ksw_extz_t));
    memset(&ez_rc, 0, sizeof(ksw_extz_t));
    memset(&ez, 0, sizeof(ksw_extz_t));

    // Extract the context from the reference
    // lc: left context
    ksw_reset_extz(&ez_lc);
    if(lcs_len > 0)
    {
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
      // Check if the extension reached the end or the query
      assert(score_only or ez_lc.reach_end);

      // std::string blc = print_BLAST_like((uint8_t*)lc,(uint8_t*)lcs,ez_lc.cigar,ez_lc.n_cigar);
      // std::cout<<blc;

      delete lc;
    }

    // rc: right context
    ksw_reset_extz(&ez_rc);
    if(rcs_len > 0)
    {
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
      // Check if the extension reached the end or the query
      assert(score_only or ez_rc.reach_end);

      // std::string brc = print_BLAST_like((uint8_t*)rc,(uint8_t*)rcs,ez_rc.cigar,ez_rc.n_cigar);
      // std::cout<<brc;
      delete rc;
    }


    // Compute the final score
    int32_t score = mem_len * smatch + score_lc + score_rc;

    if(not score_only)
    {
      // Compute starting position in reference
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

      char* tmp = (char*)calloc(max(ref_len,seq_len),1);

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

        assert(ez.score >= score);

        // Concatenate the CIGAR strings
        
        std::string cigar_s;
        for(size_t i = 0; i < ez.n_cigar; ++i)
          cigar_s += std::to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];

        // Compute the MD:Z field and thenumber of mismatches
        std::string mdz_s;
        size_t nm = write_MD_core((uint8_t*)ref,seq,ez.cigar,ez.n_cigar,tmp,0,mdz_s);

        write_sam(ez.score, score2, min_score, ref_pos, idx[ref_pos].c_str(), read, strand, out, cigar_s, mdz_s, nm, "*", 0, 0);
      }
      else
      {
        // Concatenate the CIGAR strings
        size_t n_cigar = ez_lc.n_cigar + ez_rc.n_cigar + 1;
        uint32_t *cigar = (uint32_t*)calloc(n_cigar,sizeof(uint32_t));
        size_t i = 0;

        for(size_t j = 0; j < ez_lc.n_cigar; ++j)
          cigar[i++] = ez_lc.cigar[ez_lc.n_cigar -j -1];


        if(ez_lc.n_cigar > 0 and ((cigar[i-1]& 0xf) == 0))
        { // If the previous operation is also an M then merge the two operations
          cigar[i-1] += (((uint32_t)mem_len) << 4);
          --n_cigar;
        }
        else
          cigar[i++] = (((uint32_t)mem_len) << 4);


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

        // std::string bfull = print_BLAST_like((uint8_t*)ref,seq,cigar,n_cigar);
        // std::cout << bfull;

        std::string cigar_s;
        for(size_t i = 0; i < n_cigar; ++i)
          cigar_s += std::to_string(cigar[i] >> 4) + "MID"[cigar[i] & 0xf];


        // Compute the MD:Z field and thenumber of mismatches
        std::string mdz_s;
        size_t nm = write_MD_core((uint8_t*)ref,seq,cigar,n_cigar,tmp,0,mdz_s);

        write_sam(score, score2, min_score, ref_pos, idx[ref_pos].c_str(), read, strand, out, cigar_s, mdz_s, nm, "*", 0, 0);

        delete cigar;
      }
      delete tmp;
      delete ref;
      delete seq;
    }

    if (ez_lc.m_cigar > 0)
      delete ez_lc.cigar;
    if (ez_rc.m_cigar > 0)
      delete ez_rc.cigar;
    if (ez.m_cigar > 0)
      delete ez.cigar;

    return score;
  }

  // If score_only is true we compute the score of the alignment. 
  // If score_only is false, we extend again the read and we write the result
  // in the SAM file, so we need to give the second best score. 
  int32_t fill_chain(
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
      // Check if the extension reached the end or the query
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
      // Check if the extension reached the end or the query
      assert(score_only or ez_rc.reach_end);

      // std::string brc = print_BLAST_like((uint8_t*)rc,(uint8_t*)rcs,ez_rc.cigar,ez_rc.n_cigar);
      // std::cout<<brc;
      delete rc;
    }


    // Compute the partial score score
    int32_t score = score_lc + score_rc;
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

    // TODO: fix the gap fill when MEMs overlap. The quick fix is a full alignment between the first and last MEM
    // // Fill the gaps between each mem
    // std::vector<ksw_extz_t> ez_cc(anchors.size()-1);
    // for(size_t i = 1; i < anchors.size(); ++i)
    // {
    //   flag = KSW_EZ_RIGHT;
    //   ksw_reset_extz(&ez_cc[i-1]);

    //   // size_t cc_occ = mems[anchors[i-1].first].occs[anchors[i-1].second] + mems[anchors[i-1].first].len;
    //   // size_t cc_len = mems[anchors[i].first].occs[anchors[i].second] - cc_occ; 
    //   // char *cc = (char *)malloc(cc_len);
    //   // ra.expandSubstr(cc_occ, cc_len, cc); // TODO: Replace all expand substrings with only one expand at the beginning.
    //   // // verbose("rc: " + std::string(rc));
    //   // // Convert A,C,G,T,N into 0,1,2,3,4
    //   // for (size_t i = 0; i < cc_len; ++i)
    //   //   cc[i] = seq_nt4_table[(int)cc[i]];

    //   size_t cc_occ = mems[anchors[i-1].first].occs[anchors[i-1].second] + mems[anchors[i-1].first].len;
    //   size_t cc_len = mems[anchors[i].first].occs[anchors[i].second] - cc_occ; 
    //   cc_occ -= ref_pos;
    //   char *cc = (char *)malloc(cc_len);
    //   // Convert A,C,G,T,N into 0,1,2,3,4
    //   for (size_t i = 0; i < cc_len; ++i)
    //     cc[i] = ref[cc_occ + i];


    //   size_t ccs_pos = mems[anchors[i-1].first].idx + mems[anchors[i-1].first].len;
    //   size_t ccs_len = mems[anchors[i].first].idx - ccs_pos;

    //   // Query: rcs
    //   // Target: rc
    //   // verbose("aligning rc and rcs");
    //   ksw_extz2_sse(km, ccs_len, (uint8_t*)(seq + ccs_pos), cc_len, (uint8_t*)cc, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez_cc[i-1]);
    //   // score_cc = ez_cc.mqe;
    //   // verbose("rc score: " + std::to_string(score_rc));
    //   // Check if the extension reached the end or the query
    //   // assert(score_only or ez_cc[i-1].reach_end);

    //   // Update score
    //   score +=  mems[anchors[i-1].first].len * smatch + ez_cc[i-1].mqe;

    //   // std::string brc = print_BLAST_like((uint8_t*)rc,(uint8_t*)rcs,ez_rc.cigar,ez_rc.n_cigar);
    //   // std::cout<<brc;
    //   delete cc;
    // }

    // score +=  mems[anchors.back().first].len * smatch;

//******************************************************************************
// BEGIN RAPID PROTOTYPING HACK
//******************************************************************************
    // Realign the whole sequence globally
    ksw_reset_extz(&ez);
    ksw_extz2_sse(km, seq_len, (uint8_t *)seq, ref_len, (uint8_t *)ref, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez);

    score = ez.score;

    realign = true;

//******************************************************************************
// END RAPID PROTOTYPING HACK
//******************************************************************************
    if(not score_only)
    {
      // Compute starting position in reference

      char* tmp = (char*)calloc(max(ref_len,seq_len),1);

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

        assert(ez.score >= score);

        // Concatenate the CIGAR strings
        
        // std::string cigar_s;
        sam->cigar = "";
        for(size_t i = 0; i < ez.n_cigar; ++i)
          sam->cigar += std::to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];

        // Compute the MD:Z field and the number of mismatches
        sam->nm = write_MD_core((uint8_t*)ref,seq,ez.cigar,ez.n_cigar,tmp,0,sam->md);

        sam->as = ez.score;
        sam->pos = ref_pos + 1; // ref_pos is 1 based
        sam->rname = idx[ref_pos];
        sam->rlen = ref_len;
        // write_sam(ez.score, score2, min_score, ref_pos, "human", read, strand, out, cigar_s, mdz_s, nm, "*", 0, 0);
      }
      else
      {
        // // Concatenate the CIGAR strings
        // size_t n_cigar = ez_lc.n_cigar + ez_rc.n_cigar + 1;
        // for(size_t i = 0; i < anchors.size()-1; ++i)
        //   n_cigar += ez_cc[i].n_cigar + 1;

        // uint32_t *cigar = (uint32_t*)calloc(n_cigar,sizeof(uint32_t));
        // size_t i = 0;

        // for(size_t j = 0; j < ez_lc.n_cigar; ++j)
        //   cigar[i++] = ez_lc.cigar[ez_lc.n_cigar -j -1];

        // // CIGARs between MEMs
        // for(size_t j = 0; j < anchors.size(); ++j)
        // {
        //   // CIGAR for the MEM
        //   size_t mem_len = mems[anchors[j].first].len;
        //   if(i > 0 and ((cigar[i-1]& 0xf) == 0))
        //   { // If the previous operation is also an M then merge the two operations
        //     cigar[i-1] += (((uint32_t)mem_len) << 4);
        //     --n_cigar;
        //   }
        //   else
        //     cigar[i++] = (((uint32_t)mem_len) << 4);

        //   // CIGAR of the gap between the j-th MEM and the j+1-th MEM
        //   if(j < anchors.size() - 1)
        //   {
        //     if(ez_cc[j].n_cigar > 0)
        //     {
        //       // Check the first element
        //       if((ez_cc[j].cigar[0]& 0xf) == 0)
        //       { // If the next operation is also an M then merge the two operations
        //         cigar[i-1] += ez_cc[j].cigar[0];
        //         --n_cigar;
        //       }
        //       else
        //         cigar[i++] = ez_cc[j].cigar[0];
        //     } 
                      
        //       // Copy all the other elements
        //       for(size_t k = 1; k < ez_cc[j].n_cigar; ++k)
        //         cigar[i++] = ez_cc[j].cigar[k];
        //   }
        // }

        // // CIGAR of the right context
        // if(ez_rc.n_cigar > 0)
        // {
        //   if((ez_rc.cigar[0]& 0xf) == 0)
        //   { // If the next operation is also an M then merge the two operations
        //     cigar[i-1] += ez_rc.cigar[0];
        //     --n_cigar;
        //   }
        //   else
        //     cigar[i++] = ez_rc.cigar[0];
        // }

        // for(size_t j = 1; j < ez_rc.n_cigar; ++j)
        //   cigar[i++] = ez_rc.cigar[j];
        
        // assert(i <= n_cigar);

        // // std::string bfull = print_BLAST_like((uint8_t*)ref,seq,cigar,n_cigar);
        // // std::cout << bfull;

        // std::string cigar_s;
        // for(size_t i = 0; i < n_cigar; ++i)
        //   cigar_s += std::to_string(cigar[i] >> 4) + "MID"[cigar[i] & 0xf];


        // // Compute the MD:Z field and thenumber of mismatches
        // std::pair<std::string,size_t> md_nm = write_MD_core((uint8_t*)ref,seq,cigar,n_cigar,tmp,0);

        // write_sam(score,score2,min_score,ref_pos,"human",read,strand,out,cigar_s,md_nm.first,md_nm.second);

        // delete cigar;
      }
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
      std::string res = "@HD VN:1.6 SO:unknown\n";
      res += idx.to_sam();
      res += "@PG ID:moni PN:moni VN:0.1.0\n";
      return res; 
  }

protected:
    ms_t ms;
    slp_t ra;
    seqidx idx;
    // SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd> ra;

    const size_t min_len = 0;
    const size_t ext_len = 100;   // Extension length
    size_t aligned_reads = 0;
    size_t n = 0;
    size_t top_k = 1; // report the top_k alignments

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


    //   // From https://github.com/BenLangmead/bowtie2/blob/4512b199768e562e8627ffdfd9253affc96f6fc6/unique.cpp
    //   // There is no valid second-best alignment and the best alignment has a
    //   // perfect score.
    //   const uint32_t unp_nosec_perf = 44;

    //   // There is no valid second-best alignment.  We stratify the alignment
    //   // score of the best alignment into 10 bins.
    //   const uint32_t unp_nosec[11] = {
    //     43, 42, 41, 36, 32, 27, 20, 11, 4, 1, 0
    //   };

    //   // The best alignment has a perfect score, and we stratify the distance
    //   // between best and second-best alignment scores into 10 bins.
    //   const uint32_t unp_sec_perf[11] = {
    //     2, 16, 23, 30, 31, 32, 34, 36, 38, 40, 42
    //   };

    //   // The best alignment has a non-perfect score, and we stratify both by best
    //   // alignment score (specifically, the maximum score minus the best "best")
    //   // and by the distance between the best and second-best alignment scores
    //   // ("difference").  Each is stratified into 10 bins.  Each row is a
    //   // difference (smaller elts = smaller differences) and each column is a
    //   // best score (smaller elts = higher best alignment scores).
    //   const uint32_t unp_sec[11][11] = {
    //     {  2,  2,  2,  1,  1, 0, 0, 0, 0, 0, 0},
    //     { 20, 14,  7,  3,  2, 1, 0, 0, 0, 0, 0},
    //     { 20, 16, 10,  6,  3, 1, 0, 0, 0, 0, 0},
    //     { 20, 17, 13,  9,  3, 1, 1, 0, 0, 0, 0},
    //     { 21, 19, 15,  9,  5, 2, 2, 0, 0, 0, 0},
    //     { 22, 21, 16, 11, 10, 5, 0, 0, 0, 0, 0},
    //     { 23, 22, 19, 16, 11, 0, 0, 0, 0, 0, 0},
    //     { 24, 25, 21, 30,  0, 0, 0, 0, 0, 0, 0},
    //     { 30, 26, 29,  0,  0, 0, 0, 0, 0, 0, 0},
    //     { 30, 27,  0,  0,  0, 0, 0, 0, 0, 0, 0},
    //     { 30,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0},
    //   };

    //   //
    //   // Paired mapping quality:
    //   //

    //   // There is no valid second-best alignment and the best alignment has a
    //   // perfect score.
    //   const uint32_t pair_nosec_perf = 44;
};

#endif /* end of include guard: _ALIGNER_KSW2_HH */
