/* align - Align the reads to the reference
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
   \file align.cpp
   \brief align.cpp Align the reads to the reference.
   \author Massimiliano Rossi
   \date 13/07/2020
*/

extern "C" {
#include <xerrors.h>
}

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/io.hpp>

#include <moni.hpp>

#include <malloc_count.h>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>
#include <PlainSlp.hpp>
#include <FixedBitLenCode.hpp>

#include <ksw2.h>

#include <omp.h>

#include <libgen.h>

#define _REALIGN

// KSEQ_INIT(gzFile, gzread);

//*************** Borrowed from minimap2 ***************************************
static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
	uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}
//******************************************************************************
////////////////////////////////////////////////////////////////////////////////
/// kseq extra
////////////////////////////////////////////////////////////////////////////////

static inline size_t ks_tell(kseq_t *seq)
{
  return gztell(seq->f->f) - seq->f->end + seq->f->begin;
}

void copy_kstring_t(kstring_t &l, kstring_t &r)
{
  l.l = r.l;
  l.m = r.m;
  l.s = (char *)malloc(l.m);
  for (size_t i = 0; i < r.m; ++i)
    l.s[i] = r.s[i];
}
void copy_kseq_t(kseq_t *l, kseq_t *r)
{
  copy_kstring_t(l->name, r->name);
  copy_kstring_t(l->comment, r->comment);
  copy_kstring_t(l->seq, r->seq);
  copy_kstring_t(l->qual, r->qual);
  l->last_char = r->last_char;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// helper functions
////////////////////////////////////////////////////////////////////////////////

char complement(char n)
{
  switch (n)
  {
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'G':
    return 'C';
  case 'C':
    return 'G';
  default:
    return n;
  }
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// SLP definitions
////////////////////////////////////////////////////////////////////////////////

using SelSd = SelectSdvec<>;
using DagcSd = DirectAccessibleGammaCode<SelSd>;
using Fblc = FixedBitLenCode<>;

using shaped_slp_t = SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd>;
using plain_slp_t = PlainSlp<uint32_t, Fblc, Fblc>;

template<typename slp_t>
std::string get_slp_file_extension()
{
  return std::string(".slp");
}

template <>
std::string get_slp_file_extension<shaped_slp_t>()
{
  return std::string(".slp");
}

template <>
std::string get_slp_file_extension<plain_slp_t>()
{
  return std::string(".plain.slp");
}
////////////////////////////////////////////////////////////////////////////////

template<typename slp_t>
class aligner
{
public:
  // using SelSd = SelectSdvec<>;
  // using DagcSd = DirectAccessibleGammaCode<SelSd>;
  // using SlpT = SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd>;
  using ll = long long int;

  aligner(std::string filename, 
            size_t min_len_ = 50, 
            bool forward_only_ = true): 
                min_len(min_len_), 
                forward_only(forward_only_)
  {
    verbose("Loading the matching statistics index");
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    std::string filename_ms = filename + ms.get_file_extension();

    ifstream fs_ms(filename_ms);
    ms.load(fs_ms);
    fs_ms.close();

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Matching statistics index construction complete");
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

    verbose("Matching statistics index loading complete");
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

    memset(&ez_lc, 0, sizeof(ksw_extz_t));
    memset(&ez_rc, 0, sizeof(ksw_extz_t));
    memset(&ez, 0, sizeof(ksw_extz_t));
  }

  // Destructor
  ~aligner() 
  {
    if(ez_lc.m_cigar > 0)
      delete ez_lc.cigar;
    if(ez_rc.m_cigar > 0)
      delete ez_rc.cigar;
    if(ez.m_cigar > 0)
      delete ez.cigar;
      // NtD
  }

  bool align(kseq_t *read, FILE* out, uint8_t strand, bool mem_chaining = true)
  {
    if(mem_chaining) // TODO: fix the no chaining alignment
      return (strand == 0? align_chain(read,out):false);
    else
      return align_no_chain(read,out,strand);
  }

  // bool align(kseq_t *read, FILE* out, uint8_t strand)
  bool align_no_chain(kseq_t *read, FILE* out, uint8_t strand)
  {
    size_t mem_pos = 0;
    size_t mem_len = 0;
    size_t mem_idx = 0;

    bool aligned = false;

    // Find MEMs.
    auto pointers = ms.query(read->seq.s, read->seq.l);
    std::vector<size_t> lengths(pointers.size());
    size_t l = 0;
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

      lengths[i] = l;
      l = (l == 0 ? 0 : (l - 1));

      // Update MEM
      if (lengths[i] > mem_len and n_Ns < lengths[i])
      {
        mem_len = lengths[i];
        mem_pos = pointers[i];
        mem_idx = i;
      }
    }

    // Align the read
    if (mem_len >= min_len)
    {
      // Extract all the occurrences of the MEM
      std::vector<size_t> occs;
      occs.push_back(mem_pos);

      // Phi direction
      size_t curr = mem_pos;
      size_t next = ms.Phi(curr);
      size_t lcp =  lceToRBounded(ra,curr,next,mem_len);
      while(lcp >= mem_len)
      {
        occs.push_back(next);
        
        curr = next;
        next = ms.Phi(curr);
        lcp  = lceToRBounded(ra,curr,next,mem_len);
        // verbose("Phi: " + std::to_string(lcp));
        // if(occs.size() > 100)
        //   error("More than 100 occs Phi" + std::string(read->seq.s));
      }

      // Phi_inv direction
      curr = mem_pos;
      next = ms.Phi_inv(curr);
      lcp =  lceToRBounded(ra,curr,next,mem_len);
      while(lcp >= mem_len)
      {
        occs.push_back(next);
        
        curr = next;
        next = ms.Phi_inv(curr);
        lcp  = lceToRBounded(ra,curr,next,mem_len);
        // verbose("Phi_inv: " + std::to_string(next));
        // if(occs.size() > 100)
        //   error("More than 100 occs Phi_inv" + std::string(read->seq.s));
      }

      // Extractin left and right context of the read
      // lcs: left context sequence
      size_t lcs_len = mem_idx;
      uint8_t* lcs = (uint8_t*)malloc(lcs_len);
      // verbose("lcs: " + std::string(read->seq.s).substr(0,lcs_len));
      // Convert A,C,G,T,N into 0,1,2,3,4
      // The left context is reversed
      for (size_t i = 0; i < lcs_len; ++i)
        lcs[lcs_len -i -1] = seq_nt4_table[(int)read->seq.s[i]];

      // rcs: right context sequence
      size_t rcs_occ = (mem_idx + mem_len); // The first character of the right context
      size_t rcs_len = read->seq.l - rcs_occ;
      uint8_t* rcs = (uint8_t*)malloc(rcs_len);
      // verbose("rcs: " + std::string(read->seq.s).substr(rcs_occ,rcs_len));
      // Convert A,C,G,T,N into 0,1,2,3,4
      for (size_t i = 0; i < rcs_len; ++i) 
        rcs[i] = seq_nt4_table[(int)read->seq.s[ rcs_occ + i]];

      int32_t min_score = 20 + 8 * log(read->seq.l);
      // verbose("Number of occurrences: " + std::to_string(occs.size()));
      // For all the occurrences align
      for(auto curr_mem_pos: occs)
      {
        int32_t score = extend(
          curr_mem_pos,
          mem_len,
          lcs,   // Left context of the read
          lcs_len, // Left context of the read lngth
          rcs,   // Right context of the read
          rcs_len // Right context of the read length
        );

        if(score > min_score)
        {
          extend(curr_mem_pos,mem_len,lcs,lcs_len,rcs,rcs_len,false,0,min_score,read,strand,out); 
          aligned = true;
        }

      }
      delete lcs;
      delete rcs;
    }
    return aligned;
  }

  typedef struct mem_t{
    size_t pos = 0; // Position in the reference
    size_t len = 0; // Length
    size_t idx = 0; // Position in the pattern
    std::vector<size_t> occs; // List of occurrences of the MEM

    mem_t(size_t p, size_t l, size_t i)
    {
      pos = p; // Position in the reference
      len = l; // Length of the MEM
      idx = i; // Position in the read
    }

  } mem_t;

  // Fill the vector of occurrences of the mem_t data structure
  // TODO: remove the checks for the lengths once Dominik fix the lce queries
  void find_MEM_occs(mem_t& mem)
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
  bool align_chain(kseq_t *read, FILE* out)
  {
    // std::vector<mem_t> mems;


    bool aligned = false;

    // // Find MEMs.
    // auto pointers = ms.query(read->seq.s, read->seq.l);
    // size_t l = 0;   // Current match length
    // size_t pl = 0;  // Previous match length
    // size_t n_Ns = 0;
    // for (size_t i = 0; i < pointers.size(); ++i)
    // {
    //   size_t pos = pointers[i];
    //   while ((i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
    //   {
    //     if(read->seq.s[i + l] == 'N') n_Ns++;
    //     else n_Ns = 0;
    //     ++l;
    //   }


    //   // Update MEMs
    //   if (l >= pl and n_Ns < l and l >= min_len)
    //   {
    //     mems.push_back(mem_t(pointers[i],l,i));
    //     find_MEM_occs(mems.back());
    //   }

    //   // Compute next match length
    //   pl = l;
    //   l = (l == 0 ? 0 : (l - 1));
    // }

    // if( mems.size() <= 0)
    //   return false;

    // /************* minimap2 dynamic programming for mem chaining ***************/
    // /* https://github.com/lh3/minimap2/blob/master/chain.c */

    // // Sort anchors
    // // Lamda helper to sort the anchors
    // auto cmp = [&] (std::pair<size_t, size_t> i, std::pair<size_t, size_t> j) -> bool {
    //   return (mems[i.first].occs[i.second] + mems[i.first].len - 1) > (mems[j.first].occs[j.second] + mems[j.first].len - 1);
    // };

    // // TODO: improve this initialization
    // size_t tot_mem_length = 0;

    // std::vector< std::pair< size_t, size_t > > anchors;
    // for(size_t i = 0; i < mems.size(); ++i)
    // {
    //   for(size_t j = 0; j < mems[i].occs.size(); ++j)
    //     anchors.push_back(make_pair(i,j));
    //   tot_mem_length +=  mems[i].len * mems[i].occs.size();
    // }

    // float avg_mem_length = (float)tot_mem_length / anchors.size();

    // std::sort(anchors.begin(),anchors.end(),cmp);

    // // Dynamic programming

    // // TODO: Parameters to be defined
    // const ll G = LLONG_MAX;
    // const ll max_dist_x = 1000000;//LLONG_MAX;
    // const ll max_dist_y = 1000000;//LLONG_MAX;
    // const ll max_iter = 50;
    // const ll max_pred = 50;
    // const ll min_chain_score = 1;
    // const ll min_chain_length = 1;


    // std::vector<ll> f(anchors.size(),0); // Score ending in position i
    // std::vector<ll> p(anchors.size(),0); // Position of the next anchor giving the max when chained with the one in position i
    // std::vector<ll> msc(anchors.size(),0); // Max score up to position i
    // std::vector<ll> t(anchors.size(),0); // Stores i in position p[j], if i is chained with j. See heuristics in minimap2

    // ll lb = 0;
    // // For all the anchors
    // for(size_t i = 0 ; i < anchors.size(); ++i)
    // {
    //   // Get anchor i
    //   const auto a_i = anchors[i];
    //   const size_t x_i = mems[a_i.first].occs[a_i.second] + mems[a_i.first].len - 1;
    //   const size_t y_i = mems[a_i.first].idx + mems[a_i.first].len - 1;
    //   const size_t w_i = mems[a_i.first].len;

    //   ll max_f = w_i;
    //   ll max_j = -1;
    //   size_t n_pred = 0;
    //   // For all previous anchors
    //   // Heuristics from minimap2 -> do not try more than 50 anchors
    //   if(i - lb > max_iter) lb = i - max_iter;
    //   for(ll j = i-1; j > lb; --j)
    //   {
    //     const auto a_j = anchors[j];
    //     const size_t x_j = mems[a_j.first].occs[a_j.second] + mems[a_j.first].len - 1;
    //     const size_t y_j = mems[a_j.first].idx + mems[a_j.first].len - 1;

    //     // If the current anchor is too far, exit.
    //     if(x_i > x_j + max_dist_x)
    //     {
    //       j = lb + 1;
    //       continue;
    //     }

    //     const size_t x_d = x_i - x_j;
    //     const size_t y_d = y_i - y_j;
    //     const size_t l = (y_d > x_d ? (y_d - x_d) : (x_d - y_d));
    //     const uint32_t ilog_l = (l > 0? ilog2_32(l): 0);

    //     if(y_j >= y_i or max(y_d, x_d) > G)
    //       continue;

    //     const ll alpha = min(min(y_d,x_d),w_i);
    //     const ll beta = (l > 0? (ll)(.01 * l * avg_mem_length) + ilog_l >> 1 : 0); 

    //     ll score = f[j] + (alpha - beta);

    //     if( score > max_f )
    //     {
    //       max_f = score;
    //       max_j = j;
    //       if( n_pred > 0) --n_pred;
    //     }
    //     else // minimap2: If i is chained wth j, than chaining i with a predecessor of j does not improve the score
    //       if (t[j] == i and (++n_pred > max_pred))            
    //         break;
          
    //     if(p[j] > 0) t[p[j]] = i;
    //   }

    //   f[i] = max_f;
    //   p[i] = max_j;
    //   if( max_j >= 0 and msc[max_j] > max_f)
    //     msc[i] = msc[max_j];
    //   else
    //     msc[i] = max_f;
    // }

    // // Find the end positions of chains
    // memset(t.data(),0,sizeof(size_t) * t.size());
    // for(size_t i = 0; i < anchors.size(); ++i)
    //   if(p[i] >= 0) t[p[i]] = 1;
    
    // size_t n_chains = 0;
    // for(size_t i = 0; i < anchors.size(); ++i)
    //   if(t[i] == 0 and msc[i] > min_chain_score) n_chains ++;
    
    // // TODO: check if we want to report also non aligned reads
    // if(n_chains == 0)
    //   return false;

    // // TODO: replace the vector of pairs with a lambda for the sort
    // std::vector<std::pair<size_t,size_t>> chain_starts(n_chains); // Stores uniqe chains and their index of uniqe chains in anchors
    // size_t k = 0;
    // for(size_t i = 0; i < n_chains; ++i)
    // {
    //   if(t[i] == 0 and  msc[i] > min_chain_score)
    //   {
    //     size_t j = i;
    //     while(j >= 0 and f[j] < msc[j]) j = p[j]; // Find teh peak that maximizes f
    //     if (j < 0) i = 0;
    //     chain_starts[k++] = make_pair(f[j],j);
    //   }
    // }

    // chain_starts.resize(k), chain_starts.shrink_to_fit();
    // n_chains = chain_starts.size();
    // std::sort(chain_starts.begin(), chain_starts.end(), std::greater<std::pair<size_t,size_t>>());
    
    // std::vector<std::pair<size_t, std::vector<size_t>>> chains;

    // // Backtrack
    // memset(t.data(),0,sizeof(size_t) * t.size());
    // for(size_t i = 0; i < n_chains; ++i)
    // {
    //   size_t j = chain_starts[i].second;
    //   std::vector<size_t> chain;
    //   do {
    //     chain.push_back(j); // stores th reverse of the chain
    //     t[j] = 1;
    //     j = p[j];
    //   } while (j >= 0 && t[j] == 0);
    //   if (j < 0) { // l - prev_l is the length of the chain
    //     if (chain.size() >= min_chain_length)
    //       chains.push_back(std::make_pair(chain_starts[i].first, chain));     
    //   } else if (chain_starts[i].first - f[j] >= min_chain_score) { // Two chains share a common prefix
    //     if (chain.size() >= min_chain_length)
    //       chains.push_back(std::make_pair((chain_starts[i].first - f[j]), chain));
    //   }
    // }

    // // Sort the chains by max scores.
    // std::sort(chains.begin(), chains.end(), std::greater<std::pair<size_t,std::vector<size_t>>>());

    // // Clear space
    // p.resize(0), p.shrink_to_fit();
    // t.resize(0), t.shrink_to_fit();
    // f.resize(0), f.shrink_to_fit();
    // msc.resize(0), msc.shrink_to_fit();


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
    std::vector<std::pair<ll, std::vector<size_t>>> chains;

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


    find_mems(&read_rev,mems_rev);

    // for (size_t i = 0; i < mems_rev.size(); ++i)
    // {
    //   std::cout << "MEM[" << i << "]: \n";
    //   std::cout << "    len:" << mems_rev[i].len << "\n";
    //   std::cout << "    pos:" << mems_rev[i].pos << "\n";
    //   std::cout << "    idx:" << mems_rev[i].idx << "\n";
    //   std::cout << "   occs:" << mems_rev[i].occs.size() << "\n";
    // }

    std::vector< std::pair< size_t, size_t > > anchors_rev;
    std::vector<std::pair<ll, std::vector<size_t>>> chains_rev;

    const bool rev_chains = find_chains(mems_rev, anchors_rev, chains_rev);
    
    if ((not fwd_chains) and (not rev_chains))
    {
      std::string dummy = "";
      write_sam(0, 0, 0, 0, "human", read, 0, out, dummy, dummy, 0);
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
      if(chains[i].first > chains_rev[j].first)
      {
        different_scores.insert(chains[i].first);
        if (different_scores.size() < 3)
        {
          // Align the chain
          auto chain = chains[i];
          // Reverse the chain order
          std::reverse(chain.second.begin(), chain.second.end());
          // Compute the score of a chain.
          int32_t score = chain_score(chain, anchors, mems, min_score, read);
          best_scores.push_back(std::make_pair(score, i++));
        }
      }
      else
      {
        different_scores.insert(chains_rev[j].first);
        if (different_scores.size() < 3)
        {
          // Align the chain
          auto chain = chains_rev[j];
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
      different_scores.insert(chains[i].first);
      if (different_scores.size() < 3)
      {
        // Align the chain
        auto chain = chains[i];
        // Reverse the chain order
        std::reverse(chain.second.begin(), chain.second.end());
        // Compute the score of a chain.
        int32_t score = chain_score(chain, anchors, mems, min_score, read);
        best_scores.push_back(std::make_pair(score, i++));
      }
    }
    while (different_scores.size() < 3 and j < chains_rev.size())
    {
      different_scores.insert(chains_rev[j].first);
      if(different_scores.size() < 3)
      {
        // Align the chain
        auto chain = chains_rev[j];
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
      write_sam(0, 0, 0, 0, "human", read, 0, out, dummy, dummy, 0);
      return false;
    }

    int32_t score2 = best_scores[1].first;
    int32_t score = 0;

    if(best_scores[0].second >= off)
    { // Reverse case
      j = best_scores[0].second - off;
      // Align the chain
      auto chain = chains_rev[j];
      // Reverse the chain order
      std::reverse(chain.second.begin(), chain.second.end());
      // Compute the score of a chain.
      score = chain_score(chain, anchors_rev, mems_rev, min_score, &read_rev, false, score2, 1, out);
    }
    else
    { // Forward case
      i = best_scores[0].second;
      // Align the chain
      auto chain = chains[i];
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
      write_sam(0, 0, 0, 0, "human", read, 0, out, dummy, dummy, 0);
    }

    return aligned;
  }

  // Given a set of mems, find the chains
  bool find_chains(
    const   std::vector<mem_t>& mems,
            std::vector< std::pair< size_t, size_t > >& anchors,
            std::vector<std::pair<ll, std::vector<size_t>>>& chains        
    )
  {
    /************* minimap2 dynamic programming for mem chaining ***************/
    /* https://github.com/lh3/minimap2/blob/master/chain.c */

    // Sort anchors
    // Lamda helper to sort the anchors
    auto cmp = [&] (std::pair<size_t, size_t> i, std::pair<size_t, size_t> j) -> bool {
      return (mems[i.first].occs[i.second] + mems[i.first].len - 1) < (mems[j.first].occs[j.second] + mems[j.first].len - 1);
    };

    // TODO: improve this initialization
    size_t tot_mem_length = 0;

    // std::vector< std::pair< size_t, size_t > > anchors;
    for(size_t i = 0; i < mems.size(); ++i)
    {
      for(size_t j = 0; j < mems[i].occs.size(); ++j)
        anchors.push_back(make_pair(i,j));
      tot_mem_length +=  mems[i].len * mems[i].occs.size();
    }

    float avg_mem_length = (float)tot_mem_length / anchors.size();

    std::sort(anchors.begin(),anchors.end(),cmp);

    // Dynamic programming

    // TODO: Parameters to be defined
    const ll G = LLONG_MAX;
    const ll max_dist_x = 100;//LLONG_MAX;
    const ll max_dist_y = 100;//LLONG_MAX;
    const ll max_iter = 50;
    const ll max_pred = 50;
    const ll min_chain_score = 1;
    const ll min_chain_length = 1;


    std::vector<ll> f(anchors.size(),0); // Score ending in position i
    std::vector<ll> p(anchors.size(),0); // Position of the next anchor giving the max when chained with the one in position i
    std::vector<ll> msc(anchors.size(),0); // Max score up to position i
    std::vector<ll> t(anchors.size(),0); // Stores i in position p[j], if i is chained with j. See heuristics in minimap2

    ll lb = 0;
    // For all the anchors
    for(size_t i = 0 ; i < anchors.size(); ++i)
    {
      // Get anchor i
      const auto a_i = anchors[i];
      // const ll k_i = mems[a_i.first].occs[a_i.second];
      const ll x_i = mems[a_i.first].occs[a_i.second] + mems[a_i.first].len - 1;
      // const ll z_i = mems[a_i.first].idx;
      const ll y_i = mems[a_i.first].idx + mems[a_i.first].len - 1;
      const ll w_i = mems[a_i.first].len;

      ll max_f = w_i;
      ll max_j = -1;
      size_t n_pred = 0;
      // For all previous anchors
      // Heuristics from minimap2 -> do not try more than 50 anchors
      if(i - lb > max_iter) lb = i - max_iter;
      for(ll j = i-1; j >= lb; --j)
      {
        const auto a_j = anchors[j];
        const ll x_j = mems[a_j.first].occs[a_j.second] + mems[a_j.first].len - 1;
        const ll y_j = mems[a_j.first].idx + mems[a_j.first].len - 1;

        // If the current anchor is too far, exit.
        if(x_i > x_j + max_dist_x)
        {
          // j = lb - 1;
          lb = j;
          continue;
        }
        // // skip if incompatible
        // if(k_i < x_j or z_i < y_j) continue;

        const ll x_d = x_i - x_j;
        const ll y_d = y_i - y_j;
        const int32_t l = (y_d > x_d ? (y_d - x_d) : (x_d - y_d));
        const uint32_t ilog_l = (l > 0? ilog2_32(l): 0);

        if(y_j >= y_i or max(y_d, x_d) > G)
          continue;

        const ll alpha = min(min(y_d,x_d),w_i);
        const ll beta = (l > 0? (ll)(.01 * l * avg_mem_length) + ilog_l >> 1 : 0); 

        ll score = f[j] + (alpha - beta);

        if( score > max_f )
        {
          max_f = score;
          max_j = j;
          if( n_pred > 0) --n_pred;
        }
        else // minimap2: If i is chained wth j, than chaining i with a predecessor of j does not improve the score
          if (t[j] == i and (++n_pred > max_pred))            
            break;
          
        if(p[j] > 0) t[p[j]] = i;
      }

      f[i] = max_f;
      p[i] = max_j;
      if( max_j >= 0 and msc[max_j] > max_f)
        msc[i] = msc[max_j];
      else
        msc[i] = max_f;
    }

    // Find the end positions of chains
    memset(t.data(),0,sizeof(size_t) * t.size());
    for(size_t i = 0; i < anchors.size(); ++i)
      if(p[i] >= 0) t[p[i]] = 1;
    
    size_t n_chains = 0;
    for(size_t i = 0; i < anchors.size(); ++i)
      if(t[i] == 0 and msc[i] > min_chain_score) n_chains ++;
    
    // TODO: check if we want to report also non aligned reads
    if(n_chains == 0)
      return false;

    // TODO: replace the vector of pairs with a lambda for the sort
    std::vector<std::pair<ll,size_t>> chain_starts(n_chains); // Stores uniqe chains and their index of uniqe chains in anchors
    size_t k = 0;
    for(size_t i = 0; i < anchors.size(); ++i)
    {
      if(t[i] == 0 and  msc[i] > min_chain_score)
      {
        size_t j = i;
        while(j >= 0 and f[j] < msc[j]) j = p[j]; // Find teh peak that maximizes f
        if (j < 0) i = 0;
        chain_starts[k++] = make_pair(f[j],j);
      }
    }

    chain_starts.resize(k), chain_starts.shrink_to_fit();
    n_chains = chain_starts.size();
    std::sort(chain_starts.begin(), chain_starts.end(), std::greater<std::pair<ll,size_t>>());
    
    // std::vector<std::pair<size_t, std::vector<size_t>>> chains;

    // Backtrack
    memset(t.data(),0,sizeof(size_t) * t.size());
    for(size_t i = 0; i < n_chains; ++i)
    {
      ll j = chain_starts[i].second;
      std::vector<size_t> chain;
      do {
        chain.push_back(j); // stores th reverse of the chain
        t[j] = 1;
        j = p[j];
      } while (j >= 0 && t[j] == 0);
      if (j < 0) { // l - prev_l is the length of the chain
        if (chain.size() >= min_chain_length)
          chains.push_back(std::make_pair(chain_starts[i].first, chain));     
      } else if (chain_starts[i].first - f[j] >= min_chain_score) { // Two chains share a common prefix
        if (chain.size() >= min_chain_length)
          chains.push_back(std::make_pair((chain_starts[i].first - f[j]), chain));
      }
    }

    // Sort the chains by max scores.
    std::sort(chains.begin(), chains.end(), std::greater<std::pair<ll,std::vector<size_t>>>());

    // Clear space
    p.resize(0), p.shrink_to_fit();
    t.resize(0), t.shrink_to_fit();
    f.resize(0), f.shrink_to_fit();
    msc.resize(0), msc.shrink_to_fit();

    return true;
  }

  void find_mems(
    const kseq_t *read,
    std::vector<mem_t>& mems
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
        mems.push_back(mem_t(pointers[i],l,i));
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
      FILE *out = nullptr)
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

      if(!score_only and score > min_score)
      {
        fill_chain(mems,chain_anchors,lcs,lcs_len,rcs,rcs_len,read,score_only,0,min_score,strand,out); 
      }

      delete rcs;
      delete lcs;

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
    flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT;

    if(score_only) 
      flag = KSW_EZ_SCORE_ONLY;

    int score_lc = 0;
    int score_rc = 0;

    // TODO: Update end_bonus according to the MEM contribution to the score
    
    // Extract the context from the reference
    // lc: left context
    ksw_reset_extz(&ez_lc);
    if(lcs_len > 0)
    {
      size_t lc_occ = (mem_pos > 100 ? mem_pos - 100 : 0);
      size_t lc_len = (mem_pos > 100 ? 100 : 100 - mem_pos);
      char *tmp_lc = (char *)malloc(100);
      ra.expandSubstr(lc_occ, lc_len, tmp_lc);
      // verbose("lc: " + std::string(lc));
      // Convert A,C,G,T,N into 0,1,2,3,4
      // The left context is reversed
      uint8_t *lc = (uint8_t *)malloc(100);
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
      size_t rc_len = (rc_occ < n-100 ? 100 : n - rc_occ);
      char *rc = (char *)malloc(100);
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
        std::pair<std::string,size_t> md_nm = write_MD_core((uint8_t*)ref,seq,ez.cigar,ez.n_cigar,tmp,0);

        write_sam(ez.score,score2,min_score,ref_pos,"human",read,strand,out,cigar_s,md_nm.first,md_nm.second);
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
        std::pair<std::string,size_t> md_nm = write_MD_core((uint8_t*)ref,seq,cigar,n_cigar,tmp,0);

        write_sam(score,score2,min_score,ref_pos,"human",read,strand,out,cigar_s,md_nm.first,md_nm.second);

        delete cigar;
      }
      delete tmp;
      delete ref;
      delete seq;
    }


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
    const int32_t score2 = 0,    // The score of the second best alignment
    const int32_t min_score = 0, // The minimum score to call an alignment
    int8_t strand = 0,    // 0: forward aligned ; 1: reverse complement aligned
    FILE *out = nullptr,   // The SAM file pointer
    bool realign = false   // Realign globally the read
    // const bool realign = false   // Realign globally the read
  )
  {
    flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT;

    if(score_only) 
      flag = KSW_EZ_SCORE_ONLY;

    int score_lc = 0;
    int score_rc = 0;

    // TODO: Update end_bonus according to the MEM contribution to the score
    
    // Extract the context from the reference
    // lc: left context of the first mem
    ksw_reset_extz(&ez_lc);
    if(lcs_len > 0)
    {
      size_t mem_pos = mems[anchors[0].first].occs[anchors[0].second];

      size_t lc_occ = (mem_pos > 100 ? mem_pos - 100 : 0);
      size_t lc_len = (mem_pos > 100 ? 100 : 100 - mem_pos);
      char *tmp_lc = (char *)malloc(100);
      ra.expandSubstr(lc_occ, lc_len, tmp_lc);
      // verbose("lc: " + std::string(lc));
      // Convert A,C,G,T,N into 0,1,2,3,4
      // The left context is reversed
      uint8_t *lc = (uint8_t *)malloc(100);
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
      size_t rc_len = (rc_occ < n-100 ? 100 : n - rc_occ);
      char *rc = (char *)malloc(100);
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
        
        std::string cigar_s;
        for(size_t i = 0; i < ez.n_cigar; ++i)
          cigar_s += std::to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];

        // Compute the MD:Z field and thenumber of mismatches
        std::pair<std::string,size_t> md_nm = write_MD_core((uint8_t*)ref,seq,ez.cigar,ez.n_cigar,tmp,0);

        write_sam(ez.score,score2,min_score,ref_pos,"human",read,strand,out,cigar_s,md_nm.first,md_nm.second);
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


    return score;
  }

  // Readapted from https://github.com/lh3/minimap2/blob/c9874e2dc50e32bbff4ded01cf5ec0e9be0a53dd/format.c
  // tmp is a string of length max(reference length, query length)
  static std::pair<std::string,size_t> write_MD_core(const uint8_t *tseq, const uint8_t *qseq, const uint32_t *cigar, const size_t n_cigar, char *tmp, int write_tag)
  {
    std::string mdz;
    int i, q_off, t_off, l_MD = 0, NM = 0;
    if (write_tag) mdz += "MD:Z:"; //printf("MD:Z:");
    for (i = q_off = t_off = 0; i < (int)n_cigar; ++i) {
      int j, op = cigar[i]&0xf, len = cigar[i]>>4;
      assert((op >= 0 && op <= 3) || op == 7 || op == 8);
      if (op == 0 || op == 7 || op == 8) { // match
        for (j = 0; j < len; ++j) {
          if (qseq[q_off + j] != tseq[t_off + j]) {
            mdz += std::to_string(l_MD) + "ACGTN"[tseq[t_off + j]];
            // printf("%d%c", l_MD, "ACGTN"[tseq[t_off + j]]);
            l_MD = 0; 
            ++NM;
          } else ++l_MD;
        }
        q_off += len, t_off += len;
      } else if (op == 1) { // insertion to ref
        q_off += len;
        NM += len;
      } else if (op == 2) { // deletion from ref
        for (j = 0, tmp[len] = 0; j < len; ++j)
          tmp[j] = "ACGTN"[tseq[t_off + j]];
        mdz += std::to_string(l_MD) + "^" + std::string(tmp);
        // printf("%d^%s", l_MD, tmp);
        l_MD = 0;
        t_off += len;
        NM += len;
      } else if (op == 3) { // reference skip
        t_off += len;
      }
    }
    if (l_MD > 0) mdz += std::to_string(l_MD);//printf("%d", l_MD);
    // assert(t_off == r->re - r->rs && q_off == r->qe - r->qs);
    return make_pair(mdz,NM);
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


  // Adapted from https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/master/src/main.c
  void write_sam(const int32_t score,
                        const int32_t score2,
                        const int32_t min_score,
                        size_t ref_pos,
                        const char *ref_seq_name,
                        const kseq_t *read,
                        int8_t strand,      // 0: forward aligned ; 1: reverse complement aligned
                        FILE *out,
                        std::string &cigar,
                        std::string &md,
                        size_t mismatches) 
  {
    // Sam format output
    fprintf(out, "%s\t", read->name.s);
    if (score == 0)
      fprintf(out, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
    else
    {
      int32_t c, p;
      // uint32_t mapq = -4.343 * log(1 - (double)abs(score - score2) / (double)score);
      // mapq = (uint32_t)(mapq + 4.99);
      // mapq = mapq < 254 ? mapq : 254;
      uint32_t mapq = compute_mapq(score,score2,min_score,read->seq.l); 
      if (strand)
        fprintf(out, "16\t");
      else
        fprintf(out, "0\t");
      // TODO: Find the correct reference name.
      fprintf(out, "%s\t%d\t%d\t", ref_seq_name, ref_pos + 1, mapq);
      fprintf(out, "%s", cigar.c_str());
      fprintf(out, "\t*\t0\t0\t");
      fprintf(out, "%s", read->seq.s);
      fprintf(out, "\t");
      if (read->qual.s && strand)
      {
        for (p = read->qual.l - 1; p >= 0; --p)
          fprintf(out, "%c", read->qual.s[p]);
      }
      else if (read->qual.s)
        fprintf(out, "%s", read->qual.s);
      else
        fprintf(out, "*");
      fprintf(out, "\tAS:i:%d", score);
      fprintf(out, "\tNM:i:%d", mismatches);
      if (score2 > 0)
        fprintf(out, "\tZS:i:%d", score2);
      fprintf(out, "\tMD:Z:%s\n", md.c_str());
    }
  }

  /*!
    Compute the mapping quality of the alignment
    Inspired from https://github.com/BenLangmead/bowtie2/blob/4512b199768e562e8627ffdfd9253affc96f6fc6/unique.h
  */
  size_t compute_mapq(
    const int32_t score,      // Best alignment score
    const int32_t score2,     // Second best alignemt score
    const int32_t min_score,  // Minimum alignemt score
    const size_t read_l       // Read length
  )
  {
    int32_t max_score = read_l * smatch;
    int32_t best = max_score - score;
    size_t best_bin = (size_t)((double)best * (10.0 / (double)(max_score - min_score)) + 0.5);
    if(score2 > min_score)
    {
      int32_t diff = score - score2;
      size_t diff_bin = (size_t)((double)diff * (10.0 / (double)(max_score - min_score)) + 0.5);
      if(best == max_score)
        return unp_sec_perf[best_bin];
      else
        return unp_sec[diff_bin][best_bin];
    }
    else
    {
      if(best == max_score)
        return unp_nosec_perf;
      else
        return unp_nosec[best_bin];
    }
  }

  static std::string print_BLAST_like(const uint8_t *tseq, const uint8_t *qseq, const uint32_t *cigar, const size_t n_cigar)
  {
    std::string target_o;
    std::string bars_o;
    std::string seq_o;


    int i, q_off, t_off, l_MD = 0;
    for (i = q_off = t_off = 0; i < (int)n_cigar; ++i) {
      int j, op = cigar[i]&0xf, len = cigar[i]>>4;
      assert((op >= 0 && op <= 3) || op == 7 || op == 8);
      if (op == 0 || op == 7 || op == 8) { // match
        for (j = 0; j < len; ++j) {
          if (qseq[q_off + j] != tseq[t_off + j]) {
            bars_o +="*";
          } else {
            bars_o +="|";
          }
          target_o += std::to_string(tseq[t_off + j]);
          seq_o += std::to_string(qseq[q_off + j]);
        }
        q_off += len, t_off += len;
      } else if (op == 1) { // insertion to ref
        for (j = 0; j < len; ++j) {
          target_o += " ";
          bars_o += " ";
          seq_o += std::to_string(qseq[q_off + j]);
        }
        q_off += len;
      } else if (op == 2) { // deletion from ref
        for (j = 0; j < len; ++j) {
          seq_o += " ";
          bars_o += " ";
          target_o += std::to_string(tseq[t_off + j]);
        }
        t_off += len;
      } else if (op == 3) { // reference skip
        for (j = 0; j < len; ++j) {
          seq_o += " ";
          bars_o += " ";
          target_o += std::to_string(tseq[t_off + j]);
        }
        t_off += len;
      }
    }
    return target_o + "\n" + bars_o + "\n" + seq_o + "\n";
  }

protected:
  ms_pointers<> ms;
  slp_t ra;
  // SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd> ra;

  size_t min_len = 0;
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

  int8_t smatch = 2;      // Match score default
  int8_t smismatch = 4;   // Mismatch score default
  int8_t gapo = 4;        // Gap open penalty
  int8_t gapo2 = 13;      // Gap open penalty
  int8_t gape = 2;        // Gap extension penalty
  int8_t gape2 = 1;       // Gap extension penalty
  int end_bonus = 400;    // Bonus to add at the extension score to declare the alignment
  
  int w = -1;             // Band width
  int flag = 0;
  int zdrop = -1;

	void *km = 0;           // Kalloc

  // int8_t max_rseq = 0;

  int m = 5;
  int8_t mat[25];
  // int minsc = 0, xtra = KSW_XSTART;
  // uint8_t *rseq = 0;

  bool forward_only;

  ksw_extz_t ez_lc;
  ksw_extz_t ez_rc;
  ksw_extz_t ez;

  // From https://github.com/BenLangmead/bowtie2/blob/4512b199768e562e8627ffdfd9253affc96f6fc6/unique.cpp
  // There is no valid second-best alignment and the best alignment has a
  // perfect score.
  const uint32_t unp_nosec_perf = 44;

  // There is no valid second-best alignment.  We stratify the alignment
  // score of the best alignment into 10 bins.
  const uint32_t unp_nosec[11] = {
    43, 42, 41, 36, 32, 27, 20, 11, 4, 1, 0
  };

  // The best alignment has a perfect score, and we stratify the distance
  // between best and second-best alignment scores into 10 bins.
  const uint32_t unp_sec_perf[11] = {
    2, 16, 23, 30, 31, 32, 34, 36, 38, 40, 42
  };

  // The best alignment has a non-perfect score, and we stratify both by best
  // alignment score (specifically, the maximum score minus the best "best")
  // and by the distance between the best and second-best alignment scores
  // ("difference").  Each is stratified into 10 bins.  Each row is a
  // difference (smaller elts = smaller differences) and each column is a
  // best score (smaller elts = higher best alignment scores).
  const uint32_t unp_sec[11][11] = {
    {  2,  2,  2,  1,  1, 0, 0, 0, 0, 0, 0},
    { 20, 14,  7,  3,  2, 1, 0, 0, 0, 0, 0},
    { 20, 16, 10,  6,  3, 1, 0, 0, 0, 0, 0},
    { 20, 17, 13,  9,  3, 1, 1, 0, 0, 0, 0},
    { 21, 19, 15,  9,  5, 2, 2, 0, 0, 0, 0},
    { 22, 21, 16, 11, 10, 5, 0, 0, 0, 0, 0},
    { 23, 22, 19, 16, 11, 0, 0, 0, 0, 0, 0},
    { 24, 25, 21, 30,  0, 0, 0, 0, 0, 0, 0},
    { 30, 26, 29,  0,  0, 0, 0, 0, 0, 0, 0},
    { 30, 27,  0,  0,  0, 0, 0, 0, 0, 0, 0},
    { 30,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0},
  };

  //
  // Paired mapping quality:
  //

  // There is no valid second-best alignment and the best alignment has a
  // perfect score.
  const uint32_t pair_nosec_perf = 44;

};


////////////////////////////////////////////////////////////////////////////////
/// Parallel computation
////////////////////////////////////////////////////////////////////////////////

// This should be done using buffering.
size_t next_start_fastq(gzFile fp)
{
  int c;
  // Special case when we arr at the beginning of the file.
  if ((gztell(fp) == 0) && ((c = gzgetc(fp)) != EOF) && c == '@')
    return 0;

  // Strart from the previous character
  gzseek(fp, -1, SEEK_CUR);

  std::vector<std::pair<int, size_t>> window;
  // Find the first new line
  for (size_t i = 0; i < 4; ++i)
  {
    while (((c = gzgetc(fp)) != EOF) && (c != (int)'\n'))
    {
    }
    if (c == EOF)
      return gztell(fp);
    if ((c = gzgetc(fp)) == EOF)
      return gztell(fp);
    window.push_back(std::make_pair(c, gztell(fp) - 1));
  }

  for (size_t i = 0; i < 2; ++i)
  {
    if (window[i].first == '@' && window[i + 2].first == '+')
      return window[i].second;
    if (window[i].first == '+' && window[i + 2].first == '@')
      return window[i + 2].second;
  }

  return gztell(fp);
}

// test if the file is gzipped
static inline bool is_gzipped(std::string filename)
{
  FILE *fp = fopen(filename.c_str(), "rb");
  if(fp == NULL) error("Opening file " + filename);
  int byte1 = 0, byte2 = 0;
  fread(&byte1, sizeof(char), 1, fp);
  fread(&byte2, sizeof(char), 1, fp);
  fclose(fp);
  return (byte1 == 0x1f && byte2 == 0x8b);
}

// Return the length of the file
// Assumes that the file is not compressed
static inline size_t get_file_size(std::string filename)
{
  if (is_gzipped(filename))
  {
    std::cerr << "The input is gzipped!" << std::endl;
    return -1;
  }
  FILE *fp = fopen(filename.c_str(), "r");
  fseek(fp, 0L, SEEK_END);
  size_t size = ftell(fp);
  fclose(fp);
  return size;
}

std::vector<size_t> split_fastq(std::string filename, size_t n_threads)
{
  //Precondition: the file is not gzipped
  // scan file for start positions and execute threads
  size_t size = get_file_size(filename);

  gzFile fp = gzopen(filename.c_str(), "r");
  if (fp == Z_NULL)
  {
    throw new std::runtime_error("Cannot open input file " + filename);
  }

  std::vector<size_t> starts(n_threads + 1);
  for (int i = 0; i < n_threads + 1; ++i)
  {
    size_t start = (size_t)((size * i) / n_threads);
    gzseek(fp, start, SEEK_SET);
    starts[i] = next_start_fastq(fp);
  }
  gzclose(fp);
  return starts;
}

template <typename aligner_t>
struct mt_param_t
{
  // Parameters
  aligner_t *aligner;
  std::string pattern_filename;
  std::string sam_filename;
  size_t start;
  size_t end;
  size_t wk_id;
  // Return values
  size_t n_reads;
  size_t n_aligned_reads;
};

template <typename aligner_t>
void *mt_align_worker(void *param)
{
  mt_param_t<aligner_t> *p = (mt_param_t<aligner_t>*) param;
  size_t n_reads = 0;
  size_t n_aligned_reads = 0;

  FILE *sam_fd;
  gzFile fp;

  if ((sam_fd = fopen(p->sam_filename.c_str(), "w")) == nullptr)
    error("open() file " + p->sam_filename + " failed");

  if ((fp = gzopen(p->pattern_filename.c_str(), "r")) == Z_NULL)
    error("open() file " + p->pattern_filename + " failed");

  gzseek(fp, p->start, SEEK_SET);

  kseq_t rev;
  int l;

  kseq_t *seq = kseq_init(fp);
  while ((ks_tell(seq) < p->end) && ((l = kseq_read(seq)) >= 0))
  {

    bool fwd_align = p->aligner->align(seq, sam_fd, 0);

    //copy seq
    copy_kseq_t(&rev, seq);

    for (size_t i = 0; i < seq->seq.l; ++i)
      rev.seq.s[i] = complement(seq->seq.s[seq->seq.l - i - 1]);

    if (rev.seq.m > rev.seq.l)
      rev.seq.s[rev.seq.l] = 0;

    bool rev_align = p->aligner->align(&rev, sam_fd, 1);

    if (fwd_align or rev_align)
      n_aligned_reads++;
    n_reads++;

    free(rev.name.s);
    free(rev.comment.s);
    free(rev.seq.s);
    free(rev.qual.s);
  }

  verbose("Number of aligned reads block ", p->wk_id, " : ", n_aligned_reads, "/", n_reads);
  p->n_reads = n_reads;
  p->n_aligned_reads = n_aligned_reads;
  kseq_destroy(seq);
  gzclose(fp);
  fclose(sam_fd);

  return NULL;
}
template <typename aligner_t>
size_t mt_align(aligner_t *aligner, std::string pattern_filename, std::string sam_filename, size_t n_threads)
{
  pthread_t t[n_threads] = {0};
  mt_param_t<aligner_t> params[n_threads];
  std::vector<size_t> starts = split_fastq(pattern_filename, n_threads);
  for(size_t i = 0; i < n_threads; ++i)
  {
    params[i].aligner = aligner;
    params[i].pattern_filename = pattern_filename;
    params[i].sam_filename = sam_filename + "_" + std::to_string(i) + ".sam";
    params[i].start = starts[i];
    params[i].end = starts[i+1];
    params[i].wk_id = i;
    xpthread_create(&t[i], NULL, &mt_align_worker<aligner_t>, &params[i], __LINE__, __FILE__);
  }

  size_t tot_reads = 0;
  size_t tot_aligned_reads = 0;

  for(size_t i = 0; i < n_threads; ++i)
  {
    xpthread_join(t[i],NULL,__LINE__,__FILE__);
  }

  sleep(5);
  for(size_t i = 0; i < n_threads; ++i)
  {
    tot_reads += params[i].n_reads;
    tot_aligned_reads += params[i].n_aligned_reads;
  }

  verbose("Number of aligned reads: ", tot_aligned_reads, "/", tot_reads);
  return tot_aligned_reads;
}


////////////////////////////////////////////////////////////////////////////////
/// Single Thread
////////////////////////////////////////////////////////////////////////////////
template<typename aligner_t>
size_t st_align(aligner_t *aligner, std::string pattern_filename, std::string sam_filename)
{
  size_t n_reads = 0;
  size_t n_aligned_reads = 0;
  kseq_t rev;
  int l;
  FILE *sam_fd;

  sam_filename += ".sam";

  if ((sam_fd = fopen(sam_filename.c_str(), "w")) == nullptr)
    error("open() file " + sam_filename + " failed");

  gzFile fp = gzopen(pattern_filename.c_str(), "r");
  kseq_t* seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0)
  {

    bool fwd_align = aligner->align(seq, sam_fd, 0);

    //copy seq
    copy_kseq_t(&rev, seq);

    for (size_t i = 0; i < seq->seq.l; ++i)
      rev.seq.s[i] = complement(seq->seq.s[seq->seq.l - i - 1]);

    if (rev.seq.m > rev.seq.l)
      rev.seq.s[rev.seq.l] = 0;

    bool rev_align = aligner->align(&rev, sam_fd, 1);

    if (fwd_align or rev_align)
      n_aligned_reads++;
    n_reads++;

    free(rev.name.s);
    free(rev.comment.s);
    free(rev.seq.s);
    free(rev.qual.s);
  }

  verbose("Number of aligned reads: ", n_aligned_reads, "/", n_reads);
  kseq_destroy(seq);
  gzclose(fp);
  fclose(sam_fd);

  sleep(5);

  return n_aligned_reads;
}

//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
  std::string filename = "";
  size_t w = 10; // sliding window size and its default
  bool store = false; // store the data structure in the file
  bool memo  = false; // print the memory usage
  bool csv   = false; // print stats on stderr in csv format
  bool rle   = false; // outpt RLBWT
  std::string patterns = ""; // path to patterns file
  size_t l = 25; // minumum MEM length
  size_t th = 1; // number of threads
  bool is_fasta = false; // read a fasta file
  bool shaped_slp = false; // use shaped slp
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string usage("usage: " + std::string(argv[0]) + " infile [-s store] [-m memo] [-c csv] [-p patterns] [-f fasta] [-r rle] [-t threads] [-l len] [-q shaped_slp]\n\n" +
                    "Computes the pfp data structures of infile, provided that infile.parse, infile.dict, and infile.occ exists.\n" +
                    "     wsize: [integer] - sliding window size (def. 10)\n" +
                    "     store: [boolean] - store the data structure in infile.pfp.ds. (def. false)\n" +
                    "      memo: [boolean] - print the data structure memory usage. (def. false)\n" +
                    "     fasta: [boolean] - the input file is a fasta file. (def. false)\n" +
                    "       rle: [boolean] - output run length encoded BWT. (def. false)\n" +
                    "shaped_slp: [boolean] - use shaped slp. (def. false)\n" +
                    "   pattens: [string]  - path to patterns file.\n" +
                    "       len: [integer] - minimum MEM lengt (def. 25)\n" +
                    "    thread: [integer] - number of threads (def. 1)\n" +
                    "       csv: [boolean] - print the stats in csv form on strerr. (def. false)\n");

  std::string sarg;
  while ((c = getopt(argc, argv, "w:smcfql:rhp:t:")) != -1)
  {
    switch (c)
    {
    case 'w':
      sarg.assign(optarg);
      arg.w = stoi(sarg);
      break;
    case 's':
      arg.store = true;
      break;
    case 'm':
      arg.memo = true;
      break;
    case 'c':
      arg.csv = true;
      break;
    case 'r':
      arg.rle = true;
      break;
    case 'p':
      arg.patterns.assign(optarg);
      break;
    case 'l':
      sarg.assign(optarg);
      arg.l = stoi(sarg);
      break;
    case 't':
      sarg.assign(optarg);
      arg.th = stoi(sarg);
      break;
    case 'f':
      arg.is_fasta = true;
      break;
    case 'q':
      arg.shaped_slp = true;
      break;
    case 'h':
      error(usage);
    case '?':
      error("Unknown option.\n", usage);
      exit(1);
    }
  }
  // the only input parameter is the file name
  if (argc == optind + 1)
  {
    arg.filename.assign(argv[optind]);
  }
  else
  {
    error("Invalid number of arguments\n", usage);
  }
}

//********** end argument options ********************

template<typename aligner_t>
void dispatcher(Args &args){

  verbose("Construction of the aligner");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  aligner_t aligner(args.filename, args.l);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();

  std::string base_name = basename(args.filename.data());
  std::string sam_filename = args.patterns + "_" + base_name + "_" + std::to_string(args.l);

  if (is_gzipped(args.patterns))
  {
    verbose("The input is gzipped - forcing single thread alignment.");
    args.th = 1;
  }

  if (args.th == 1)
    st_align<aligner_t>(&aligner, args.patterns, sam_filename);
  else
    mt_align<aligner_t>(&aligner, args.patterns, sam_filename, args.th);

  // TODO: Merge the SAM files.

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

  size_t space = 0;
  if (args.memo)
  {
  }

  if (args.store)
  {
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;
}

int main(int argc, char *const argv[])
{

  Args args;
  parseArgs(argc, argv, args);

  if(args.shaped_slp){
    dispatcher<aligner<shaped_slp_t>>(args);
  }else{
    dispatcher<aligner<plain_slp_t>>(args);
  }
  

  return 0;
}