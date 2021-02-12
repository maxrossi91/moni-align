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

#include <ksw2.h>

#include <omp.h>

#include <libgen.h>

#define _REALIGN

// KSEQ_INIT(gzFile, gzread);

class aligner_t
{
public:
  using SelSd = SelectSdvec<>;
  using DagcSd = DirectAccessibleGammaCode<SelSd>;
  // using SlpT = SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd>;

  aligner_t(std::string filename, 
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

    verbose("Loading random access");
    t_insert_start = std::chrono::high_resolution_clock::now();

    std::string filename_slp = filename + ".slp";

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
  ~aligner_t() 
  {
    if(ez_lc.m_cigar > 0)
      delete ez_lc.cigar;
    if(ez_rc.m_cigar > 0)
      delete ez_rc.cigar;
    if(ez.m_cigar > 0)
      delete ez.cigar;
      // NtD
  }


  bool align(kseq_t *read, FILE* out, uint8_t strand)
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
          extend(curr_mem_pos,mem_len,lcs,lcs_len,rcs,rcs_len,false,0,read,strand,out); 
          aligned = true;
        }

        // verbose("Processing occurrence: " + std::to_string(curr_mem_pos));
        // // Extract the context from the reference
        // // lc: left context
        // size_t lc_occ = (curr_mem_pos > 100 ? curr_mem_pos - 100 : 0);
        // size_t lc_len = (curr_mem_pos > 100 ? 100 : 100 - curr_mem_pos);
        // char *lc = (char *)malloc(100);
        // ra.expandSubstr(lc_occ, lc_len, lc);
        // verbose("lc: " + std::string(lc));
        // // Convert A,C,G,T,N into 0,1,2,3,4
        // // The left context is reversed
        // for (size_t i = 0; i < lc_len; ++i)
        //   lc[lc_len -i -1] = seq_nt4_table[(int)lc[i]];

        // // rc: right context
        // size_t rc_occ = curr_mem_pos + mem_len;
        // size_t rc_len = (rc_occ < n-100 ? 100 : n - rc_occ);
        // char *rc = (char *)malloc(100);
        // ra.expandSubstr(rc_occ, rc_len, rc);
        // verbose("rc: " + std::string(rc));
        // // Convert A,C,G,T,N into 0,1,2,3,4
        // for (size_t i = 0; i < rc_len; ++i)
        //   rc[i] = seq_nt4_table[(int)rc[i]];



        // // TODO: Update end_bonus according to the MEM contribution to the score

        // // Queries: lcs and rcs
        // // Targets: lc  and rc
	      // ksw_extz_t ez_lc;
        // ksw_reset_extz(&ez_lc);
        // // TODO: Decide if we want only the score or also the CIGAR
        // flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT;
        // verbose("aligning lc and lcs");
        // ksw_extz2_sse(km, lcs_len, (uint8_t*)lcs, lc_len, (uint8_t*)lc, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez_lc);
        // verbose("lc score: " + std::to_string(ez_lc.mqe));

	      // ksw_extz_t ez_rc;
        // ksw_reset_extz(&ez_rc);

        // flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT;
        // ksw_extz2_sse(km, rcs_len, (uint8_t*)rcs, rc_len, (uint8_t*)rc, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez_rc);
        // verbose("rc score: " + std::to_string(ez_rc.mqe));
        
        // // Check if both left and right extensions reach the end or the query
        // assert(ez_lc.reach_end and ez_rc.reach_end);
        // // Compute the final score
        // int32_t score = mem_len * smatch + ez_lc.mqe + ez_rc.mqe;
        // int32_t score2 = 0;

        // // Concatenate the CIGAR strings
        // std::string cigar;

        // for(size_t i = 0; i < ez_lc.n_cigar; ++i)
        //   cigar += std::to_string(ez_lc.cigar[ez_lc.n_cigar -i -1] >> 4) + "MID"[ez_lc.cigar[ez_lc.n_cigar -i -1] & 0xf];
        // cigar += std::to_string(mem_len) + "M";
        // for(size_t i = 0; i < ez_rc.n_cigar; ++i)
        //   cigar += std::to_string(ez_rc.cigar[i] >> 4) + "MID"[ez_rc.cigar[i] & 0xf];

        // // TODO: Compute the MD:Z field
        // std::string md = "MD:Z";

        // // TODO: Compute number of mismatches
        // size_t mismatch = 0;
        // // Compute starting position in reference
        // size_t ref_pos = curr_mem_pos - ez_lc.mqe_t;

        // if(score >= min_score)
        // {
        //   write_sam(score,score2,ref_pos,"human",read,strand,out,cigar,md,mismatch);
        //   aligned = true;
        // }


        // delete lc;
        // delete rc;
      }
      delete lcs;
      delete rcs;
    }
    return aligned;
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
    const uint8_t* lcs,   // Left context of the read
    const size_t lcs_len, // Left context of the read lngth
    const uint8_t* rcs,   // Right context of the read
    const size_t rcs_len, // Right context of the read length
    const bool score_only = true, // Report only the score
    const int32_t score2 = 0,    // The score of the second best alignment
    const kseq_t *read = nullptr, // The read that has been aligned
    int8_t strand = 0,    // 0: forward aligned ; 1: reverse complement aligned
    FILE *out = nullptr,   // The SAM file pointer
    const bool realign = false   // Realign globally the read
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
      ksw_reset_extz(&ez_lc);
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
      ksw_reset_extz(&ez_rc);
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
        std::pair<std::string,size_t> md_nm = write_MD_core((uint8_t*)ref,seq,ez.cigar,ez.n_cigar,tmp,1);

        write_sam(ez.score,score2,ref_pos,"human",read,strand,out,cigar_s,md_nm.first,md_nm.second);
      }
      else
      {
        // Concatenate the CIGAR strings
        size_t n_cigar = ez_lc.n_cigar + ez_rc.n_cigar + 1;
        uint32_t *cigar = (uint32_t*)calloc(n_cigar,sizeof(uint32_t));
        size_t i = 0;

        for(size_t j = 0; j < ez_lc.n_cigar; ++j)
          cigar[i++] = ez_lc.cigar[ez_lc.n_cigar -j -1];


        if(ez_lc.n_cigar > 0 and (cigar[i-1]& 0xf == 0))
        { // If the previous operation is also an M then merge the two operations
          cigar[i-1] += (((uint32_t)mem_len) << 4);
          --n_cigar;
        }
        else
          cigar[i++] = (((uint32_t)mem_len) << 4);


        if(ez_rc.n_cigar > 0)
        {
          if(ez_rc.cigar[0]& 0xf == 0)
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

        std::string cigar_s;
        for(size_t i = 0; i < n_cigar; ++i)
          cigar_s += std::to_string(cigar[i] >> 4) + "MID"[cigar[i] & 0xf];


        // Compute the MD:Z field and thenumber of mismatches
        std::pair<std::string,size_t> md_nm = write_MD_core((uint8_t*)ref,seq,cigar,n_cigar,tmp,1);

        write_sam(score,score2,ref_pos,"human",read,strand,out,cigar_s,md_nm.first,md_nm.second);

        delete cigar;
      }
      delete tmp;
      delete ref;
      delete seq;
    }


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
        mdz += std::to_string(l_MD) + std::string(tmp);
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
  static void write_sam(const int32_t score,
                        const int32_t score2,
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
      uint32_t mapq = -4.343 * log(1 - (double)abs(score - score2) / (double)score);
      mapq = (uint32_t)(mapq + 4.99);
      mapq = mapq < 254 ? mapq : 254;
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
  SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd> ra;

  size_t min_len = 0;
  size_t aligned_reads = 0;
  size_t n = 0;

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
};

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

typedef struct{
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
} mt_param;

void *mt_align_worker(void *param)
{
  mt_param *p = (mt_param*) param;
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

size_t mt_align(aligner_t *aligner, std::string pattern_filename, std::string sam_filename, size_t n_threads)
{
  pthread_t t[n_threads] = {0};
  mt_param params[n_threads];
  std::vector<size_t> starts = split_fastq(pattern_filename, n_threads);
  for(size_t i = 0; i < n_threads; ++i)
  {
    params[i].aligner = aligner;
    params[i].pattern_filename = pattern_filename;
    params[i].sam_filename = sam_filename + "_" + std::to_string(i) + ".sam";
    params[i].start = starts[i];
    params[i].end = starts[i+1];
    params[i].wk_id = i;
    xpthread_create(&t[i], NULL, &mt_align_worker, &params[i], __LINE__, __FILE__);
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

int main(int argc, char *const argv[])
{

  Args args;
  parseArgs(argc, argv, args);

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

  if(args.th == 1)
    st_align(&aligner,args.patterns,sam_filename);
  else
    mt_align(&aligner,args.patterns,sam_filename,args.th);

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

  return 0;
}