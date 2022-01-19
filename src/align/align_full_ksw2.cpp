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
#define MTIME


#include <sdsl/io.hpp>

#define _REALIGN
#include <aligner_ksw2.hpp>
#include <align_reads_dispatcher.hpp>

#include <common.hpp>
#include <malloc_count.h>

// #include <SelfShapedSlp.hpp>
// #include <DirectAccessibleGammaCode.hpp>
// #include <SelectType.hpp>
// #include <PlainSlp.hpp>
// #include <FixedBitLenCode.hpp>

// #include <ksw2.h>

// #include <omp.h>

// #include <libgen.h>
// #include <kpbseq.h>

// #include <seqidx.hpp>

// MTIME_INIT(3);

// // KSEQ_INIT(gzFile, gzread);

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
// // ////////////////////////////////////////////////////////////////////////////////
// // /// kseq extra
// // ////////////////////////////////////////////////////////////////////////////////

// // static inline size_t ks_tell(kseq_t *seq)
// // {
// //   return gztell(seq->f->f) - seq->f->end + seq->f->begin;
// // }

// // void copy_kstring_t(kstring_t &l, kstring_t &r)
// // {
// //   l.l = r.l;
// //   l.m = r.m;
// //   l.s = (char *)malloc(l.m);
// //   for (size_t i = 0; i < r.m; ++i)
// //     l.s[i] = r.s[i];
// // }
// // void copy_kseq_t(kseq_t *l, kseq_t *r)
// // {
// //   copy_kstring_t(l->name, r->name);
// //   copy_kstring_t(l->comment, r->comment);
// //   copy_kstring_t(l->seq, r->seq);
// //   copy_kstring_t(l->qual, r->qual);
// //   l->last_char = r->last_char;
// // }
// // ////////////////////////////////////////////////////////////////////////////////

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

// template <typename slp_t,
//           typename ms_t>
// class aligner
// {
// public:
//   // using SelSd = SelectSdvec<>;
//   // using DagcSd = DirectAccessibleGammaCode<SelSd>;
//   // using SlpT = SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd>;
//   using ll = long long int;

//   aligner(std::string filename, 
//             size_t min_len_ = 50, 
//             bool forward_only_ = true): 
//                 min_len(min_len_), 
//                 forward_only(forward_only_)
//   {
//     verbose("Loading the matching statistics index");
//     std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

//     std::string filename_ms = filename + ms.get_file_extension();

//     ifstream fs_ms(filename_ms);
//     ms.load(fs_ms);
//     fs_ms.close();

//     std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

//     verbose("Matching statistics index loading complete");
//     verbose("Memory peak: ", malloc_count_peak());
//     verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

//     std::string filename_slp = filename + get_slp_file_extension<slp_t>();
//     verbose("Loading random access file: " + filename_slp);
//     t_insert_start = std::chrono::high_resolution_clock::now();


//     ifstream fs(filename_slp);
//     ra.load(fs);
//     fs.close();

//     n = ra.getLen();

//     t_insert_end = std::chrono::high_resolution_clock::now();

//     verbose("Random access loading complete");
//     verbose("Memory peak: ", malloc_count_peak());
//     verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

//     std::string filename_idx = filename + idx.get_file_extension();
//     verbose("Loading fasta index file: " + filename_idx);
//     t_insert_start = std::chrono::high_resolution_clock::now();


//     ifstream fs_idx(filename_idx);
//     idx.load(fs_idx);
//     fs_idx.close();

//     t_insert_end = std::chrono::high_resolution_clock::now();

//     verbose("Fasta index loading complete");
//     verbose("Memory peak: ", malloc_count_peak());
//     verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

//     verbose("Initialize the local aligner");
//     t_insert_start = std::chrono::high_resolution_clock::now();

//     ksw_gen_simple_mat(m,mat,smatch,-smismatch);

//     t_insert_end = std::chrono::high_resolution_clock::now();

//     verbose("Local aligner initialization complete");
//     verbose("Memory peak: ", malloc_count_peak());
//     verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

//     verbose("Minimum MEM length: ", min_len);

//     // memset(&ez_lc, 0, sizeof(ksw_extz_t));
//     // memset(&ez_rc, 0, sizeof(ksw_extz_t));
//     // memset(&ez, 0, sizeof(ksw_extz_t));
//   }

//   // Destructor
//   ~aligner() 
//   {
//       // NtD
//   }

//   bool align(kseq_t *read, FILE* out, uint8_t strand, bool mem_chaining = true)
//   {
//     if(mem_chaining) // TODO: fix the no chaining alignment
//       return (strand == 0? align_chain(read,out):false);
//     else
//       return align_no_chain(read,out,strand);
//   }

//   // bool align(kseq_t *read, FILE* out, uint8_t strand)
//   bool align_no_chain(kseq_t *read, FILE* out, uint8_t strand)
//   {
//     size_t mem_pos = 0;
//     size_t mem_len = 0;
//     size_t mem_idx = 0;

//     bool aligned = false;

//     // Find MEMs.
//     auto pointers = ms.query(read->seq.s, read->seq.l);
//     std::vector<size_t> lengths(pointers.size());
//     size_t l = 0;
//     size_t n_Ns = 0;
//     for (size_t i = 0; i < pointers.size(); ++i)
//     {
//       size_t pos = pointers[i];
//       while ((i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
//       {
//         if(read->seq.s[i + l] == 'N') n_Ns++;
//         else n_Ns = 0;
//         ++l;
//       }

//       lengths[i] = l;
//       l = (l == 0 ? 0 : (l - 1));

//       // Update MEM
//       if (lengths[i] > mem_len and n_Ns < lengths[i])
//       {
//         mem_len = lengths[i];
//         mem_pos = pointers[i];
//         mem_idx = i;
//       }
//     }

//     // Align the read
//     if (mem_len >= min_len)
//     {
//       // Extract all the occurrences of the MEM
//       std::vector<size_t> occs;
//       occs.push_back(mem_pos);

//       // Phi direction
//       size_t curr = mem_pos;
//       size_t next = ms.Phi(curr);
//       size_t lcp =  lceToRBounded(ra,curr,next,mem_len);
//       while(lcp >= mem_len)
//       {
//         occs.push_back(next);
        
//         curr = next;
//         next = ms.Phi(curr);
//         lcp  = lceToRBounded(ra,curr,next,mem_len);
//         // verbose("Phi: " + std::to_string(lcp));
//         // if(occs.size() > 100)
//         //   error("More than 100 occs Phi" + std::string(read->seq.s));
//       }

//       // Phi_inv direction
//       curr = mem_pos;
//       next = ms.Phi_inv(curr);
//       lcp =  lceToRBounded(ra,curr,next,mem_len);
//       while(lcp >= mem_len)
//       {
//         occs.push_back(next);
        
//         curr = next;
//         next = ms.Phi_inv(curr);
//         lcp  = lceToRBounded(ra,curr,next,mem_len);
//         // verbose("Phi_inv: " + std::to_string(next));
//         // if(occs.size() > 100)
//         //   error("More than 100 occs Phi_inv" + std::string(read->seq.s));
//       }

//       // Extractin left and right context of the read
//       // lcs: left context sequence
//       size_t lcs_len = mem_idx;
//       uint8_t* lcs = (uint8_t*)malloc(lcs_len);
//       // verbose("lcs: " + std::string(read->seq.s).substr(0,lcs_len));
//       // Convert A,C,G,T,N into 0,1,2,3,4
//       // The left context is reversed
//       for (size_t i = 0; i < lcs_len; ++i)
//         lcs[lcs_len -i -1] = seq_nt4_table[(int)read->seq.s[i]];

//       // rcs: right context sequence
//       size_t rcs_occ = (mem_idx + mem_len); // The first character of the right context
//       size_t rcs_len = read->seq.l - rcs_occ;
//       uint8_t* rcs = (uint8_t*)malloc(rcs_len);
//       // verbose("rcs: " + std::string(read->seq.s).substr(rcs_occ,rcs_len));
//       // Convert A,C,G,T,N into 0,1,2,3,4
//       for (size_t i = 0; i < rcs_len; ++i) 
//         rcs[i] = seq_nt4_table[(int)read->seq.s[ rcs_occ + i]];

//       int32_t min_score = 20 + 8 * log(read->seq.l);
//       // verbose("Number of occurrences: " + std::to_string(occs.size()));
//       // For all the occurrences align
//       for(auto curr_mem_pos: occs)
//       {
//         int32_t score = extend(
//           curr_mem_pos,
//           mem_len,
//           lcs,   // Left context of the read
//           lcs_len, // Left context of the read lngth
//           rcs,   // Right context of the read
//           rcs_len // Right context of the read length
//         );

//         if(score >= min_score)
//         {
//           extend(curr_mem_pos,mem_len,lcs,lcs_len,rcs,rcs_len,false,0,min_score,read,strand,out); 
//           aligned = true;
//         }

//       }
//       delete lcs;
//       delete rcs;
//     }
//     return aligned;
//   }

//   #define MATE_1 0 // The MEM cames from mate_1
//   #define MATE_2 1 // The MEM cames from mate_2
//   #define MATE_F 0 // The MEM cames from the forward strand
//   #define MATE_RC 2 // The MEM cames from the reverse-complement strand

//   typedef struct mem_t{
//     size_t pos = 0;  // Position in the reference
//     size_t len = 0;  // Length
//     size_t idx = 0;  // Position in the pattern
//     size_t mate = 0; // Left mate (0) or Right mate (1)
//     size_t rpos = 0; // Position in the read for chaining
//                      // With a Forward-Reverse library
//                      // If the mem is in the FWD strand it is the position of the last character in the read
//                      // If the mem is in the REV strand it is the position of the first character in the read
//     std::vector<size_t> occs; // List of occurrences of the MEM

//     mem_t(size_t p, size_t l, size_t i)
//     {
//       pos = p;  // Position in the reference
//       len = l;  // Length of the MEM
//       idx = i;  // Position in the read
//     }
//     mem_t(size_t p, size_t l, size_t i, size_t m, size_t r)
//     {
//       pos = p;  // Position in the reference
//       len = l;  // Length of the MEM
//       idx = i;  // Position in the read
//       mate = m; // Left mate (0) or Right mate (1)
//       rpos = r; // Position in the read for chaining
//     }

//   } mem_t;

//   typedef struct sam_t{
//     bool reverse = false; // The read is the reverse complement of the original 
//     const kseq_t *read = nullptr; // The read of the SAM entry
//                           // Contains: QNAME, SEQ, and QUAL
//     size_t flag   = 4;    // FLAG: bitwise FLAG
//     size_t pos    = 0;    // POS: 1-based leftmost mapping POSition
//     size_t mapq   = 255;  // MAPQ: MAPping Quality
//     size_t pnext  = 0;    // PNEXT: Position of the mate/next read
//     ll tlen   = 0;        // TLEN: Position of the mate/next read

//     std::string rname = "*"; // RNAME: Reference sequence NAME
//     std::string cigar = "*"; // CIGAR: CIGAR string
//     std::string rnext = "*"; // RNEXT: Reference name of the mate/next read

//     size_t as = 0; // AS: Alignment score generated by aligner
//     size_t nm = 0; // NM: Edit distance to the reference
//     size_t zs = 0; // ZS: Second best score

//     std::string md = ""; // MD: String encoding mismatched and deleted reference bases

//     size_t rlen = 0; // Length of the match in the referenc. Requiredd to compute TLEN

//     // Quind the sam_t struct for a read
//     sam_t(const kseq_t* read_)
//     {
//       read = read_;
//     }

//     // void clear()
//     // {
//     //   reverse = false;   // The read is the reverse complement of the original
//     //   read = nullptr; // The read of the SAM entry
//     //                           // Contains: QNAME, SEQ, and QUAL
//     //   flag = 4;        // FLAG: bitwise FLAG
//     //   pos = 0;         // POS: 1-based leftmost mapping POSition
//     //   mapq = 255;      // MAPQ: MAPping Quality
//     //   pnext = 0;       // PNEXT: Position of the mate/next read
//     //   tlen = 0;        // TLEN: Position of the mate/next read

//     //   rname = "*"; // RNAME: Reference sequence NAME
//     //   cigar = "*"; // CIGAR: CIGAR string
//     //   rnext = "*"; // RNEXT: Reference name of the mate/next read

//     //   as = 0; // AS: Alignment score generated by aligner
//     //   nm = 0; // NM: Edit distance to the reference
//     //   zs = 0; // ZS: Second best score

//     //   md = ""; // MD: String encoding mismatched and deleted reference bases
//     // }

//   } sam_t;


//   inline void write_sam(FILE *out, const sam_t s)
//   {
//     fprintf(out, "%s\t", s.read->name.s);         // QNAME
//     fprintf(out, "%d\t", s.flag);                 // FLAG
//     fprintf(out, "%s\t", s.rname.c_str());        // RNAME
//     fprintf(out, "%d\t", s.pos);                  // POS
//     fprintf(out, "%d\t", s.mapq);                 // MAPQ
//     fprintf(out, "%s\t", s.cigar.c_str());        // CIGAR
//     fprintf(out, "%s\t", s.rnext.c_str());        // RNEXT
//     fprintf(out, "%d\t", s.pnext);                // PNEXT
//     fprintf(out, "%d\t", s.tlen);                 // TLEN
//     fprintf(out, "%s\t", s.read->seq.s);          // SEQ

//     if (s.read->qual.s)                           // QUAL
//       fprintf(out, "%s", s.read->qual.s);
//     else
//       fprintf(out, "*");

//     // Optional TAGs
//     if (!(s.flag & SAM_UNMAPPED))
//     {
//       fprintf(out, "\tAS:i:%d", s.as);            // AS
//       fprintf(out, "\tNM:i:%d", s.nm);            // NM
//       if (s.zs > 0)
//         fprintf(out, "\tZS:i:%d", s.zs);          // ZS
//       fprintf(out, "\tMD:Z:%s\n", s.md.c_str());  // MD
//     }else
//       fprintf(out, "\n");
//   }


//   // Fill the vector of occurrences of the mem_t data structure
//   // TODO: remove the checks for the lengths once Dominik fix the lce queries
//   void find_MEM_occs(mem_t &mem)
//   {
//     mem.occs.push_back(mem.pos);

//     const auto sa_first = ms.get_first_run_sample();
//     const auto sa_last = ms.get_last_run_sample();

//     // Phi direction
//     size_t curr = mem.pos;
//     if(curr != sa_first)
//     {
//       size_t next = ms.Phi(curr);
//       if((n-curr) >= mem.len and (n-next) >= mem.len)
//       {
//         size_t lcp =  lceToRBounded(ra,curr,next,mem.len);
//         while(lcp >= mem.len)
//         {
//           mem.occs.push_back(next);
          
//           curr = next;
//           if(curr != sa_first)
//           {
//             next = ms.Phi(curr);
//             if((n-curr) >= mem.len and (n-next) >= mem.len)
//               lcp  = lceToRBounded(ra,curr,next,mem.len);
//             else
//               lcp = 0;
//           }else lcp = 0;
//         }
//       }
//     }

//     // Phi_inv direction
//     curr = mem.pos;
//     if(curr != sa_last)
//     {
//       size_t next = ms.Phi_inv(curr);
//       if((n-curr) >= mem.len and (n-next) >= mem.len)
//       {
//         size_t lcp =  lceToRBounded(ra,curr,next,mem.len);
//         while(lcp >= mem.len)
//         {
//           mem.occs.push_back(next);
          
//           curr = next;
//           if(curr != sa_last)
//           {
//             next = ms.Phi_inv(curr);
//             if((n-curr) >= mem.len and (n-next) >= mem.len)
//               lcp  = lceToRBounded(ra,curr,next,mem.len);
//             else
//               lcp = 0;
//           }else lcp = 0;
//         }
//       }
//     }

//   }


//   // bool align_chains(kseq_t *read, FILE* out, uint8_t strand)
//   bool align_chain(kseq_t *read, FILE* out)
//   {
//     // std::vector<mem_t> mems;


//     bool aligned = false;

//     // Find MEMs in the forward direction
//     std::vector<mem_t> mems;

//     find_mems(read,mems);

//     // for(size_t i = 0; i < mems.size(); ++i){
//     //   std::cout << "MEM[" << i <<"]: \n";
//     //   std::cout << "    len:" << mems[i].len <<"\n";
//     //   std::cout << "    pos:" << mems[i].pos <<"\n";
//     //   std::cout << "    idx:" << mems[i].idx <<"\n";
//     //   std::cout << "   occs:" << mems[i].occs.size() <<"\n";
//     // }

//     std::vector< std::pair< size_t, size_t > > anchors;
//     std::vector<chain_t> chains;

//     // if(!find_chains(mems,anchors,chains))
//     //   return false;
//     const bool fwd_chains = find_chains(mems, anchors, chains);

//     // Find MEMs in reverse direction
//     std::vector<mem_t> mems_rev;

//     //copy seq
//     kseq_t read_rev;
//     copy_kseq_t(&read_rev, read);

//     for (size_t i = 0; i < read->seq.l; ++i)
//       read_rev.seq.s[i] = complement(read->seq.s[read->seq.l - i - 1]);

//     if (read_rev.seq.m > read_rev.seq.l)
//       read_rev.seq.s[read_rev.seq.l] = 0;


//     find_mems(&read_rev,mems_rev,1);

//     // for (size_t i = 0; i < mems_rev.size(); ++i)
//     // {
//     //   std::cout << "MEM[" << i << "]: \n";
//     //   std::cout << "    len:" << mems_rev[i].len << "\n";
//     //   std::cout << "    pos:" << mems_rev[i].pos << "\n";
//     //   std::cout << "    idx:" << mems_rev[i].idx << "\n";
//     //   std::cout << "   occs:" << mems_rev[i].occs.size() << "\n";
//     // }

//     std::vector< std::pair< size_t, size_t > > anchors_rev;
//     std::vector<chain_t> chains_rev;

//     const bool rev_chains = find_chains(mems_rev, anchors_rev, chains_rev);
    
//     if ((not fwd_chains) and (not rev_chains))
//     {
//       std::string dummy = "";
//       write_sam(0, 0, 0, 0, "*", read, 0, out, dummy, dummy, 0, "*", 0, 0);
//       return false;
//     } 

//     int32_t min_score = 20 + 8 * log(read->seq.l);

//     // Compute the second best score
//     std::vector<std::pair<int32_t, size_t>> best_scores;
//     // get the occurrences of the top 4 best scores
//     std::set<size_t> different_scores;
//     size_t i = 0;
//     size_t j = 0;
//     size_t off = chains.size();
//     while (i < chains.size() and j < chains_rev.size() and different_scores.size() < 3)
//     {
//       if(chains[i].score > chains_rev[j].score)
//       {
//         different_scores.insert(chains[i].score);
//         if (different_scores.size() < 3)
//         {
//           // Align the chain
//           auto chain = std::make_pair(chains[i].score, chains[i].anchors);
//           // Reverse the chain order
//           std::reverse(chain.second.begin(), chain.second.end());
//           // Compute the score of a chain.
//           int32_t score = chain_score(chain, anchors, mems, min_score, read);
//           best_scores.push_back(std::make_pair(score, i++));
//         }
//       }
//       else
//       {
//         different_scores.insert(chains_rev[j].score);
//         if (different_scores.size() < 3)
//         {
//           // Align the chain
//           auto chain = std::make_pair(chains_rev[j].score, chains_rev[j].anchors);
//           // Reverse the chain order
//           std::reverse(chain.second.begin(), chain.second.end());
//           // Compute the score of a chain.
//           int32_t score = chain_score(chain, anchors_rev, mems_rev, min_score, &read_rev);
//           best_scores.push_back(std::make_pair(score, off + j++));
//         }
//       }
//     }
//     while (different_scores.size() < 3 and i < chains.size())
//     {
//       different_scores.insert(chains[i].score);
//       if (different_scores.size() < 3)
//       {
//         // Align the chain
//         auto chain = std::make_pair(chains[i].score, chains[i].anchors);
//         // Reverse the chain order
//         std::reverse(chain.second.begin(), chain.second.end());
//         // Compute the score of a chain.
//         int32_t score = chain_score(chain, anchors, mems, min_score, read);
//         best_scores.push_back(std::make_pair(score, i++));
//       }
//     }
//     while (different_scores.size() < 3 and j < chains_rev.size())
//     {
//       different_scores.insert(chains_rev[j].score);
//       if(different_scores.size() < 3)
//       {
//         // Align the chain
//         auto chain = std::make_pair(chains_rev[j].score, chains_rev[j].anchors);
//         // Reverse the chain order
//         std::reverse(chain.second.begin(), chain.second.end());
//         // Compute the score of a chain.
//         int32_t score = chain_score(chain, anchors_rev, mems_rev, min_score, &read_rev);
//         best_scores.push_back(std::make_pair(score, off + j++));
//       }
//     }

//     if(best_scores.size() < 2)
//       best_scores.push_back(std::make_pair(0,chains.size() + chains_rev.size()));

//     std::sort(best_scores.begin(), best_scores.end(), std::greater<std::pair<int32_t, size_t>>());

//     assert(best_scores.size() > 1);

//     if(best_scores[0].first < min_score)
//     {
//       std::string dummy = "";
//       write_sam(0, 0, 0, 0, "*", read, 0, out, dummy, dummy, 0, "*", 0, 0);
//       return false;
//     }

//     int32_t score2 = best_scores[1].first;
//     int32_t score = 0;

//     if(best_scores[0].second >= off)
//     { // Reverse case
//       j = best_scores[0].second - off;
//       // Align the chain
//       auto chain = std::make_pair(chains_rev[j].score, chains_rev[j].anchors);
//       // Reverse the chain order
//       std::reverse(chain.second.begin(), chain.second.end());
//       // Compute the score of a chain.
//       score = chain_score(chain, anchors_rev, mems_rev, min_score, &read_rev, false, score2, 1, out);
//     }
//     else
//     { // Forward case
//       i = best_scores[0].second;
//       // Align the chain
//       auto chain = std::make_pair(chains[i].score, chains[i].anchors);
//       // Reverse the chain order
//       std::reverse(chain.second.begin(), chain.second.end());
//       // Compute the score of a chain.
//       score = chain_score(chain, anchors, mems, min_score, read, false, score2, 0, out);
//     }

//     if (score >= min_score) aligned = true;
//     // TODO: Implement the topk retrival

//     // std::vector<std::pair<size_t, size_t>> top4;
//     // for(size_t i = 0; i < min(chains.size(),(size_t)2); ++i)
//     //   top4.push_back(std::make_pair(chains[i].first,i));
//     // for(size_t i = 0; i < min(chains_rev.size(),(size_t)2); ++i)
//     //   top4.push_back(std::make_pair(chains_rev[i].first,2 + i));

//     // std::sort(top4.begin(), top4.end(),std::greater<std::pair<size_t,size_t>>());

//     // int32_t score2 = 0;
//     // int32_t score = 0;

//     // //TODO: Rewrite this part that is so ugly
//     // if(top4.size() > 1)
//     // {
//     //   if(top4[1].second > 1)
//     //   {
//     //     auto chain = chains_rev[top4[1].second-2];
//     //     // Reverse the chain order
//     //     std::reverse(chain.second.begin(), chain.second.end());
//     //     // Compute the score of a chain.
//     //     score2 = chain_score(chain,anchors_rev, mems_rev, min_score,(not (top_k > 1)),-1,&read_rev,16,out);
//     //   }
//     //   else
//     //   {
//     //     auto chain = chains[top4[1].second];
//     //     // Reverse the chain order
//     //     std::reverse(chain.second.begin(), chain.second.end());
//     //     // Compute the score of a chain.
//     //     score2 = chain_score(chain,anchors, mems, min_score,(not (top_k > 1)),-1,read,0,out);
//     //   }
//     // }

//     // // Report the high-score chain
//     // if(top4[0].second > 1)
//     // {
//     //   auto chain = chains_rev[top4[0].second-2];
//     //   // Reverse the chain order
//     //   std::reverse(chain.second.begin(), chain.second.end());
//     //   // Compute the score of a chain.
//     //   score = chain_score(chain,anchors_rev, mems_rev, min_score,false,score2,&read_rev,16,out);
//     // }else{
//     //   auto chain = chains[top4[0].second];
//     //   // Reverse the chain order
//     //   std::reverse(chain.second.begin(), chain.second.end());
//     //   // Compute the score of a chain.
//     //   score = chain_score(chain,anchors, mems, min_score,false,score2,read,0,out);
      
//     // }
//     // if(score > min_score) aligned = true;
    



//     // // Compute the alignment score for the top_k chains
//     // // NOTE: for seconadry alignments the MAPQ score is set to 255 (as in Bowtie2)
//     // // QUESTION: Should I compute the scores for all the chains?
//     // size_t j = 0;
//     // size_t j_rev = 0;
//     // for(size_t i = 0; i < min(chains.size() + chains_rev.size(),top_k); ++i)
//     // {
//     //   // QUESTION: Should I compute the score for the top_k distinct alignments?

//     //   if(j < chains.size() and j_rev < chains_rev.size())
//     //   {
//     //     if(chains[j].first > chains_rev[j_rev].first)
//     //     {
//     //       if(i > 1){
//     //         auto chain = chains[j];
//     //         // Reverse the chain order
//     //         std::reverse(chain.second.begin(), chain.second.end());
//     //         // Compute the score of a chain.
//     //         chain_score(chain,anchors, mems, min_score,false,-1,read,0,out);
//     //       }
//     //       j++;
//     //     }
//     //     else
//     //     {
//     //       if(i > 1){
//     //         auto chain = chains_rev[j_rev];
//     //         // Reverse the chain order
//     //         std::reverse(chain.second.begin(), chain.second.end());
//     //         // Compute the score of a chain.
//     //         score2 = chain_score(chain,anchors_rev, mems_rev, min_score,false,-1,&read_rev,16,out);
//     //       }
//     //       j_rev++;
//     //     }
//     //   }
//     //   else if (j < chains.size())
//     //   {
//     //     if(i > 1)
//     //     {
//     //         auto chain = chains[j];
//     //         // Reverse the chain order
//     //         std::reverse(chain.second.begin(), chain.second.end());
//     //         // Compute the score of a chain.
//     //         chain_score(chain,anchors, mems, min_score,false,-1,read,0,out);
//     //     }
//     //     j++;
//     //   }
//     //   else if (j_rev < chains_rev.size())
//     //   {
//     //     if(i > 1)
//     //     {
//     //       auto chain = chains_rev[j_rev];
//     //       // Reverse the chain order
//     //       std::reverse(chain.second.begin(), chain.second.end());
//     //       // Compute the score of a chain.
//     //       chain_score(chain,anchors_rev, mems_rev, min_score,false,-1,&read_rev,16,out);
          
//     //     }
//     //     j_rev++;
//     //   }

//     //   // // Extract the anchors
//     //   // std::vector<std::pair<size_t,size_t>> chain_anchors(chain.second.size());
//     //   // for(size_t i = 0; i < chain_anchors.size(); ++i)
//     //   //   chain_anchors[i] = anchors[chain.second[i]];
//     //   // // Extracting left and right context of the read
//     //   // // lcs: left context sequence
//     //   // size_t lcs_len = mems[chain_anchors[0].first].idx;
//     //   // uint8_t* lcs = (uint8_t*)malloc(lcs_len);
//     //   // // verbose("lcs: " + std::string(read->seq.s).substr(0,lcs_len));
//     //   // // Convert A,C,G,T,N into 0,1,2,3,4
//     //   // // The left context is reversed
//     //   // for (size_t i = 0; i < lcs_len; ++i)
//     //   //   lcs[lcs_len -i -1] = seq_nt4_table[(int)read->seq.s[i]];

//     //   // // rcs: right context sequence
//     //   // size_t rcs_occ = (mems[chain_anchors.back().first].idx + mems[chain_anchors.back().first].len); // The first character of the right context
//     //   // size_t rcs_len = read->seq.l - rcs_occ;
//     //   // uint8_t* rcs = (uint8_t*)malloc(rcs_len);
//     //   // // verbose("rcs: " + std::string(read->seq.s).substr(rcs_occ,rcs_len));
//     //   // // Convert A,C,G,T,N into 0,1,2,3,4
//     //   // for (size_t i = 0; i < rcs_len; ++i) 
//     //   //   rcs[i] = seq_nt4_table[(int)read->seq.s[ rcs_occ + i]];

//     //   // int32_t min_score = 20 + 8 * log(read->seq.l);
//     //   // // Fill between MEMs
//     //   // int32_t score = fill_chain(
//     //   //     mems,
//     //   //     chain_anchors,
//     //   //     lcs,   // Left context of the read
//     //   //     lcs_len, // Left context of the read lngth
//     //   //     rcs,   // Right context of the read
//     //   //     rcs_len, // Right context of the read length
//     //   //     read
//     //   //   );

//     //   //   if(score > min_score)
//     //   //   {
//     //   //     fill_chain(mems,chain_anchors,lcs,lcs_len,rcs,rcs_len,read,false,0,strand,out); 
//     //   //     aligned = true;
//     //   //   }

//     //   //   delete rcs;
//     //   //   delete lcs;
      
//     // }

//     /****************************************************************************/
//     /**** TODO: Tailored dynamic programming
//     // Sort MEM occurrences
//     for(auto mem: mems)
//       std::sort(mem.occs.begin(), mem.occs.end());

//     // Find non-overlapping MEM chains
//     std::vector<std::vector<size_t>> chains;
//     // first: chain idx
//     // second: mems idx
//     std::queue<std::pair<size_t,size_t>> next;
//     for(size_t i = 0; i < mems.size(); ++i)
//     {
//       chains.push_back(std::vector<size_t>());
//       next.push(make_pair(i,i)); 
//     }

//     while(not next.empty())
//     {
//       auto e = next.front(); next.pop();
//       chains[e.first].push_back(e.second);
//       // Find next non-overlapping MEM
//       size_t mem_end = mems[e.second].idx + mems[e.second].len;
//       size_t i = e.second + 1;
//       while( i < mems.size())
//       {
//         if(mems[i].idx >= mem_end)
//         {
//           chains.push_back(std::vector<size_t>(chains[e.first]));
//           next.push(make_pair(chains.size()-1,i));
//         }
//         ++i;
//       }
//     }
//     // Debug
//     verbose("New chains: ",chains.size());
//     for(auto chain: chains)
//     {
//       for(auto mem: chain)
//         verbose(mems[mem].idx, mems[mem].len, mems[mem].occs.size());
//       verbose("new chain");
//     }

//     // Reverse chains order
//     std::reverse(chains.begin(),chains.end());

//     // Find compatible occurrences of the MEMs in the chain
//     size_t max_gap = 10;

//     for(size_t i = 0; i< chains.size(); ++i)
//     {
//       if(chains[i].size() < 2)
//         continue;

//       auto& chain = chains[i];

//       std::vector<size_t> j(chain.size(),0); // Index of the occurrences

//       std:vector<std::vector<size_t>> valid();
//       // Lambda helpers
//       // Compute the distance between the i-th and i+1-th MEM in the chain
//       auto mem_distance = [] (size_t i) -> size_t {
//         return mems[chain[i+1]].idx - mems[chain[i]].idx;
//       };
//       // Compute the distance between the i-th and i+1-th MEM occurrences in j
//       auto occ_distance = [] (size_t i) -> size_t {
//         return mems[chain[i+1]].occs[j[i+1]] - mems[chain[i]].occs[j[i]];
//       };
//       // Increment j indexes
//       auto inc_j = [&] () -> bool {
        
//       };

      

//     }
//     ********/




//     // The distance between the MEMs in the read
    

//     // // Align the read
//     // if (mem_len >= min_len)
//     // {
//     //   // // Extract all the occurrences of the MEM
//     //   // std::vector<size_t> occs;
//     //   // occs.push_back(mem_pos);

//     //   // // Phi direction
//     //   // size_t curr = mem_pos;
//     //   // size_t next = ms.Phi(curr);
//     //   // size_t lcp =  lceToRBounded(ra,curr,next,mem_len);
//     //   // while(lcp >= mem_len)
//     //   // {
//     //   //   occs.push_back(next);
        
//     //   //   curr = next;
//     //   //   next = ms.Phi(curr);
//     //   //   lcp  = lceToRBounded(ra,curr,next,mem_len);
//     //   //   // verbose("Phi: " + std::to_string(lcp));
//     //   //   // if(occs.size() > 100)
//     //   //   //   error("More than 100 occs Phi" + std::string(read->seq.s));
//     //   // }

//     //   // // Phi_inv direction
//     //   // curr = mem_pos;
//     //   // next = ms.Phi_inv(curr);
//     //   // lcp =  lceToRBounded(ra,curr,next,mem_len);
//     //   // while(lcp >= mem_len)
//     //   // {
//     //   //   occs.push_back(next);
        
//     //   //   curr = next;
//     //   //   next = ms.Phi_inv(curr);
//     //   //   lcp  = lceToRBounded(ra,curr,next,mem_len);
//     //   //   // verbose("Phi_inv: " + std::to_string(next));
//     //   //   // if(occs.size() > 100)
//     //   //   //   error("More than 100 occs Phi_inv" + std::string(read->seq.s));
//     //   // }

//     //   // Extractin left and right context of the read
//     //   // lcs: left context sequence
//     //   size_t lcs_len = mem_idx;
//     //   uint8_t* lcs = (uint8_t*)malloc(lcs_len);
//     //   // verbose("lcs: " + std::string(read->seq.s).substr(0,lcs_len));
//     //   // Convert A,C,G,T,N into 0,1,2,3,4
//     //   // The left context is reversed
//     //   for (size_t i = 0; i < lcs_len; ++i)
//     //     lcs[lcs_len -i -1] = seq_nt4_table[(int)read->seq.s[i]];

//     //   // rcs: right context sequence
//     //   size_t rcs_occ = (mem_idx + mem_len); // The first character of the right context
//     //   size_t rcs_len = read->seq.l - rcs_occ;
//     //   uint8_t* rcs = (uint8_t*)malloc(rcs_len);
//     //   // verbose("rcs: " + std::string(read->seq.s).substr(rcs_occ,rcs_len));
//     //   // Convert A,C,G,T,N into 0,1,2,3,4
//     //   for (size_t i = 0; i < rcs_len; ++i) 
//     //     rcs[i] = seq_nt4_table[(int)read->seq.s[ rcs_occ + i]];

//     //   int32_t min_score = 20 + 8 * log(read->seq.l);
//     //   // verbose("Number of occurrences: " + std::to_string(occs.size()));
//     //   // For all the occurrences align
//     //   for(auto curr_mem_pos: occs)
//     //   {
//     //     int32_t score = extend(
//     //       curr_mem_pos,
//     //       mem_len,
//     //       lcs,   // Left context of the read
//     //       lcs_len, // Left context of the read lngth
//     //       rcs,   // Right context of the read
//     //       rcs_len // Right context of the read length
//     //     );

//     //     if(score > min_score)
//     //     {
//     //       extend(curr_mem_pos,mem_len,lcs,lcs_len,rcs,rcs_len,false,0,min_score,read,strand,out); 
//     //       aligned = true;
//     //     }

//     //   }
//     //   delete lcs;
//     //   delete rcs;
//     // }

//     if(not aligned)
//     {
//       std::string dummy = "";
//       write_sam(0, 0, 0, 0, "*", read, 0, out, dummy, dummy, 0, "*", 0, 0);
//     }

//     return aligned;
//   }

//   typedef struct{
//     bool aligned = false;
//     enum type_t {Unpaired, Paired} type;

//     kseq_t* read;
//     kseq_t* mate1;
//     kseq_t* mate2;




//   } alignment_t;

//   // Aligning pair-ended sequences
//   bool align(kseq_t *mate1, kseq_t *mate2, FILE *out)
//   {

//     bool aligned = false;
    
//     MTIME_START(0); // Timing helper

//     // Generate rc reads
//     kseq_t mate1_rev, mate2_rev;
//     rc_copy_kseq_t(&mate1_rev, mate1);
//     rc_copy_kseq_t(&mate2_rev, mate2);

//     // Find MEMs
//     std::vector<mem_t> mems;

//     find_mems(mate1, mems, 0, MATE_1 | MATE_F);
//     find_mems(&mate1_rev, mems, mate2->seq.l, MATE_1 | MATE_RC );
//     find_mems(mate2, mems, 0, MATE_2 | MATE_F);
//     find_mems(&mate2_rev, mems, mate1->seq.l, MATE_2 | MATE_RC);

//     MTIME_END(0); //Timing helper
//     MTIME_START(1); //Timing helper

//     std::vector<std::pair<size_t, size_t>> anchors;
//     std::vector<chain_t> chains;

//     // Chain MEMs
//     const bool chained = find_chains(mems, anchors, chains);

//     MTIME_END(1);   //Timing helper
//     MTIME_START(2); //Timing helper

//     if (not chained)
//     {
//       sam_t sam_m1(mate1);
//       sam_t sam_m2(mate2);

//       // Fill sam fields RNEXT, PNEXT and TLEN
//       sam_m1.rnext = std::string(mate2->name.s);
//       sam_m2.rnext = std::string(mate1->name.s);

//       sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;
      
//       write_sam(out, sam_m1);
//       write_sam(out, sam_m2);
//       return false;
//     }

//     int32_t min_score = 20 + 8 * log(mate1->seq.l);

//     // // Compute the second best score
//     std::vector<std::pair<int32_t, size_t>> best_scores;
//     // get the occurrences of the top 4 best scores
//     std::set<size_t> different_scores;
//     size_t i = 0;
//     while (i < chains.size() and different_scores.size() < 3)
//     {
//         different_scores.insert(chains[i].score);
//         if (different_scores.size() < 3)
//         {
//           // Align the chain
//           auto chain = chains[i];
//           // Reverse the chain order
//           std::reverse(chain.anchors.begin(), chain.anchors.end());
//           // Compute the score of a chain.
//           paired_score_t score;
//           if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
//             score = paired_chain_score(chain, anchors, mems, min_score, mate1, &mate2_rev);
//           else
//             score = paired_chain_score(chain, anchors, mems, min_score, &mate1_rev, mate2);
          
//           best_scores.push_back(std::make_pair(score.tot, i++));
//         }
//     }
    
//     if (best_scores.size() < 2)
//       best_scores.push_back(std::make_pair(0, chains.size()));

//     std::sort(best_scores.begin(), best_scores.end(), std::greater<std::pair<int32_t, size_t>>());

//     assert(best_scores.size() > 1);

//     if (best_scores[0].first < min_score)
//     {
//       sam_t sam_m1(mate1);
//       sam_t sam_m2(mate2);

//       // Fill sam fields RNEXT, PNEXT and TLEN
//       sam_m1.rnext = std::string(mate2->name.s);
//       sam_m2.rnext = std::string(mate1->name.s);

//       sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;

//       write_sam(out, sam_m1);
//       write_sam(out, sam_m2);
//       return false;
//     }

//     int32_t score2 = best_scores[1].first;
//     paired_score_t score;


//     { // Forward case
//       i = best_scores[0].second;
//       // Align the chain
//       auto chain = chains[i];
//       // Reverse the chain order
//       std::reverse(chain.anchors.begin(), chain.anchors.end());
//       // Compute the score of a chain.
//       if( (chain.mate == 0) || ((chain.mate & MATE_RC) and (chain.mate & MATE_2)) )
//         score = paired_chain_score(chain, anchors, mems, min_score, mate1, &mate2_rev, false, score2, 0, out);
//       else
//         score = paired_chain_score(chain, anchors, mems, min_score, &mate1_rev, mate2, false, score2, 1, out);
//     }

//     if (score.tot >= min_score)
//       aligned = true;

//     if (not aligned)
//     {
//       sam_t sam_m1(mate1);
//       sam_t sam_m2(mate2);

//       // Fill sam fields RNEXT, PNEXT and TLEN
//       sam_m1.rnext = std::string(mate2->name.s);
//       sam_m2.rnext = std::string(mate1->name.s);

//       sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;

//       write_sam(out, sam_m1);
//       write_sam(out, sam_m2);
//       return false;
//     }

//     MTIME_END(2); //Timing helper

//     return aligned;
//   }
//   // // Aligning pair-ended sequences
//   // bool align(kseq_t *mate1, kseq_t *mate2, FILE *out)
//   // {

//   //   bool aligned = false;
    
//   //   MTIME_START(0); // Timing helper

//   //   // Find MEMs of mate1 and reverse of mate2
//   //   std::vector<mem_t> mems;

//   //   find_mems(mate1, mems);
//   //   //copy seq
//   //   kseq_t mate2_rev;
//   //   rc_copy_kseq_t(&mate2_rev, mate2);
//   //   // copy_kseq_t(&mate2_rev, mate2);

//   //   // for (size_t i = 0; i < mate2->seq.l; ++i)
//   //   //   mate2_rev.seq.s[i] = complement(mate2->seq.s[mate2->seq.l - i - 1]);

//   //   // if (mate2_rev.seq.m > mate2_rev.seq.l)
//   //   //   mate2_rev.seq.s[mate2_rev.seq.l] = 0;

//   //   find_mems(&mate2_rev, mems, 1, mate1->seq.l, 1);

//   //   MTIME_END(0); //Timing helper
//   //   MTIME_START(1); //Timing helper

//   //   // for(size_t i = 0; i < mems.size(); ++i){
//   //   //   std::cout << "MEM[" << i <<"]: \n";
//   //   //   std::cout << "    len:" << mems[i].len <<"\n";
//   //   //   std::cout << "    pos:" << mems[i].pos <<"\n";
//   //   //   std::cout << "    idx:" << mems[i].idx <<"\n";
//   //   //   std::cout << "   occs:" << mems[i].occs.size() <<"\n";
//   //   // }

//   //   std::vector<std::pair<size_t, size_t>> anchors;
//   //   std::vector<chain_t> chains;

//   //   // Chain MEMs of mate1 and reverse of mate2
//   //   const bool fwd_chains = find_chains(mems, anchors, chains);

//   //   MTIME_END(1);   //Timing helper
//   //   MTIME_START(0); //Timing helper

//   //   // Find MEMs of mate2 and reverse of mate1
//   //   std::vector<mem_t> mems_rev;

//   //   find_mems(mate2, mems_rev, 0, 0, 1);

//   //   //copy seq
//   //   kseq_t mate1_rev;
//   //   rc_copy_kseq_t(&mate1_rev, mate1);
//   //   // copy_kseq_t(&mate1_rev, mate1);

//   //   // for (size_t i = 0; i < mate1->seq.l; ++i)
//   //   //   mate1_rev.seq.s[i] = complement(mate1->seq.s[mate1->seq.l - i - 1]);

//   //   // if (mate1_rev.seq.m > mate1_rev.seq.l)
//   //   //   mate1_rev.seq.s[mate1_rev.seq.l] = 0;

//   //   find_mems(&mate1_rev, mems_rev, 1, mate1->seq.l);

//   //   MTIME_END(0);   //Timing helper
//   //   MTIME_START(1); //Timing helper

//   //   // for (size_t i = 0; i < mems_rev.size(); ++i)
//   //   // {
//   //   //   std::cout << "MEM[" << i << "]: \n";
//   //   //   std::cout << "    len:" << mems_rev[i].len << "\n";
//   //   //   std::cout << "    pos:" << mems_rev[i].pos << "\n";
//   //   //   std::cout << "    idx:" << mems_rev[i].idx << "\n";
//   //   //   std::cout << "   occs:" << mems_rev[i].occs.size() << "\n";
//   //   // }

//   //   std::vector<std::pair<size_t, size_t>> anchors_rev;
//   //   std::vector<chain_t> chains_rev;

//   //   // Chain MEMs of mate2 and reverse of mate1
//   //   const bool rev_chains = find_chains(mems_rev, anchors_rev, chains_rev);

//   //   MTIME_END(1); //Timing helper
//   //   MTIME_START(2); //Timing helper
//   //   if ((not fwd_chains) and (not rev_chains))
//   //   {
//   //     sam_t sam_m1(mate1);
//   //     sam_t sam_m2(mate2);

//   //     // Fill sam fields RNEXT, PNEXT and TLEN
//   //     sam_m1.rnext = std::string(mate2->name.s);
//   //     sam_m2.rnext = std::string(mate1->name.s);

//   //     sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;
      
//   //     write_sam(out, sam_m1);
//   //     write_sam(out, sam_m2);
//   //     return false;
//   //   }

//   //   int32_t min_score = 20 + 8 * log(mate1->seq.l);



//   //   // // Compute the second best score
//   //   std::vector<std::pair<int32_t, size_t>> best_scores;
//   //   // get the occurrences of the top 4 best scores
//   //   std::set<size_t> different_scores;
//   //   size_t i = 0;
//   //   size_t j = 0;
//   //   size_t off = chains.size();
//   //   while (i < chains.size() and j < chains_rev.size() and different_scores.size() < 3)
//   //   {
//   //     if (chains[i].score > chains_rev[j].score)
//   //     {
//   //       different_scores.insert(chains[i].score);
//   //       if (different_scores.size() < 3)
//   //       {
//   //         // Align the chain
//   //         auto chain = chains[i];
//   //         // Reverse the chain order
//   //         std::reverse(chain.anchors.begin(), chain.anchors.end());
//   //         // Compute the score of a chain.
//   //         paired_score_t score = paired_chain_score(chain, anchors, mems, min_score, mate1, &mate2_rev);
//   //         best_scores.push_back(std::make_pair(score.tot, i++));
//   //       }
//   //     }
//   //     else
//   //     {
//   //       different_scores.insert(chains_rev[j].score);
//   //       if (different_scores.size() < 3)
//   //       {
//   //         // Align the chain
//   //         auto chain = chains_rev[j];
//   //         // Reverse the chain order
//   //         std::reverse(chain.anchors.begin(), chain.anchors.end());
//   //         // Compute the score of a chain.
//   //         paired_score_t score = paired_chain_score(chain, anchors_rev, mems_rev, min_score, &mate1_rev, mate2);
//   //         best_scores.push_back(std::make_pair(score.tot, off + j++));
//   //       }
//   //     }
//   //   }
//   //   while (different_scores.size() < 3 and i < chains.size())
//   //   {
//   //     different_scores.insert(chains[i].score);
//   //     if (different_scores.size() < 3)
//   //     {
//   //       // Align the chain
//   //       auto chain = chains[i];
//   //       // Reverse the chain order
//   //       std::reverse(chain.anchors.begin(), chain.anchors.end());
//   //       // Compute the score of a chain.
//   //       paired_score_t score = paired_chain_score(chain, anchors, mems, min_score, mate1, &mate2_rev);
//   //       best_scores.push_back(std::make_pair(score.tot, i++));
//   //     }
//   //   }
//   //   while (different_scores.size() < 3 and j < chains_rev.size())
//   //   {
//   //     different_scores.insert(chains_rev[j].score);
//   //     if (different_scores.size() < 3)
//   //     {
//   //       // Align the chain
//   //       auto chain = chains_rev[j];
//   //       // Reverse the chain order
//   //       std::reverse(chain.anchors.begin(), chain.anchors.end());
//   //       // Compute the score of a chain.
//   //       paired_score_t score = paired_chain_score(chain, anchors_rev, mems_rev, min_score, &mate1_rev, mate2);
//   //       best_scores.push_back(std::make_pair(score.tot, off + j++));
//   //     }
//   //   }

//   //   if (best_scores.size() < 2)
//   //     best_scores.push_back(std::make_pair(0, chains.size() + chains_rev.size()));

//   //   std::sort(best_scores.begin(), best_scores.end(), std::greater<std::pair<int32_t, size_t>>());

//   //   assert(best_scores.size() > 1);

//   //   if (best_scores[0].first < min_score)
//   //   {
//   //     sam_t sam_m1(mate1);
//   //     sam_t sam_m2(mate2);

//   //     // Fill sam fields RNEXT, PNEXT and TLEN
//   //     sam_m1.rnext = std::string(mate2->name.s);
//   //     sam_m2.rnext = std::string(mate1->name.s);

//   //     sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;

//   //     write_sam(out, sam_m1);
//   //     write_sam(out, sam_m2);
//   //     return false;
//   //   }

//   //   int32_t score2 = best_scores[1].first;
//   //   paired_score_t score;

//   //   // TODO: Check if the mates are paired or not.
//   //   if (best_scores[0].second >= off)
//   //   { // Reverse case
//   //     j = best_scores[0].second - off;
//   //     // Align the chain
//   //     auto chain = chains_rev[j];
//   //     // Reverse the chain order
//   //     std::reverse(chain.anchors.begin(), chain.anchors.end());
//   //     // Compute the score of a chain.
//   //     score = paired_chain_score(chain, anchors_rev, mems_rev, min_score, &mate1_rev, mate2, false, score2, 1, out);
//   //   }
//   //   else
//   //   { // Forward case
//   //     i = best_scores[0].second;
//   //     // Align the chain
//   //     auto chain = chains[i];
//   //     // Reverse the chain order
//   //     std::reverse(chain.anchors.begin(), chain.anchors.end());
//   //     // Compute the score of a chain.
//   //     score = paired_chain_score(chain, anchors, mems, min_score, mate1, &mate2_rev, false, score2, 0, out);
//   //   }

//   //   if (score.tot >= min_score)
//   //     aligned = true;

//   //   if (not aligned)
//   //   {
//   //     sam_t sam_m1(mate1);
//   //     sam_t sam_m2(mate2);

//   //     // Fill sam fields RNEXT, PNEXT and TLEN
//   //     sam_m1.rnext = std::string(mate2->name.s);
//   //     sam_m2.rnext = std::string(mate1->name.s);

//   //     sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;

//   //     write_sam(out, sam_m1);
//   //     write_sam(out, sam_m2);
//   //     return false;
//   //   }

//   //   MTIME_END(2); //Timing helper

//   //   return aligned;
//   // }

//   typedef struct{
//     ll score = 0;
//     size_t mate = 2;
//     bool paired = false;
//     std::vector<size_t> anchors;
//   } chain_t;


//   // Given a set of mems, find the chains
//   bool find_chains(
//     const   std::vector<mem_t>& mems,
//             std::vector< std::pair< size_t, size_t > >& anchors,
//             std::vector<chain_t>& chains        
//     )
//   {
//     /************* minimap2 dynamic programming for mem chaining ***************/
//     /* https://github.com/lh3/minimap2/blob/master/chain.c */

//     // Sort anchors
//     // Lamda helper to sort the anchors
//     auto cmp = [&] (std::pair<size_t, size_t> i, std::pair<size_t, size_t> j) -> bool {
//       return (mems[i.first].occs[i.second] + mems[i.first].len - 1) < (mems[j.first].occs[j.second] + mems[j.first].len - 1);
//     };

//     // TODO: improve this initialization
//     size_t tot_mem_length = 0;

//     // std::vector< std::pair< size_t, size_t > > anchors;
//     for(size_t i = 0; i < mems.size(); ++i)
//     {
//       for(size_t j = 0; j < mems[i].occs.size(); ++j)
//         anchors.push_back(make_pair(i,j));
//       tot_mem_length +=  mems[i].len * mems[i].occs.size();
//     }

//     float avg_mem_length = (float)tot_mem_length / anchors.size();

//     std::sort(anchors.begin(),anchors.end(),cmp);

//     // Dynamic programming

//     // TODO: Parameters to be defined
//     const ll G = LLONG_MAX;
//     const ll max_dist_x = 500;//LLONG_MAX;
//     // const ll max_dist_x = 100;//LLONG_MAX;
//     const ll max_dist_y = 100;//LLONG_MAX;
//     const ll max_iter = 50;
//     const ll max_pred = 50;
//     const ll min_chain_score = 1;
//     const ll min_chain_length = 1;


//     std::vector<ll> f(anchors.size(),0); // Score ending in position i
//     std::vector<ll> p(anchors.size(),0); // Position of the next anchor giving the max when chained with the one in position i
//     std::vector<ll> msc(anchors.size(),0); // Max score up to position i
//     std::vector<ll> t(anchors.size(),0); // Stores i in position p[j], if i is chained with j. See heuristics in minimap2

//     ll lb = 0;
//     // For all the anchors
//     for(size_t i = 0 ; i < anchors.size(); ++i)
//     {
//       // Get anchor i
//       const auto a_i = anchors[i];
//       const mem_t mem_i = mems[a_i.first];
//       // const ll k_i = mem_i.occs[a_i.second];
//       const ll x_i = mem_i.occs[a_i.second] + mem_i.len - 1;
//       // const ll z_i = mem_i.idx;
//       const ll y_i = mem_i.rpos;
//       // const ll y_i = mem_i.idx + mem_i.len - 1;
//       const ll w_i = mem_i.len;
//       const size_t mate_i = mem_i.mate;

//       ll max_f = w_i;
//       ll max_j = -1;
//       size_t n_pred = 0;
//       // For all previous anchors
//       // Heuristics from minimap2 -> do not try more than 50 anchors
//       if(i - lb > max_iter) lb = i - max_iter;
//       for(ll j = i-1; j >= lb; --j)
//       {
//         const auto a_j = anchors[j];
//         const mem_t mem_j = mems[a_j.first];
//         const ll x_j = mem_j.occs[a_j.second] + mem_j.len - 1;
//         const ll y_j = mem_j.rpos;
//         // const ll y_j = mem_j.idx + mem_j.len - 1;
//         const size_t mate_j = mem_j.mate;

//         // Check if anchors are compatible
//         // if mate_i == mate_j they are from the same mate and same orientation
//         // if mate_i ^ mate_j == 3 they are from different mate and withdifferent orientation
//         // TODO: we can change to same orientation replacing 3 with 1
//         if((mate_i != mate_j) and ((mate_i ^ mate_j)!= 3)) continue;

//         // If the current anchor is too far, exit.
//         if(x_i > x_j + max_dist_x)
//         {
//           // j = lb - 1;
//           lb = j;
//           continue;
//         }
//         // // skip if incompatible
//         // if(k_i < x_j or z_i < y_j) continue;

//         const ll x_d = x_i - x_j;
//         const ll y_d = y_i - y_j;
//         const int32_t l = (y_d > x_d ? (y_d - x_d) : (x_d - y_d));
//         const uint32_t ilog_l = (l > 0? ilog2_32(l): 0);

//         if((mate_i == mate_j and y_j >= y_i) or max(y_d, x_d) > G)
//           continue;

//         const ll alpha = min(min(y_d,x_d),w_i);
//         // const ll beta = (l > 0? (ll)(.01 * l * avg_mem_length) + ilog_l >> 1 : 0); 
//         ll beta = 0;
//         if(mate_i != mate_j){
//           if (x_d == 0) ++beta; // possibly due to overlapping paired ends; give a minor bonus
//           else {
//             const int c_lin = (int)(l * .01 * avg_mem_length);
//             beta = c_lin < ilog_l ? c_lin : ilog_l;
//           }
//         }else{
//           beta = (l > 0? (ll)(.01 * l * avg_mem_length) + ilog_l >> 1 : 0);
//         }



//         // No gap scale as in minimap2
//         ll score = f[j] + (alpha - beta);

//         if( score > max_f )
//         {
//           max_f = score;
//           max_j = j;
//           if( n_pred > 0) --n_pred;
//         }
//         else // minimap2: If i is chained wth j, than chaining i with a predecessor of j does not improve the score
//           if (t[j] == i and (++n_pred > max_pred))            
//             break;
          
//         if(p[j] > 0) t[p[j]] = i;
//       }

//       f[i] = max_f;
//       p[i] = max_j;
//       if( max_j >= 0 and msc[max_j] > max_f)
//         msc[i] = msc[max_j];
//       else
//         msc[i] = max_f;
//     }

//     // Find the end positions of chains
//     memset(t.data(),0,sizeof(size_t) * t.size());
//     for(size_t i = 0; i < anchors.size(); ++i)
//       if(p[i] >= 0) t[p[i]] = 1;
    
//     size_t n_chains = 0;
//     for(size_t i = 0; i < anchors.size(); ++i)
//       if(t[i] == 0 and msc[i] > min_chain_score) n_chains ++;
    
//     // TODO: check if we want to report also non aligned reads
//     if(n_chains == 0)
//       return false;

//     // TODO: replace the vector of pairs with a lambda for the sort
//     std::vector<std::pair<ll,size_t>> chain_starts(n_chains); // Stores uniqe chains and their index of uniqe chains in anchors
//     size_t k = 0;
//     for(size_t i = 0; i < anchors.size(); ++i)
//     {
//       if(t[i] == 0 and  msc[i] > min_chain_score)
//       {
//         size_t j = i;
//         while(j >= 0 and f[j] < msc[j]) j = p[j]; // Find teh peak that maximizes f
//         if (j < 0) i = 0;
//         chain_starts[k++] = make_pair(f[j],j);
//       }
//     }

//     chain_starts.resize(k), chain_starts.shrink_to_fit();
//     n_chains = chain_starts.size();
//     std::sort(chain_starts.begin(), chain_starts.end(), std::greater<std::pair<ll,size_t>>());
    
//     // std::vector<std::pair<size_t, std::vector<size_t>>> chains;

//     // Backtrack
//     memset(t.data(),0,sizeof(size_t) * t.size());
//     for(size_t i = 0; i < n_chains; ++i)
//     {
//       ll j = chain_starts[i].second;
//       chain_t chain;
//       chain.mate = mems[anchors[j].first].mate;
//       chain.score = chain_starts[i].first;
//       do {
//         size_t tmp_mate = mems[anchors[j].first].mate;
//         chain.paired = chain.paired or (chain.mate != tmp_mate);
//         chain.anchors.push_back(j); // stores th reverse of the chain
//         t[j] = 1;
//         j = p[j];
//       } while (j >= 0 && t[j] == 0);
//       if (j < 0) { // l - prev_l is the length of the chain
//         if (chain.anchors.size() >= min_chain_length)
//           chains.push_back(chain);     
//       } else if (chain_starts[i].first - f[j] >= min_chain_score) { // Two chains share a common prefix
//         if (chain.anchors.size() >= min_chain_length)
//           chains.push_back(chain);
//       }
//     }

//     // Lamda helper to sort the anchors
//     auto chain_t_cmp = [](chain_t i, chain_t j) -> bool
//     {
//       return i.score > j.score;
//     };

//     // Sort the chains by max scores.
//     std::sort(chains.begin(), chains.end(), chain_t_cmp);
 
//     // // Backtrack
//     // memset(t.data(),0,sizeof(size_t) * t.size());
//     // for(size_t i = 0; i < n_chains; ++i)
//     // {
//     //   ll j = chain_starts[i].second;
//     //   std::vector<size_t> chain;
//     //   bool paired = true;
//     //   size_t mate = mems[anchors[j].first].mate; 
//     //   do {
//     //     size_t tmp_mate = mems[anchors[j].first].mate;
//     //     paired = paired and (mate == tmp_mate);
//     //     chain.push_back(j); // stores th reverse of the chain
//     //     t[j] = 1;
//     //     j = p[j];
//     //   } while (j >= 0 && t[j] == 0);
//     //   if (j < 0) { // l - prev_l is the length of the chain
//     //     if (chain.size() >= min_chain_length)
//     //       chains.push_back(std::make_pair(chain_starts[i].first, chain));     
//     //   } else if (chain_starts[i].first - f[j] >= min_chain_score) { // Two chains share a common prefix
//     //     if (chain.size() >= min_chain_length)
//     //       chains.push_back(std::make_pair((chain_starts[i].first - f[j]), chain));
//     //   }
//     // }

//     // // Sort the chains by max scores.
//     // std::sort(chains.begin(), chains.end(), std::greater<std::pair<ll,std::vector<size_t>>>());

//     // Clear space
//     p.resize(0), p.shrink_to_fit();
//     t.resize(0), t.shrink_to_fit();
//     f.resize(0), f.shrink_to_fit();
//     msc.resize(0), msc.shrink_to_fit();

//     return true;
//   }

//   void find_mems(
//     const kseq_t *read,
//     std::vector<mem_t>& mems,
//     size_t r_offset = 0,
//     size_t mate = 0
//     ) 
//   {
//     auto pointers = ms.query(read->seq.s, read->seq.l);
//     size_t l = 0;   // Current match length
//     size_t pl = 0;  // Previous match length
//     size_t n_Ns = 0;
//     for (size_t i = 0; i < pointers.size(); ++i)
//     {
//       size_t pos = pointers[i];
//       while ((i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
//       {
//         if(read->seq.s[i + l] == 'N') n_Ns++;
//         else n_Ns = 0;
//         ++l;
//       }

//       // Update MEMs
//       if (l >= pl and n_Ns < l and l >= min_len)
//       {
//         size_t r = r_offset + (i + l - 1); // compatible with minimap2 chaining algorithm
//         // size_t r = r_offset + ((reverse) ? (i) : (i + l - 1)); // compatible with minimap2 chaining algorithm
//         // size_t r = i + l - 1;
//         // r = (reverse)? (read->seq.l - (r + 1 - l) - 1) : (r); // compatible with minimap2 chaining algorithm
//         mems.push_back(mem_t(pointers[i],l,i,mate,r));
//         find_MEM_occs(mems.back());
//       }

//       // Compute next match length
//       pl = l;
//       l = (l == 0 ? 0 : (l - 1));
//     }

//   }

//   int32_t chain_score(
//       const std::pair<size_t, std::vector<size_t>> &chain,
//       const std::vector<std::pair<size_t, size_t>> &anchors,
//       const std::vector<mem_t> &mems,
//       const int32_t min_score,
//       const kseq_t *read,
//       const bool score_only = true,
//       const int32_t score2 = 0,
//       const uint8_t strand = 0,
//       FILE *out = nullptr,
//       sam_t* sam = nullptr)     // The SAM information pointer)
//   {
//     bool aligned = false;
//     // Extract the anchors
//     std::vector<std::pair<size_t,size_t>> chain_anchors(chain.second.size());
//     for(size_t i = 0; i < chain_anchors.size(); ++i)
//       chain_anchors[i] = anchors[chain.second[i]];
//     // Extracting left and right context of the read
//     // lcs: left context sequence
//     size_t lcs_len = mems[chain_anchors[0].first].idx;
//     uint8_t* lcs = (uint8_t*)malloc(lcs_len);
//     // verbose("lcs: " + std::string(read->seq.s).substr(0,lcs_len));
//     // Convert A,C,G,T,N into 0,1,2,3,4
//     // The left context is reversed
//     for (size_t i = 0; i < lcs_len; ++i)
//       lcs[lcs_len -i -1] = seq_nt4_table[(int)read->seq.s[i]];

//     // rcs: right context sequence
//     size_t rcs_occ = (mems[chain_anchors.back().first].idx + mems[chain_anchors.back().first].len); // The first character of the right context
//     size_t rcs_len = read->seq.l - rcs_occ;
//     uint8_t* rcs = (uint8_t*)malloc(rcs_len);
//     // verbose("rcs: " + std::string(read->seq.s).substr(rcs_occ,rcs_len));
//     // Convert A,C,G,T,N into 0,1,2,3,4
//     for (size_t i = 0; i < rcs_len; ++i) 
//       rcs[i] = seq_nt4_table[(int)read->seq.s[ rcs_occ + i]];

//     // Fill between MEMs
//     int32_t score = fill_chain(
//         mems,
//         chain_anchors,
//         lcs,   // Left context of the read
//         lcs_len, // Left context of the read lngth
//         rcs,   // Right context of the read
//         rcs_len, // Right context of the read length
//         read
//       );

//       if(!score_only and score >= min_score)
//       {
//         bool output = (sam == nullptr);

//         if(output)
//           sam = new sam_t(read);

//         fill_chain(mems,chain_anchors,lcs,lcs_len,rcs,rcs_len,read,score_only,sam); 
        
//         sam->flag = (strand?16:0);
//         sam->zs = score2;
//         sam->mapq = compute_mapq(sam->as, sam->zs, min_score, sam->read->seq.l);
//         sam->reverse = (strand != 0);
//         sam->rname = idx[sam->pos - 1];
//         if(output)
//         {
//           write_sam(out,*sam);
//           delete sam;
//         }
        
//         // fill_chain(mems,chain_anchors,lcs,lcs_len,rcs,rcs_len,read,score_only,0,min_score,strand,out); 
//       }

//       delete rcs;
//       delete lcs;

//       return score;
//   }

//   typedef struct{
//     int32_t tot = 0;
//     int32_t m1 = 0;
//     int32_t m2 = 0; 
//     bool paired = false;
//   } paired_score_t;

//   paired_score_t paired_chain_score(
//       const chain_t &chain,
//       const std::vector<std::pair<size_t, size_t>> &anchors,
//       const std::vector<mem_t> &mems,
//       const int32_t min_score,
//       const kseq_t *mate1,
//       const kseq_t *mate2,
//       const bool score_only = true,
//       const int32_t score2 = 0,
//       const uint8_t strand = 0,
//       FILE *out = nullptr)
//   {
//     paired_score_t score;
//     if(chain.paired){

//       // Extract the anchors
//       std::pair<ll, std::vector<size_t>> mate1_chain;
//       std::pair<ll, std::vector<size_t>> mate2_chain;
//       for(size_t i = 0; i < chain.anchors.size(); ++i)
//       {
//         size_t anchor_id = chain.anchors[i];
//         if ((mems[anchors[anchor_id].first].mate & MATE_2) == 0)
//           mate1_chain.second.push_back(anchor_id);
//         else
//           mate2_chain.second.push_back(anchor_id);
//       }

//       score.paired = chain.paired;

//       if(score_only)
//       {
//           score.m1 = chain_score(mate1_chain, anchors, mems, min_score, mate1);
//           score.m2 = chain_score(mate2_chain, anchors, mems, min_score, mate2);
//       }
//       else
//       {
//         sam_t sam_m1(mate1);
//         sam_t sam_m2(mate2);


//         score.m1 = chain_score(mate1_chain, anchors, mems, min_score, mate1, false, score2, strand, nullptr, &sam_m1);
//         score.m2 = chain_score(mate2_chain, anchors, mems, min_score, mate2, false, score2, strand, nullptr, &sam_m2);

//         // Fill sam fields RNEXT, PNEXT and TLEN
//         sam_m1.rnext = std::string(mate2->name.s);
//         sam_m2.rnext = std::string(mate1->name.s);

//         if(score.m1 >= min_score and score.m2 >= min_score)
//         {

//           sam_m1.pnext = sam_m2.pos;
//           sam_m2.pnext = sam_m1.pos;

//           ll tlen = (sam_m2.pos + mate2->seq.l) - sam_m1.pos;

//           sam_m1.tlen = tlen;
//           sam_m2.tlen = -tlen;

//           sam_m1.rname = idx[sam_m1.pos - 1];
//           sam_m2.rname = idx[sam_m2.pos - 1];

//           sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_MAPPED_PAIRED;
//           if(strand)
//           {
//             sam_m1.flag |= SAM_REVERSED | SAM_FIRST_IN_PAIR;
//             sam_m2.flag |= SAM_MATE_REVERSED | SAM_SECOND_IN_PAIR;
//           }
//           else
//           {
//             sam_m1.flag |= SAM_MATE_REVERSED | SAM_FIRST_IN_PAIR;
//             sam_m2.flag |= SAM_REVERSED | SAM_SECOND_IN_PAIR;
//           }
//         }else if(score.m1 >= min_score) {
//           sam_m1.rname = idx[sam_m1.pos - 1];

//           sam_m1.flag = SAM_PAIRED | SAM_MATE_UNMAPPED | SAM_FIRST_IN_PAIR;
//           sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_SECOND_IN_PAIR;
//           if(strand)
//             sam_m1.flag |= SAM_REVERSED;
//         }else if(score.m2 >= min_score) {
//           sam_m2.rname = idx[sam_m2.pos - 1];

//           sam_m1.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_FIRST_IN_PAIR;
//           sam_m2.flag = SAM_PAIRED | SAM_MATE_UNMAPPED | SAM_SECOND_IN_PAIR;
//           if(not strand)
//             sam_m2.flag |= SAM_REVERSED;
//         }else {
//           sam_m1.flag = sam_m2.flag = SAM_PAIRED | SAM_UNMAPPED | SAM_MATE_UNMAPPED;
//         }

//         write_sam(out,sam_m1);
//         write_sam(out,sam_m2);
//       }
//     }
//     else
//     {
//       if(chain.mate == 0)
//       {
//         // score.m1 = chain_score(mate1_chain, anchors, mems, min_score, mate1, false, score2, strand, nullptr, &sam_m1);
//       }
//       else
//       {
//         // score.m2 = chain_score(mate2_chain, anchors, mems, min_score, mate2, false, score2, strand, nullptr, &sam_m2);

//       }
//     }


//     score.tot = score.m1 + score.m2;
//     return score;
//   }

//   size_t get_aligned_reads()
//   {
//     return aligned_reads;
//   }

//   // If score_only is true we compute the score of the alignment. 
//   // If score_only is false, we extend again the read and we write the result
//   // in the SAM file, so we need to give the second best score.
//   int32_t extend(
//       const size_t mem_pos,
//       const size_t mem_len,
//       const uint8_t *lcs,           // Left context of the read
//       const size_t lcs_len,         // Left context of the read lngth
//       const uint8_t *rcs,           // Right context of the read
//       const size_t rcs_len,         // Right context of the read length
//       const bool score_only = true, // Report only the score
//       const int32_t score2 = 0,     // The score of the second best alignment
//       const int32_t min_score = 0,  // The minimum score to call an alignment
//       const kseq_t *read = nullptr, // The read that has been aligned
//       int8_t strand = 0,            // 0: forward aligned ; 1: reverse complement aligned
//       FILE *out = nullptr,          // The SAM file pointer
//       const bool realign = false    // Realign globally the read
//   )
//   {
//     int flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT;

//     if(score_only) 
//       flag = KSW_EZ_SCORE_ONLY;

//     int score_lc = 0;
//     int score_rc = 0;

//     // TODO: Update end_bonus according to the MEM contribution to the score
//     ksw_extz_t ez_lc;
//     ksw_extz_t ez_rc;
//     ksw_extz_t ez;
//     memset(&ez_lc, 0, sizeof(ksw_extz_t));
//     memset(&ez_rc, 0, sizeof(ksw_extz_t));
//     memset(&ez, 0, sizeof(ksw_extz_t));

//     // Extract the context from the reference
//     // lc: left context
//     ksw_reset_extz(&ez_lc);
//     if(lcs_len > 0)
//     {
//       size_t lc_occ = (mem_pos > 100 ? mem_pos - 100 : 0);
//       size_t lc_len = (mem_pos > 100 ? 100 : 100 - mem_pos);
//       char *tmp_lc = (char *)malloc(100);
//       ra.expandSubstr(lc_occ, lc_len, tmp_lc);
//       // verbose("lc: " + std::string(lc));
//       // Convert A,C,G,T,N into 0,1,2,3,4
//       // The left context is reversed
//       uint8_t *lc = (uint8_t *)malloc(100);
//       for (size_t i = 0; i < lc_len; ++i)
//         lc[lc_len -i -1] = seq_nt4_table[(int)tmp_lc[i]];
      
//       delete tmp_lc;

//       // Query: lcs
//       // Target: lc
//       // verbose("aligning lc and lcs");
//       ksw_extz2_sse(km, lcs_len, (uint8_t*)lcs, lc_len, (uint8_t*)lc, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez_lc);
//       score_lc =  ez_lc.mqe;
//       // verbose("lc score: " + std::to_string(score_lc));
//       // Check if the extension reached the end or the query
//       assert(score_only or ez_lc.reach_end);

//       // std::string blc = print_BLAST_like((uint8_t*)lc,(uint8_t*)lcs,ez_lc.cigar,ez_lc.n_cigar);
//       // std::cout<<blc;

//       delete lc;
//     }

//     // rc: right context
//     ksw_reset_extz(&ez_rc);
//     if(rcs_len > 0)
//     {
//       size_t rc_occ = mem_pos + mem_len;
//       size_t rc_len = (rc_occ < n-100 ? 100 : n - rc_occ);
//       char *rc = (char *)malloc(100);
//       ra.expandSubstr(rc_occ, rc_len, rc);
//       // verbose("rc: " + std::string(rc));
//       // Convert A,C,G,T,N into 0,1,2,3,4
//       for (size_t i = 0; i < rc_len; ++i)
//         rc[i] = seq_nt4_table[(int)rc[i]];

//       // Query: rcs
//       // Target: rc
//       // verbose("aligning rc and rcs");
//       ksw_extz2_sse(km, rcs_len, (uint8_t*)rcs, rc_len, (uint8_t*)rc, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez_rc);
//       score_rc = ez_rc.mqe;
//       // verbose("rc score: " + std::to_string(score_rc));
//       // Check if the extension reached the end or the query
//       assert(score_only or ez_rc.reach_end);

//       // std::string brc = print_BLAST_like((uint8_t*)rc,(uint8_t*)rcs,ez_rc.cigar,ez_rc.n_cigar);
//       // std::cout<<brc;
//       delete rc;
//     }


//     // Compute the final score
//     int32_t score = mem_len * smatch + score_lc + score_rc;

//     if(not score_only)
//     {
//       // Compute starting position in reference
//       size_t ref_pos = mem_pos - (lcs_len > 0 ? ez_lc.mqe_t + 1 : 0);
//       size_t ref_len = (lcs_len > 0 ? ez_lc.mqe_t + 1 : 0) + mem_len + (rcs_len > 0 ? ez_rc.mqe_t + 1: 0);
//       char *ref = (char *)malloc(ref_len);
//       ra.expandSubstr(ref_pos, ref_len, ref);
//       // Convert A,C,G,T,N into 0,1,2,3,4
//       for (size_t i = 0; i < ref_len; ++i)
//         ref[i] = seq_nt4_table[(int)ref[i]];

//       // Convert the read
//       size_t seq_len = read->seq.l;
//       uint8_t* seq = (uint8_t*) malloc(seq_len);
//       for (size_t i = 0; i < seq_len; ++i)
//         seq[i] = seq_nt4_table[(int)read->seq.s[i]];

//       char* tmp = (char*)calloc(max(ref_len,seq_len),1);

//       if(realign)
//       {
//         // Realign the whole sequence globally
//         flag = KSW_EZ_RIGHT;
//         ksw_reset_extz(&ez);
//         ksw_extz2_sse(km, seq_len, (uint8_t*)seq, ref_len, (uint8_t*)ref, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez);


//         // std::string bfull = print_BLAST_like((uint8_t*)ref,seq,ez.cigar,ez.n_cigar);
//         // std::cout << bfull;

//         // Example were ez.score is lrger than score:
//         // Left context alignment
//         // 22333022233022233302223302223
//         // ||||  ||||||||| ||||||*|||||*
//         // 2233  222330222 3302220302220
//         // Right context alignment
//         // 33022233022233022233022233022233      0222334
//         // *||||||||||||||||*||||||||||||||      ||||||*
//         // 130222330222330220330222330222332222330222330
//         // [INFO] 16:26:16 - Message: old score:  130  new score:  140
//         // Global alignment
//         // 2233    3022233022233302223  30222330222330222330222330222  330222  330222330222330222330222330222330222334
//         // ||||    |||||||||||*| ||||*  |||||||||||||||||||||||||||||  *|||||  |||||||||||*||||||||||||||*|||||||||||*
//         // 223322233022233022203 02220  30222330222330222330222330222  130222  330222330220330222330222332222330222330
//         // The original occurrence of the MEM has been shifted to the left by 6 positions, 
//         // reducing the gap in the right context, and moving in to the left context.

//         assert(ez.score >= score);

//         // Concatenate the CIGAR strings
        
//         std::string cigar_s;
//         for(size_t i = 0; i < ez.n_cigar; ++i)
//           cigar_s += std::to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];

//         // Compute the MD:Z field and thenumber of mismatches
//         std::string mdz_s;
//         size_t nm = write_MD_core((uint8_t*)ref,seq,ez.cigar,ez.n_cigar,tmp,0,mdz_s);

//         write_sam(ez.score, score2, min_score, ref_pos, idx[ref_pos].c_str(), read, strand, out, cigar_s, mdz_s, nm, "*", 0, 0);
//       }
//       else
//       {
//         // Concatenate the CIGAR strings
//         size_t n_cigar = ez_lc.n_cigar + ez_rc.n_cigar + 1;
//         uint32_t *cigar = (uint32_t*)calloc(n_cigar,sizeof(uint32_t));
//         size_t i = 0;

//         for(size_t j = 0; j < ez_lc.n_cigar; ++j)
//           cigar[i++] = ez_lc.cigar[ez_lc.n_cigar -j -1];


//         if(ez_lc.n_cigar > 0 and ((cigar[i-1]& 0xf) == 0))
//         { // If the previous operation is also an M then merge the two operations
//           cigar[i-1] += (((uint32_t)mem_len) << 4);
//           --n_cigar;
//         }
//         else
//           cigar[i++] = (((uint32_t)mem_len) << 4);


//         if(ez_rc.n_cigar > 0)
//         {
//           if((ez_rc.cigar[0]& 0xf) == 0)
//           { // If the next operation is also an M then merge the two operations
//             cigar[i-1] += ez_rc.cigar[0];
//             --n_cigar;
//           }
//           else
//             cigar[i++] = ez_rc.cigar[0];
//         }

//         for(size_t j = 1; j < ez_rc.n_cigar; ++j)
//           cigar[i++] = ez_rc.cigar[j];
        
//         assert(i <= n_cigar);

//         // std::string bfull = print_BLAST_like((uint8_t*)ref,seq,cigar,n_cigar);
//         // std::cout << bfull;

//         std::string cigar_s;
//         for(size_t i = 0; i < n_cigar; ++i)
//           cigar_s += std::to_string(cigar[i] >> 4) + "MID"[cigar[i] & 0xf];


//         // Compute the MD:Z field and thenumber of mismatches
//         std::string mdz_s;
//         size_t nm = write_MD_core((uint8_t*)ref,seq,cigar,n_cigar,tmp,0,mdz_s);

//         write_sam(score, score2, min_score, ref_pos, idx[ref_pos].c_str(), read, strand, out, cigar_s, mdz_s, nm, "*", 0, 0);

//         delete cigar;
//       }
//       delete tmp;
//       delete ref;
//       delete seq;
//     }

//     if (ez_lc.m_cigar > 0)
//       delete ez_lc.cigar;
//     if (ez_rc.m_cigar > 0)
//       delete ez_rc.cigar;
//     if (ez.m_cigar > 0)
//       delete ez.cigar;

//     return score;
//   }

//   // If score_only is true we compute the score of the alignment. 
//   // If score_only is false, we extend again the read and we write the result
//   // in the SAM file, so we need to give the second best score. 
//   int32_t fill_chain(
//     const std::vector<mem_t>& mems, // pairs of pos/lengths of the mems
//     const std::vector<std::pair<size_t,size_t>>& anchors, // pairs of mem index/occ index
//     const uint8_t* lcs,   // Left context of the read
//     const size_t lcs_len, // Left context of the read lngth
//     const uint8_t* rcs,   // Right context of the read
//     const size_t rcs_len, // Right context of the read length
//     const kseq_t *read, // The read that has been aligned
//     const bool score_only = true, // Report only the score
//     sam_t* sam = nullptr,     // The SAM information pointer
//     // const int32_t score2 = 0,    // The score of the second best alignment
//     // const int32_t min_score = 0, // The minimum score to call an alignment
//     // const int8_t strand = 0,    // 0: forward aligned ; 1: reverse complement aligned
//     // const FILE *out = nullptr,   // The SAM file pointer
//     bool realign = false   // Realign globally the read
//     // const bool realign = false   // Realign globally the read
//   )
//   {
//     int flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT;

//     if(score_only) 
//       flag = KSW_EZ_SCORE_ONLY;

//     int score_lc = 0;
//     int score_rc = 0;

//     ksw_extz_t ez_lc;
//     ksw_extz_t ez_rc;
//     ksw_extz_t ez;
//     memset(&ez_lc, 0, sizeof(ksw_extz_t));
//     memset(&ez_rc, 0, sizeof(ksw_extz_t));
//     memset(&ez, 0, sizeof(ksw_extz_t));
//     // TODO: Update end_bonus according to the MEM contribution to the score
    
//     // Extract the context from the reference
//     // lc: left context of the first mem
//     ksw_reset_extz(&ez_lc);
//     if(lcs_len > 0)
//     {
//       size_t mem_pos = mems[anchors[0].first].occs[anchors[0].second];

//       size_t lc_occ = (mem_pos > 100 ? mem_pos - 100 : 0);
//       size_t lc_len = (mem_pos > 100 ? 100 : 100 - mem_pos);
//       char *tmp_lc = (char *)malloc(100);
//       ra.expandSubstr(lc_occ, lc_len, tmp_lc);
//       // verbose("lc: " + std::string(lc));
//       // Convert A,C,G,T,N into 0,1,2,3,4
//       // The left context is reversed
//       uint8_t *lc = (uint8_t *)malloc(100);
//       for (size_t i = 0; i < lc_len; ++i)
//         lc[lc_len -i -1] = seq_nt4_table[(int)tmp_lc[i]];
      
//       delete tmp_lc;

//       // Query: lcs
//       // Target: lc
//       // verbose("aligning lc and lcs");
//       ksw_extz2_sse(km, lcs_len, (uint8_t*)lcs, lc_len, (uint8_t*)lc, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez_lc);
//       score_lc =  ez_lc.mqe;
//       // verbose("lc score: " + std::to_string(score_lc));
//       // Check if the extension reached the end or the query
//       assert(score_only or ez_lc.reach_end);

//       // std::string blc = print_BLAST_like((uint8_t*)lc,(uint8_t*)lcs,ez_lc.cigar,ez_lc.n_cigar);
//       // std::cout<<blc;

//       delete lc;
//     }

//     // rc: right context of the last mem
//     ksw_reset_extz(&ez_rc);
//     if(rcs_len > 0)
//     {
//       size_t mem_pos =  mems[anchors.back().first].occs[anchors.back().second];
//       size_t mem_len =  mems[anchors.back().first].len;

//       size_t rc_occ = mem_pos + mem_len;
//       size_t rc_len = (rc_occ < n-100 ? 100 : n - rc_occ);
//       char *rc = (char *)malloc(100);
//       ra.expandSubstr(rc_occ, rc_len, rc);
//       // verbose("rc: " + std::string(rc));
//       // Convert A,C,G,T,N into 0,1,2,3,4
//       for (size_t i = 0; i < rc_len; ++i)
//         rc[i] = seq_nt4_table[(int)rc[i]];

//       // Query: rcs
//       // Target: rc
//       // verbose("aligning rc and rcs");
//       ksw_extz2_sse(km, rcs_len, (uint8_t*)rcs, rc_len, (uint8_t*)rc, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez_rc);
//       score_rc = ez_rc.mqe;
//       // verbose("rc score: " + std::to_string(score_rc));
//       // Check if the extension reached the end or the query
//       assert(score_only or ez_rc.reach_end);

//       // std::string brc = print_BLAST_like((uint8_t*)rc,(uint8_t*)rcs,ez_rc.cigar,ez_rc.n_cigar);
//       // std::cout<<brc;
//       delete rc;
//     }


//     // Compute the partial score score
//     int32_t score = score_lc + score_rc;
//     // int32_t score = mem_len * smatch + score_lc + score_rc;

//     // Compute starting position in reference
//     size_t mem_pos = mems[anchors[0].first].occs[anchors[0].second];
//     size_t mem_len = mems[anchors.back().first].occs[anchors.back().second] + mems[anchors.back().first].len - mem_pos; // from the strart of the first MEM to the end of the last MEM.

//     size_t ref_pos = mem_pos - (lcs_len > 0 ? ez_lc.mqe_t + 1 : 0);
//     size_t ref_len = (lcs_len > 0 ? ez_lc.mqe_t + 1 : 0) + mem_len + (rcs_len > 0 ? ez_rc.mqe_t + 1: 0);
//     char *ref = (char *)malloc(ref_len);
//     ra.expandSubstr(ref_pos, ref_len, ref);

//     // Convert A,C,G,T,N into 0,1,2,3,4
//     for (size_t i = 0; i < ref_len; ++i)
//       ref[i] = seq_nt4_table[(int)ref[i]];

//     // Convert the read
//     size_t seq_len = read->seq.l;
//     uint8_t* seq = (uint8_t*) malloc(seq_len);
//     for (size_t i = 0; i < seq_len; ++i)
//       seq[i] = seq_nt4_table[(int)read->seq.s[i]];

//     // TODO: fix the gap fill when MEMs overlap. The quick fix is a full alignment between the first and last MEM
//     // // Fill the gaps between each mem
//     // std::vector<ksw_extz_t> ez_cc(anchors.size()-1);
//     // for(size_t i = 1; i < anchors.size(); ++i)
//     // {
//     //   flag = KSW_EZ_RIGHT;
//     //   ksw_reset_extz(&ez_cc[i-1]);

//     //   // size_t cc_occ = mems[anchors[i-1].first].occs[anchors[i-1].second] + mems[anchors[i-1].first].len;
//     //   // size_t cc_len = mems[anchors[i].first].occs[anchors[i].second] - cc_occ; 
//     //   // char *cc = (char *)malloc(cc_len);
//     //   // ra.expandSubstr(cc_occ, cc_len, cc); // TODO: Replace all expand substrings with only one expand at the beginning.
//     //   // // verbose("rc: " + std::string(rc));
//     //   // // Convert A,C,G,T,N into 0,1,2,3,4
//     //   // for (size_t i = 0; i < cc_len; ++i)
//     //   //   cc[i] = seq_nt4_table[(int)cc[i]];

//     //   size_t cc_occ = mems[anchors[i-1].first].occs[anchors[i-1].second] + mems[anchors[i-1].first].len;
//     //   size_t cc_len = mems[anchors[i].first].occs[anchors[i].second] - cc_occ; 
//     //   cc_occ -= ref_pos;
//     //   char *cc = (char *)malloc(cc_len);
//     //   // Convert A,C,G,T,N into 0,1,2,3,4
//     //   for (size_t i = 0; i < cc_len; ++i)
//     //     cc[i] = ref[cc_occ + i];


//     //   size_t ccs_pos = mems[anchors[i-1].first].idx + mems[anchors[i-1].first].len;
//     //   size_t ccs_len = mems[anchors[i].first].idx - ccs_pos;

//     //   // Query: rcs
//     //   // Target: rc
//     //   // verbose("aligning rc and rcs");
//     //   ksw_extz2_sse(km, ccs_len, (uint8_t*)(seq + ccs_pos), cc_len, (uint8_t*)cc, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez_cc[i-1]);
//     //   // score_cc = ez_cc.mqe;
//     //   // verbose("rc score: " + std::to_string(score_rc));
//     //   // Check if the extension reached the end or the query
//     //   // assert(score_only or ez_cc[i-1].reach_end);

//     //   // Update score
//     //   score +=  mems[anchors[i-1].first].len * smatch + ez_cc[i-1].mqe;

//     //   // std::string brc = print_BLAST_like((uint8_t*)rc,(uint8_t*)rcs,ez_rc.cigar,ez_rc.n_cigar);
//     //   // std::cout<<brc;
//     //   delete cc;
//     // }

//     // score +=  mems[anchors.back().first].len * smatch;

// //******************************************************************************
// // BEGIN RAPID PROTOTYPING HACK
// //******************************************************************************
//     // Realign the whole sequence globally
//     ksw_reset_extz(&ez);
//     ksw_extz2_sse(km, seq_len, (uint8_t *)seq, ref_len, (uint8_t *)ref, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez);

//     score = ez.score;

//     realign = true;

// //******************************************************************************
// // END RAPID PROTOTYPING HACK
// //******************************************************************************
//     if(not score_only)
//     {
//       // Compute starting position in reference

//       char* tmp = (char*)calloc(max(ref_len,seq_len),1);

//       if(realign)
//       {
//         // Realign the whole sequence globally
//         flag = KSW_EZ_RIGHT;
//         ksw_reset_extz(&ez);
//         ksw_extz2_sse(km, seq_len, (uint8_t*)seq, ref_len, (uint8_t*)ref, m, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez);


//         // std::string bfull = print_BLAST_like((uint8_t*)ref,seq,ez.cigar,ez.n_cigar);
//         // std::cout << bfull;

//         // Example were ez.score is lrger than score:
//         // Left context alignment
//         // 22333022233022233302223302223
//         // ||||  ||||||||| ||||||*|||||*
//         // 2233  222330222 3302220302220
//         // Right context alignment
//         // 33022233022233022233022233022233      0222334
//         // *||||||||||||||||*||||||||||||||      ||||||*
//         // 130222330222330220330222330222332222330222330
//         // [INFO] 16:26:16 - Message: old score:  130  new score:  140
//         // Global alignment
//         // 2233    3022233022233302223  30222330222330222330222330222  330222  330222330222330222330222330222330222334
//         // ||||    |||||||||||*| ||||*  |||||||||||||||||||||||||||||  *|||||  |||||||||||*||||||||||||||*|||||||||||*
//         // 223322233022233022203 02220  30222330222330222330222330222  130222  330222330220330222330222332222330222330
//         // The original occurrence of the MEM has been shifted to the left by 6 positions, 
//         // reducing the gap in the right context, and moving in to the left context.

//         assert(ez.score >= score);

//         // Concatenate the CIGAR strings
        
//         // std::string cigar_s;
//         sam->cigar = "";
//         for(size_t i = 0; i < ez.n_cigar; ++i)
//           sam->cigar += std::to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];

//         // Compute the MD:Z field and the number of mismatches
//         sam->nm = write_MD_core((uint8_t*)ref,seq,ez.cigar,ez.n_cigar,tmp,0,sam->md);

//         sam->as = ez.score;
//         sam->pos = ref_pos + 1; // ref_pos is 1 based
//         sam->rname = idx[ref_pos];
//         sam->rlen = ref_len;
//         // write_sam(ez.score, score2, min_score, ref_pos, "human", read, strand, out, cigar_s, mdz_s, nm, "*", 0, 0);
//       }
//       else
//       {
//         // // Concatenate the CIGAR strings
//         // size_t n_cigar = ez_lc.n_cigar + ez_rc.n_cigar + 1;
//         // for(size_t i = 0; i < anchors.size()-1; ++i)
//         //   n_cigar += ez_cc[i].n_cigar + 1;

//         // uint32_t *cigar = (uint32_t*)calloc(n_cigar,sizeof(uint32_t));
//         // size_t i = 0;

//         // for(size_t j = 0; j < ez_lc.n_cigar; ++j)
//         //   cigar[i++] = ez_lc.cigar[ez_lc.n_cigar -j -1];

//         // // CIGARs between MEMs
//         // for(size_t j = 0; j < anchors.size(); ++j)
//         // {
//         //   // CIGAR for the MEM
//         //   size_t mem_len = mems[anchors[j].first].len;
//         //   if(i > 0 and ((cigar[i-1]& 0xf) == 0))
//         //   { // If the previous operation is also an M then merge the two operations
//         //     cigar[i-1] += (((uint32_t)mem_len) << 4);
//         //     --n_cigar;
//         //   }
//         //   else
//         //     cigar[i++] = (((uint32_t)mem_len) << 4);

//         //   // CIGAR of the gap between the j-th MEM and the j+1-th MEM
//         //   if(j < anchors.size() - 1)
//         //   {
//         //     if(ez_cc[j].n_cigar > 0)
//         //     {
//         //       // Check the first element
//         //       if((ez_cc[j].cigar[0]& 0xf) == 0)
//         //       { // If the next operation is also an M then merge the two operations
//         //         cigar[i-1] += ez_cc[j].cigar[0];
//         //         --n_cigar;
//         //       }
//         //       else
//         //         cigar[i++] = ez_cc[j].cigar[0];
//         //     } 
                      
//         //       // Copy all the other elements
//         //       for(size_t k = 1; k < ez_cc[j].n_cigar; ++k)
//         //         cigar[i++] = ez_cc[j].cigar[k];
//         //   }
//         // }

//         // // CIGAR of the right context
//         // if(ez_rc.n_cigar > 0)
//         // {
//         //   if((ez_rc.cigar[0]& 0xf) == 0)
//         //   { // If the next operation is also an M then merge the two operations
//         //     cigar[i-1] += ez_rc.cigar[0];
//         //     --n_cigar;
//         //   }
//         //   else
//         //     cigar[i++] = ez_rc.cigar[0];
//         // }

//         // for(size_t j = 1; j < ez_rc.n_cigar; ++j)
//         //   cigar[i++] = ez_rc.cigar[j];
        
//         // assert(i <= n_cigar);

//         // // std::string bfull = print_BLAST_like((uint8_t*)ref,seq,cigar,n_cigar);
//         // // std::cout << bfull;

//         // std::string cigar_s;
//         // for(size_t i = 0; i < n_cigar; ++i)
//         //   cigar_s += std::to_string(cigar[i] >> 4) + "MID"[cigar[i] & 0xf];


//         // // Compute the MD:Z field and thenumber of mismatches
//         // std::pair<std::string,size_t> md_nm = write_MD_core((uint8_t*)ref,seq,cigar,n_cigar,tmp,0);

//         // write_sam(score,score2,min_score,ref_pos,"human",read,strand,out,cigar_s,md_nm.first,md_nm.second);

//         // delete cigar;
//       }
//       delete tmp;
//     }
//     delete ref;
//     delete seq;

//     if (ez_lc.m_cigar > 0)
//       delete ez_lc.cigar;
//     if (ez_rc.m_cigar > 0)
//       delete ez_rc.cigar;
//     if (ez.m_cigar > 0)
//       delete ez.cigar;

//     return score;
//   }

//   // Readapted from https://github.com/lh3/minimap2/blob/c9874e2dc50e32bbff4ded01cf5ec0e9be0a53dd/format.c
//   // tmp is a string of length max(reference length, query length)
//   // it return the mdz string in the msdz parameter, and the number of mismatches as return value
//   static size_t write_MD_core(const uint8_t *tseq, const uint8_t *qseq, const uint32_t *cigar, const size_t n_cigar, char *tmp, int write_tag, std::string& mdz)
//   {
//     // std::string mdz;
//     int i, q_off, t_off, l_MD = 0, NM = 0;
//     if (write_tag) mdz += "MD:Z:"; //printf("MD:Z:");
//     for (i = q_off = t_off = 0; i < (int)n_cigar; ++i) {
//       int j, op = cigar[i]&0xf, len = cigar[i]>>4;
//       assert((op >= 0 && op <= 3) || op == 7 || op == 8);
//       if (op == 0 || op == 7 || op == 8) { // match
//         for (j = 0; j < len; ++j) {
//           if (qseq[q_off + j] != tseq[t_off + j]) {
//             mdz += std::to_string(l_MD) + "ACGTN"[tseq[t_off + j]];
//             // printf("%d%c", l_MD, "ACGTN"[tseq[t_off + j]]);
//             l_MD = 0; 
//             ++NM;
//           } else ++l_MD;
//         }
//         q_off += len, t_off += len;
//       } else if (op == 1) { // insertion to ref
//         q_off += len;
//         NM += len;
//       } else if (op == 2) { // deletion from ref
//         for (j = 0, tmp[len] = 0; j < len; ++j)
//           tmp[j] = "ACGTN"[tseq[t_off + j]];
//         mdz += std::to_string(l_MD) + "^" + std::string(tmp);
//         // printf("%d^%s", l_MD, tmp);
//         l_MD = 0;
//         t_off += len;
//         NM += len;
//       } else if (op == 3) { // reference skip
//         t_off += len;
//       }
//     }
//     if (l_MD > 0) mdz += std::to_string(l_MD);//printf("%d", l_MD);
//     // assert(t_off == r->re - r->rs && q_off == r->qe - r->qs);
//     return NM;
//     // return make_pair(mdz,NM);
//   }

//   // From https://github.com/lh3/ksw2/blob/master/cli.c
//   static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
//   {
//     int i, j;
//     a = a < 0? -a : a;
//     b = b > 0? -b : b;
//     for (i = 0; i < m - 1; ++i) {
//       for (j = 0; j < m - 1; ++j)
//         mat[i * m + j] = i == j? a : b;
//       mat[i * m + m - 1] = 0;
//     }
//     for (j = 0; j < m; ++j)
//       mat[(m - 1) * m + j] = 0;
//   }


//   // Adapted from https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/master/src/main.c
//   void write_sam(const int32_t score,
//                  const int32_t score2,
//                  const int32_t min_score,
//                  size_t ref_pos,
//                  const char *ref_seq_name,
//                  const kseq_t *read,
//                  int8_t strand, // 0: forward aligned ; 1: reverse complement aligned
//                  FILE *out,
//                  std::string &cigar,
//                  std::string &md,
//                  size_t mismatches,
//                  const char *r_next,  // Reference sequence name of the primary alignment of the NEXT read in the template.
//                  const size_t p_next, // 0-based Position of the primary alignment of the NEXT read in the template.
//                  const int32_t t_len  // TLEN: signed observed Template LENgth.
//   )
//   {
//     // Sam format output
//     fprintf(out, "%s\t", read->name.s);
//     if (score == 0)
//       fprintf(out, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
//     else
//     {
//       int32_t c, p;
//       // uint32_t mapq = -4.343 * log(1 - (double)abs(score - score2) / (double)score);
//       // mapq = (uint32_t)(mapq + 4.99);
//       // mapq = mapq < 254 ? mapq : 254;
//       uint32_t mapq = compute_mapq(score,score2,min_score,read->seq.l); 
//       if (strand)
//         fprintf(out, "16\t");
//       else
//         fprintf(out, "0\t");
//       // TODO: Find the correct reference name.
//       fprintf(out, "%s\t%d\t%d\t", ref_seq_name, ref_pos + 1, mapq);
//       fprintf(out, "%s\t", cigar.c_str());
//       fprintf(out, "%s\t%d\t%d\t", r_next, p_next + 1, t_len);
//       // fprintf(out, "\t*\t0\t0\t");
//       fprintf(out, "%s", read->seq.s);
//       fprintf(out, "\t");
//       if (read->qual.s && strand)
//       {
//         for (p = read->qual.l - 1; p >= 0; --p)
//           fprintf(out, "%c", read->qual.s[p]);
//       }
//       else if (read->qual.s)
//         fprintf(out, "%s", read->qual.s);
//       else
//         fprintf(out, "*");
//       fprintf(out, "\tAS:i:%d", score);
//       fprintf(out, "\tNM:i:%d", mismatches);
//       if (score2 > 0)
//         fprintf(out, "\tZS:i:%d", score2);
//       fprintf(out, "\tMD:Z:%s\n", md.c_str());
//     }
//   }

//   /*!
//     Compute the mapping quality of the alignment
//     Inspired from https://github.com/BenLangmead/bowtie2/blob/4512b199768e562e8627ffdfd9253affc96f6fc6/unique.h
//   */
//   size_t compute_mapq(
//     const int32_t score,      // Best alignment score
//     const int32_t score2,     // Second best alignemt score
//     const int32_t min_score,  // Minimum alignemt score
//     const size_t read_l       // Read length
//   )
//   {
//     int32_t max_score = read_l * smatch;
//     int32_t best = max_score - score;
//     size_t best_bin = (size_t)((double)best * (10.0 / (double)(max_score - min_score)) + 0.5);
//     if(score2 >= min_score)
//     {
//       int32_t diff = score - score2;
//       size_t diff_bin = (size_t)((double)diff * (10.0 / (double)(max_score - min_score)) + 0.5);
//       if(best == max_score)
//         return unp_sec_perf[best_bin];
//       else
//         return unp_sec[diff_bin][best_bin];
//     }
//     else
//     {
//       if(best == max_score)
//         return unp_nosec_perf;
//       else
//         return unp_nosec[best_bin];
//     }
//   }

//   static std::string print_BLAST_like(const uint8_t *tseq, const uint8_t *qseq, const uint32_t *cigar, const size_t n_cigar)
//   {
//     std::string target_o;
//     std::string bars_o;
//     std::string seq_o;


//     int i, q_off, t_off, l_MD = 0;
//     for (i = q_off = t_off = 0; i < (int)n_cigar; ++i) {
//       int j, op = cigar[i]&0xf, len = cigar[i]>>4;
//       assert((op >= 0 && op <= 3) || op == 7 || op == 8);
//       if (op == 0 || op == 7 || op == 8) { // match
//         for (j = 0; j < len; ++j) {
//           if (qseq[q_off + j] != tseq[t_off + j]) {
//             bars_o +="*";
//           } else {
//             bars_o +="|";
//           }
//           target_o += std::to_string(tseq[t_off + j]);
//           seq_o += std::to_string(qseq[q_off + j]);
//         }
//         q_off += len, t_off += len;
//       } else if (op == 1) { // insertion to ref
//         for (j = 0; j < len; ++j) {
//           target_o += " ";
//           bars_o += " ";
//           seq_o += std::to_string(qseq[q_off + j]);
//         }
//         q_off += len;
//       } else if (op == 2) { // deletion from ref
//         for (j = 0; j < len; ++j) {
//           seq_o += " ";
//           bars_o += " ";
//           target_o += std::to_string(tseq[t_off + j]);
//         }
//         t_off += len;
//       } else if (op == 3) { // reference skip
//         for (j = 0; j < len; ++j) {
//           seq_o += " ";
//           bars_o += " ";
//           target_o += std::to_string(tseq[t_off + j]);
//         }
//         t_off += len;
//       }
//     }
//     return target_o + "\n" + bars_o + "\n" + seq_o + "\n";
//   }

//   std::string to_sam()
//   {
//       std::string res = "@HD VN:1.6 SO:unknown\n";
//       res += idx.to_sam();
//       res += "@PG ID:moni PN:moni VN:0.1.0\n";
//       return res; 
//   }

// protected:
//   ms_t ms;
//   slp_t ra;
//   seqidx idx;
//   // SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd> ra;

//   size_t min_len = 0;
//   size_t aligned_reads = 0;
//   size_t n = 0;
//   size_t top_k = 1; // report the top_k alignments

//   unsigned char seq_nt4_table[256] = {
//       0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

//   int8_t smatch = 2;      // Match score default
//   int8_t smismatch = 4;   // Mismatch score default
//   int8_t gapo = 4;        // Gap open penalty
//   int8_t gapo2 = 13;      // Gap open penalty
//   int8_t gape = 2;        // Gap extension penalty
//   int8_t gape2 = 1;       // Gap extension penalty
//   int end_bonus = 400;    // Bonus to add at the extension score to declare the alignment
  
//   int w = -1;             // Band width
//   int zdrop = -1;

// 	void *km = 0;           // Kalloc

//   // int8_t max_rseq = 0;

//   int m = 5;
//   int8_t mat[25];
//   // int minsc = 0, xtra = KSW_XSTART;
//   // uint8_t *rseq = 0;

//   bool forward_only;


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
// };

// ////////////////////////////////////////////////////////////////////////////////
// /// Merge SAMs
// ////////////////////////////////////////////////////////////////////////////////

// // Merges te file in filename in the file pointed by fp
// void append_file(const std::string filename, FILE *fp)
// {
//   const size_t buff_size = 16384;

//   uint8_t buff[buff_size];
//   size_t size = 0;

//   struct stat filestat;
//   FILE *fd;

//   if ((fd = fopen(filename.c_str(), "r")) == nullptr)
//     error("open() file " + std::string(filename) + " failed");

//   // int fn = fileno(fd);
//   // if (fstat(fn, &filestat) < 0)
//   //     error("stat() file " + std::string(filename) + " failed");

//   // size_t length = filestat.st_size;
//   size_t length = 0;

//   while ((length = fread(buff, sizeof(uint8_t), buff_size, fd)) == buff_size)
//     if ((fwrite(buff, sizeof(uint8_t), buff_size, fp)) != buff_size)
//       error("fwrite() file " + std::string(filename) + " failed");

//   assert(length < buff_size);
//   if (length > 0)
//     if ((fwrite(buff, sizeof(uint8_t), length, fp)) != length)
//       error("fwrite() file " + std::string(filename) + " failed");

//   fclose(fd);
// }

// ////////////////////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////////////
// /// Parallel computation
// ////////////////////////////////////////////////////////////////////////////////

// pthread_mutex_t mutex_reads_dispatcher;

// inline static size_t mt_kbseq_read(kbseq_t *b, kseq_t *seq, const size_t n)
// {
//   size_t l = 0;
//   // Update the number of active threads
//   xpthread_mutex_lock(&mutex_reads_dispatcher, __LINE__, __FILE__);
//   {
//     l = kbseq_read(b, seq, n);
//   }
//   xpthread_mutex_unlock(&mutex_reads_dispatcher, __LINE__, __FILE__);
//   return l;
// }

// inline static size_t mt_kpbseq_read(kpbseq_t *b, kseq_t *mate1, kseq_t *mate2, const size_t n)
// {
//   size_t l = 0;
//   // Update the number of active threads
//   xpthread_mutex_lock(&mutex_reads_dispatcher, __LINE__, __FILE__);
//   {
//     l = kpbseq_read(b, mate1, mate2, n);
//   }
//   xpthread_mutex_unlock(&mutex_reads_dispatcher, __LINE__, __FILE__);
//   return l;
// }

// template <typename aligner_t>
// struct mt_param_t
// {
//   // Parameters
//   aligner_t *aligner;
//   std::string pattern_filename;
//   std::string sam_filename;
//   size_t start;
//   size_t end;
//   size_t wk_id;
//   size_t b_size;
//   kseq_t *seq;
//   kseq_t *mate1;
//   kseq_t *mate2;
//   // Return values
//   size_t n_reads;
//   size_t n_aligned_reads;
// };

// // template <typename aligner_t>
// // void *mt_align_worker(void *param)
// // {
// //   mt_param_t<aligner_t> *p = (mt_param_t<aligner_t>*) param;
// //   size_t n_reads = 0;
// //   size_t n_aligned_reads = 0;

// //   FILE *sam_fd;
// //   gzFile fp;

// //   if ((sam_fd = fopen(p->sam_filename.c_str(), "w")) == nullptr)
// //     error("open() file " + p->sam_filename + " failed");

// //   if ((fp = gzopen(p->pattern_filename.c_str(), "r")) == Z_NULL)
// //     error("open() file " + p->pattern_filename + " failed");

// //   gzseek(fp, p->start, SEEK_SET);

// //   kseq_t rev;
// //   int l;

// //   kseq_t *seq = kseq_init(fp);
// //   while ((ks_tell(seq) < p->end) && ((l = kseq_read(seq)) >= 0))
// //   {

// //     bool fwd_align = p->aligner->align(seq, sam_fd, 0);

// //     //copy seq
// //     copy_kseq_t(&rev, seq);

// //     for (size_t i = 0; i < seq->seq.l; ++i)
// //       rev.seq.s[i] = complement(seq->seq.s[seq->seq.l - i - 1]);

// //     if (rev.seq.m > rev.seq.l)
// //       rev.seq.s[rev.seq.l] = 0;

// //     bool rev_align = p->aligner->align(&rev, sam_fd, 1);

// //     if (fwd_align or rev_align)
// //       n_aligned_reads++;
// //     n_reads++;

// //     free(rev.name.s);
// //     free(rev.comment.s);
// //     free(rev.seq.s);
// //     free(rev.qual.s);
// //   }

// //   verbose("Number of aligned reads block ", p->wk_id, " : ", n_aligned_reads, "/", n_reads);
// //   p->n_reads = n_reads;
// //   p->n_aligned_reads = n_aligned_reads;
// //   kseq_destroy(seq);
// //   gzclose(fp);
// //   fclose(sam_fd);

// //   return NULL;
// // }
// // template <typename aligner_t>
// // size_t mt_align(aligner_t *aligner, std::string pattern_filename, std::string sam_filename, size_t n_threads)
// // {
// //   pthread_t t[n_threads] = {0};
// //   mt_param_t<aligner_t> params[n_threads];
// //   std::vector<size_t> starts = split_fastq(pattern_filename, n_threads);
// //   for(size_t i = 0; i < n_threads; ++i)
// //   {
// //     params[i].aligner = aligner;
// //     params[i].pattern_filename = pattern_filename;
// //     params[i].sam_filename = sam_filename + "_" + std::to_string(i) + ".sam";
// //     params[i].start = starts[i];
// //     params[i].end = starts[i+1];
// //     params[i].wk_id = i;
// //     xpthread_create(&t[i], NULL, &mt_align_worker<aligner_t>, &params[i], __LINE__, __FILE__);
// //   }

// //   size_t tot_reads = 0;
// //   size_t tot_aligned_reads = 0;

// //   for(size_t i = 0; i < n_threads; ++i)
// //   {
// //     xpthread_join(t[i],NULL,__LINE__,__FILE__);
// //   }

// //   sleep(5);
// //   for(size_t i = 0; i < n_threads; ++i)
// //   {
// //     tot_reads += params[i].n_reads;
// //     tot_aligned_reads += params[i].n_aligned_reads;
// //   }

// //   verbose("Number of aligned reads: ", tot_aligned_reads, "/", tot_reads);
// //   return tot_aligned_reads;
// // }

// template <typename aligner_t>
// void *mt_align_worker(void *param)
// {
//   mt_param_t<aligner_t> *p = (mt_param_t<aligner_t> *)param;
//   size_t n_reads = 0;
//   size_t n_aligned_reads = 0;

//   FILE *sam_fd;
//   gzFile fp;

//   if ((sam_fd = fopen(p->sam_filename.c_str(), "w")) == nullptr)
//     error("open() file " + p->sam_filename + " failed");

//   size_t b_size = p->b_size;
//   if (p->seq != nullptr)
//   {
//     kseq_t *seq = p->seq;
//     kbseq_t *b = kbseq_init();
//     int l = 0;
//     while ((l = mt_kbseq_read(b, seq, b_size)) > 0)
//     {
//       for (size_t i = 0; i < l; ++i)
//       {
//         if (p->aligner->align(&b->buf[i], sam_fd, 0))
//           n_aligned_reads++;        
//         n_reads++;
//       }
//       // std::cout << "Block p " << p->wk_id << " end!" << std::endl;
//     }
//     kbseq_destroy(b);
//   }
//   else
//   {
//     kseq_t *mate1 = p->mate1;
//     kseq_t *mate2 = p->mate2;
//     kpbseq_t *b = kpbseq_init();
//     int l = 0;
//     while ((l = mt_kpbseq_read(b, mate1, mate2, b_size)) > 0)
//     {
//       for (size_t i = 0; i < l; ++i)
//       {
//         if (p->aligner->align(&b->mate1->buf[i], &b->mate2->buf[i],sam_fd))
//           n_aligned_reads++;
//         n_reads++;
//       }
//       // std::cout << "Block p " << p->wk_id << " end!" << std::endl;
//     }
//     kpbseq_destroy(b);
//   }

//   verbose("Number of aligned reads block ", p->wk_id, " : ", n_aligned_reads, "/", n_reads);
//   p->n_reads = n_reads;
//   p->n_aligned_reads = n_aligned_reads;
//   fclose(sam_fd);

//   return NULL;
// }

// template <typename aligner_t>
// size_t mt_align(aligner_t *aligner, std::string pattern_filename, std::string sam_filename, size_t n_threads, size_t b_size, std::string mate2_filename = "")
// {
//   xpthread_mutex_init(&mutex_reads_dispatcher, NULL, __LINE__, __FILE__);
//   kseq_t *seq = nullptr;
//   kseq_t *mate1 = nullptr;
//   kseq_t *mate2 = nullptr;

//   gzFile fp_mate2 = nullptr;
//   gzFile fp = gzopen(pattern_filename.c_str(), "r");
//   if (mate2_filename != "")
//   {
//     fp_mate2 = gzopen(mate2_filename.c_str(), "r");
//     mate1 = kseq_init(fp);
//     mate2 = kseq_init(fp_mate2);
//   }
//   else
//   {
//     seq = kseq_init(fp);
//   }

//   pthread_t t[n_threads] = {0};
//   mt_param_t<aligner_t> params[n_threads];
//   for (size_t i = 0; i < n_threads; ++i)
//   {
//     // Create a new thread
//     params[i].aligner = aligner;
//     params[i].pattern_filename = pattern_filename;
//     params[i].sam_filename = sam_filename + "_" + std::to_string(i) + ".sam";
//     params[i].b_size = b_size;
//     params[i].seq = seq;
//     params[i].mate1 = mate1;
//     params[i].mate2 = mate2;
//     params[i].wk_id = i;
//     xpthread_create(&t[i], NULL, &mt_align_worker<aligner_t>, &params[i], __LINE__, __FILE__);
//   }

//   size_t tot_reads = 0;
//   size_t tot_aligned_reads = 0;

//   for (size_t i = 0; i < n_threads; ++i)
//   {
//     xpthread_join(t[i], NULL, __LINE__, __FILE__);
//   }

//   if (fp_mate2 != nullptr)
//   {
//     kseq_destroy(mate1);
//     kseq_destroy(mate2);
//     gzclose(fp_mate2);
//   }
//   else
//   {
//     kseq_destroy(seq);
//   }
//   gzclose(fp);

//   // sleep(5);
//   verbose("Merging temporary SAM files");

//   FILE *fd;

//   if ((fd = fopen(std::string(sam_filename + ".sam").c_str(), "w")) == nullptr)
//     error("open() file " + std::string(sam_filename + ".sam") + " failed");

//   fprintf(fd, "%s", aligner->to_sam().c_str());

//   for (size_t i = 0; i < n_threads; ++i)
//   {
//     tot_reads += params[i].n_reads;
//     tot_aligned_reads += params[i].n_aligned_reads;

//     append_file(params[i].sam_filename, fd);
//     if (std::remove(params[i].sam_filename.c_str()) != 0)
//       error("remove() file " + params[i].sam_filename + " failed");
//   }

//   xpthread_mutex_destroy(&mutex_reads_dispatcher, __LINE__, __FILE__);

//   verbose("Number of aligned reads: ", tot_aligned_reads, "/", tot_reads);
//   return tot_aligned_reads;
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// Single Thread
// ////////////////////////////////////////////////////////////////////////////////
// template <typename aligner_t>
// size_t st_align(aligner_t *aligner, std::string pattern_filename, std::string sam_filename, size_t b_size, std::string mate2_filename = "")
// {
//   size_t n_reads = 0;
//   size_t n_aligned_reads = 0;

//   kseq_t *seq = nullptr;
//   kseq_t *mate1 = nullptr;
//   kseq_t *mate2 = nullptr;

//   gzFile fp_mate2 = nullptr;
//   gzFile fp = gzopen(pattern_filename.c_str(), "r");
//   if (mate2_filename != "")
//   {
//     fp_mate2 = gzopen(mate2_filename.c_str(), "r");
//     mate1 = kseq_init(fp);
//     mate2 = kseq_init(fp_mate2);
//   }
//   else
//   {
//     seq = kseq_init(fp);
//   }

//   int l;
//   FILE *sam_fd;

//   sam_filename += ".sam";

//   if ((sam_fd = fopen(sam_filename.c_str(), "w")) == nullptr)
//     error("open() file " + sam_filename + " failed");

//   fprintf(sam_fd, "%s", aligner->to_sam().c_str());

//   if (seq != nullptr)
//   {
//     kbseq_t *b = kbseq_init();
//     int l = 0;
//     while ((l = kbseq_read(b, seq, b_size)) > 0)
//     {
//       for (size_t i = 0; i < l; ++i)
//       {
//         if (aligner->align(&b->buf[i], sam_fd, 0))
//           n_aligned_reads++;
//         n_reads++;
//       }
//     }
//     kbseq_destroy(b);
//   } 
//   else
//   {
//     kpbseq_t *b = kpbseq_init();
//     int l = 0;
//     while ((l = kpbseq_read(b, mate1, mate2, b_size)) > 0)
//     {
//       for (size_t i = 0; i < l; ++i)
//       {
//         if (aligner->align(&b->mate1->buf[i], &b->mate2->buf[i], sam_fd))
//           n_aligned_reads++;
//         n_reads++;
//       }
//       // std::cout << "Block p " << p->wk_id << " end!" << std::endl;
//     }
//     kpbseq_destroy(b);
//   }
//   verbose("Number of aligned reads: ", n_aligned_reads, "/", n_reads);
//   kseq_destroy(seq);
//   gzclose(fp);
//   fclose(sam_fd);

//   sleep(5);

//   return n_aligned_reads;
// }

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
  std::string mate1 = "";    // path to file with #1 mates paired with mate2.
  std::string mate2 = "";    // path to file with #2 mates paired with mate1.
  std::string output   = ""; // output file prefix
  size_t b = 32768;       // size of the batch of read to be processed
  size_t l = 25; // minumum MEM length
  size_t th = 1; // number of threads
  bool is_fasta = false; // read a fasta file
  bool shaped_slp = false; // use shaped slp
  size_t ext_len = 100;      // Extension length
  // size_t top_k = 1;       // Report the top_k alignments

  // ksw2 parameters
  int8_t smatch = 2;      // Match score default
  int8_t smismatch = 4;   // Mismatch score default
  int8_t gapo = 4;        // Gap open penalty
  int8_t gapo2 = 13;      // Gap open penalty
  int8_t gape = 2;        // Gap extension penalty
  int8_t gape2 = 1;       // Gap extension penalty
  // int end_bonus = 400;    // Bonus to add at the extension score to declare the alignment

  // int w = -1;             // Band width
  // int zdrop = -1;         // Zdrop enable
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string usage("usage: " + std::string(argv[0]) + " infile [-p patterns] [-1 mate1] [-2 mate2] [-o output] [-t threads] [-b batch] [-l len] [-q shaped_slp]  [-L ext_l] [-A smatch] [-B smismatc] [-O gapo] [-E gape]\n\n" +
                    "Align the reads in the pattern against the reference index in infile.\n" +
                    "   pattens: [string]  - path to patterns file.\n" +
                    "     mate1: [string]  - path to file with #1 mates paired with mate2.\n" +
                    "     mate2: [string]  - path to file with #2 mates paired with mate1.\n" +
                    "    output: [string]  - output file prefix.\n" +
                    "   threads: [integer] - number of threads (def. " + std::to_string(arg.th) + ")\n" +
                    "     batch: [integer] - number of batches per therad pool (def. " + std::to_string(arg.b) + ")\n" + 
                    "       len: [integer] - minimum MEM lengt (def. " + std::to_string(arg.l) + ")\n" +
                    "shaped_slp: [boolean] - use shaped slp. (def. " + std::to_string(arg.shaped_slp) + ")\n" +
                    "     ext_l: [integer] - length of reference substring for extension (def. " + std::to_string(arg.ext_len) + ")\n" +
                    "    smatch: [integer] - match score value (def. " + std::to_string(arg.smatch) + ")\n" +
                    " smismatch: [integer] - mismatch penalty value (def. " + std::to_string(arg.smismatch) + ")\n" +
                    "      gapo: [integer] - gap open penalty value (def. " + std::to_string(arg.gapo) + "," + std::to_string(arg.gapo2) + ")\n" +
                    "      gape: [integer] - gap extension penalty value (def. " + std::to_string(arg.gape) + "," + std::to_string(arg.gape2) + ")\n");

  std::string sarg;
  char* s;
  while ((c = getopt(argc, argv, "ql:hp:o:t:1:2:b:A:B:O:E:L:")) != -1)
  {
    switch (c)
    {
    case 'p':
      arg.patterns.assign(optarg);
      break;
    case 'o':
      arg.output.assign(optarg);
      break;
    case '1':
      arg.mate1.assign(optarg);
      break;
    case '2':
      arg.mate2.assign(optarg);
      break;
    case 'l':
      sarg.assign(optarg);
      arg.l = stoi(sarg);
      break;
    case 'b':
      sarg.assign(optarg);
      arg.b = stoi(sarg);
      break;
    case 't':
      sarg.assign(optarg);
      arg.th = stoi(sarg);
      break;
    case 'L':
      sarg.assign(optarg);
      arg.ext_len = stoi(sarg);
      break;
    case 'A':
      sarg.assign(optarg);
      arg.smatch = stoi(sarg);
      break;
    case 'B':
      sarg.assign(optarg);
      arg.smismatch = stoi(sarg);
      break;
    case 'O':
      arg.gapo = arg.gapo2 = strtol(optarg, &s, 10);
      if (*s == ',') arg.gapo2 = strtol(s+1, &s, 10);
      break;
    case 'E':
      arg.gape = arg.gape2 = strtol(optarg, &s, 10);
      if (*s == ',') arg.gape2 = strtol(s+1, &s, 10);
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
typename aligner_t::config_t configurer(Args &args){
  typename aligner_t::config_t config;
  
  config.min_len    = args.l;           // Minimum MEM length
  config.ext_len    = args.ext_len;     // Extension length

  // ksw2 parameters
  config.smatch     = args.smatch;      // Match score default
  config.smismatch  = args.smismatch;   // Mismatch score default
  config.gapo       = args.gapo;        // Gap open penalty
  config.gapo2      = args.gapo2;       // Gap open penalty
  config.gape       = args.gape;        // Gap extension penalty
  config.gape2      = args.gape2;       // Gap extension penalty

  return config;
}

template<typename aligner_t>
void dispatcher(Args &args){

  verbose("Construction of the aligner");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  aligner_t aligner(args.filename, configurer<aligner_t>(args));

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();

  std::string base_name = basename(args.filename.data());


  if(args.patterns != "")
  {
    std::string sam_filename = args.patterns + "_" + base_name + "_" + std::to_string(args.l) + ".sam";
    if(args.output != "")
      sam_filename = args.output;

    verbose("Output file: ", sam_filename);

    if (is_gzipped(args.patterns))
    {
      verbose("The input is gzipped - forcing single thread alignment.");
      args.th = 1;
    }
    if (args.th == 1)
      st_align<aligner_t>(&aligner, args.patterns, sam_filename, args.b);
    else
      mt_align<aligner_t>(&aligner, args.patterns, sam_filename, args.th, args.b);
  }
  else
  {
    std::string sam_filename = args.mate1 + "_" + base_name + "_" + std::to_string(args.l) + ".sam";
    if(args.output != "")
      sam_filename = args.output;

    verbose("Output file: ", sam_filename);

    if (is_gzipped(args.mate1) or is_gzipped(args.mate2))
    {
      verbose("The input is gzipped - forcing single thread alignment.");
      args.th = 1;
    }
    if (args.th == 1)
      st_align<aligner_t>(&aligner, args.mate1, sam_filename, args.b, args.mate2 );
    else
      mt_align<aligner_t>(&aligner, args.mate1, sam_filename, args.th, args.b, args.mate2);
  }

  // TODO: Merge the SAM files.

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

}

int main(int argc, char *const argv[])
{

  Args args;
  parseArgs(argc, argv, args);

  if(args.shaped_slp){
    dispatcher<aligner<shaped_slp_t, ms_pointers<>>>(args);
  }else{
    dispatcher<aligner<plain_slp_t, ms_pointers<>>>(args);
  }

  MTIME_REPORT_ALL;

  return 0;
}