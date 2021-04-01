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

#include <ms_pointers.hpp>

#include <malloc_count.h>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>

#include <ksw.h>
#include <ssw.h>

#include <omp.h>

#include <libgen.h>

// KSEQ_INIT(gzFile, gzread);

class aligner_t
{
public:
  using SelSd = SelectSdvec<>;
  using DagcSd = DirectAccessibleGammaCode<SelSd>;

  aligner_t(std::string filename, 
            size_t min_len_ = 50, 
            bool forward_only_ = true): 
                min_len(min_len_), 
                forward_only(forward_only_)
  {
    verbose("Loading the matching statistics index");
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    std::string filename_ms = filename  + ms.get_file_extension();

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

    if (minsc > 0xffff)
      minsc = 0xffff;
    xtra |= KSW_XSUBO | minsc;
    // initialize scoring matrix
    for (i = k = 0; i < 4; ++i)
    {
      for (j = 0; j < 4; ++j)
        mat[k++] = i == j ? sa : -sb;
      mat[k++] = 0; // ambiguous base
    }
    for (j = 0; j < 5; ++j)
      mat[k++] = 0;

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Local aligner initialization complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    verbose("Minimum MEM length: ", min_len);
  }


  bool align(kseq_t *read, FILE* out, uint8_t strand)
  {
    size_t mem_pos = 0;
    size_t mem_len = 0;
    size_t mem_idx = 0;

    bool aligned = false;

    auto pointers = ms.query(read->seq.s, read->seq.l);
    std::vector<size_t> lengths(pointers.size());
    size_t l = 0;
    for (size_t i = 0; i < pointers.size(); ++i)
    {
      size_t pos = pointers[i];
      while ((i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
        ++l;

      lengths[i] = l;
      l = (l == 0 ? 0 : (l - 1));

      // Update MEM
      if (lengths[i] > mem_len)
      {
        mem_len = lengths[i];
        mem_pos = pointers[i];
        mem_idx = i;
      }
    }

    // Align the read
    if (mem_len >= min_len)
    {
      char *str = (char *)malloc(400);

      int32_t maskLen = read->seq.l / 2;
      maskLen = maskLen < 15 ? 15 : maskLen;

      // Extract the context from the reference
      size_t left_occ = (mem_pos > 100 ? mem_pos - 100 : 0);
      size_t len = mem_len + 100 + (mem_pos > 100 ? 100 : 100 - mem_pos);
      ra.expandSubstr(left_occ, len, str);

      size_t min_score = 20 + 8 * log(read->seq.l);

      uint8_t* seq = (uint8_t*)malloc(read->seq.l);
      // Convert A,C,G,T,N into 0,1,2,3,4
      for (i = 0; i < (int)read->seq.l; ++i)
        seq[i] = seq_nt4_table[(int)read->seq.s[i]];
      // for (i = 0; i < (int)read->seq.l; ++i)
      //   read->seq.s[i] = seq_nt4_table[(int)read->seq.s[i]];

      for (i = 0; i < (int)len; ++i)
        str[i] = seq_nt4_table[(int)str[i]];

      int score;


      kswq_t *q = 0;
      kswr_t r;
      
      r = ksw_align(read->seq.l, (uint8_t *)seq, len, (uint8_t *)str, 5, mat, gapo, gape, xtra, &q);
      // score = ksw_global(read->seq.l, (uint8_t *)read->seq.s, len, (uint8_t *)str, 5, mat, gapo, gape, w, &n_cigar, &cigar);

      int n_cigar;
      uint32_t * cigar;

      size_t new_seq_len = r.qe - r.qb;
      size_t new_ref_len = r.te - r.tb;
      uint8_t *new_seq = (uint8_t *)(seq + r.qb);
      // uint8_t *new_seq = (uint8_t *)(read->seq.s + r.qb);
      uint8_t *new_ref = (uint8_t *)(str + r.tb);

      score = ksw_global(new_seq_len, (uint8_t *) new_seq, new_ref_len, new_ref, 5, mat, gapo, gape, w, &n_cigar, &cigar);

      std::string cig;

      // for(size_t i = 0; i < n_cigar; ++i)
      // {
      //   // for (i = 0; i < ez->n_cigar; ++i)
      //   //   printf("%d%c", ez->cigar[i] >> 4, "MID"[ez->cigar[i] & 0xf]);
      //   cig += std::to_string(cigar[i] >> 4) + "MID"[cigar[i] & 0xf];
      // }

      size_t mismatch = mark_mismatch(r.tb, r.qb, r.qe, (int8_t *)str, (int8_t *)seq, read->seq.l, &cigar, &n_cigar);
      for (c = 0; c < (n_cigar); ++c)
      {
        char letter = cigar_int_to_op(cigar[c]);
        uint32_t length = cigar_int_to_len(cigar[c]);
        // fprintf(out, "%lu%c", (unsigned long)length, letter);
        cig += std::to_string((unsigned long)length) + letter;
      }

      // if(r.score > 0)
      //   printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n", "human", r.tb, r.te + 1, read->name.s, r.qb, r.qe + 1, r.score, r.score2, r.te2);
      //   // std::cout << "\rCurrent score... "<< r.score << std::flush;

      // // Declares a default Aligner
      // StripedSmithWaterman::Aligner aligner;
      // // Declares a default filter
      // StripedSmithWaterman::Filter filter;
      // // StripedSmithWaterman::Filter filter(true, true, min_score, 32767);
      // // Declares an alignment that stores the result
      // StripedSmithWaterman::Alignment alignment;
      // // Aligns the query to the ref
      // aligner.Align(read->seq.s, str, len, filter, &alignment, maskLen);

      // // Update alignment method
      r.tb += left_occ;
      r.te += left_occ;
      r.te2 += left_occ;

      if(r.score >= min_score)
      {
        ssw_write_sam(r,"human",read,strand,out,cig,mismatch);
        aligned = true;
      }

      // aligned_reads++;
      free(cigar);
      free(q);
      delete str;
      delete seq;
    }
    return aligned;
  }

  size_t get_aligned_reads()
  {
    return aligned_reads;
  }

  // Adapted from SSW
  static void ssw_write_sam(kswr_t &a,
                            const char *ref_seq_name,
                            const kseq_t *read,
                            int8_t strand,
                            FILE *out,
                            std::string cigar,
                            size_t mismatches) // 0: forward aligned ; 1: reverse complement aligned
  {
    // Sam format output
    fprintf(out, "%s\t", read->name.s);
    if (a.score == 0)
      fprintf(out, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
    else
    {
      int32_t c, p;
      uint32_t mapq = -4.343 * log(1 - (double)abs(a.score - a.score2) / (double)a.score);
      mapq = (uint32_t)(mapq + 4.99);
      mapq = mapq < 254 ? mapq : 254;
      if (strand)
        fprintf(out, "16\t");
      else
        fprintf(out, "0\t");
      // TODO: Find the correct reference name.
      fprintf(out, "%s\t%d\t%d\t", ref_seq_name, a.tb + 1, mapq);
      // size_t mismatch = mark_mismatch(a.tb, a.qb, a.qe, (int8_t*)ref, (int8_t*)read_, read->seq.l, cigar, cigarLen);
      // for (c = 0; c < (*cigarLen); ++c)
      // {
      //   char letter = cigar_int_to_op((*cigar)[c]);
      //   uint32_t length = cigar_int_to_len((*cigar)[c]);
      //   fprintf(out, "%lu%c", (unsigned long)length, letter);
      // }
      // fprintf(out, "\t*\t");
      // fprintf(out, "%s", a.cigar_string.c_str());
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
      fprintf(out, "\tAS:i:%d", a.score);
      fprintf(out, "\tNM:i:%d\t", mismatches);
      // fprintf(out, "\tNM:i:%d\t", a.mismatches);
      if (a.score2 > 0)
        fprintf(out, "ZS:i:%d\n", a.score2);
      else
        fprintf(out, "\n");
    }
  }


protected:
  ms_pointers<> ms;
  SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd> ra;

  size_t min_len = 0;
  size_t aligned_reads = 0;
  size_t n = 0;

  unsigned char seq_nt4_table[256] = {
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
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

  int c, sa = 2, sb = 2, i, j, k, max_rseq = 0;
  int w = 4000;
  int8_t mat[25];
  int gapo = 5, gape = 2, minsc = 0, xtra = KSW_XSTART;
  uint8_t *rseq = 0;

  bool forward_only;
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