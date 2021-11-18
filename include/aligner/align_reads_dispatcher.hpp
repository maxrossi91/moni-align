/* align_reads_dispatcher - Dispatches the reads in single and multithread.
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
   \file align_reads_dispatcher.cpp
   \brief align_reads_dispatcher.cpp Dispatches the reads in single and multithread.
   \author Massimiliano Rossi
   \date 29/04/2021
*/

#ifndef _READS_DISPATCHER_HH
#define _READS_DISPATCHER_HH


extern "C"{
#include <xerrors.h>
}

#include <common.hpp>
#include <kseq.h>
#include <zlib.h>

////////////////////////////////////////////////////////////////////////////////
/// Merge SAMs
////////////////////////////////////////////////////////////////////////////////

// Merges te file in filename in the file pointed by fp
void append_file(const std::string filename, FILE *fp)
{
  const size_t buff_size = 16384;

  uint8_t buff[buff_size];
  size_t size = 0;

  struct stat filestat;
  FILE *fd;

  if ((fd = fopen(filename.c_str(), "r")) == nullptr)
    error("open() file " + std::string(filename) + " failed");

  size_t length = 0;

  while ((length = fread(buff, sizeof(uint8_t), buff_size, fd)) == buff_size)
    if ((fwrite(buff, sizeof(uint8_t), buff_size, fp)) != buff_size)
      error("fwrite() file " + std::string(filename) + " failed");

  assert(length < buff_size);
  if (length > 0)
    if ((fwrite(buff, sizeof(uint8_t), length, fp)) != length)
      error("fwrite() file " + std::string(filename) + " failed");

  fclose(fd);
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Parallel computation
////////////////////////////////////////////////////////////////////////////////

pthread_mutex_t mutex_reads_dispatcher;

inline static size_t mt_kbseq_read(kbseq_t *b, kseq_t *seq, const size_t n)
{
  size_t l = 0;
  // Update the number of active threads
  xpthread_mutex_lock(&mutex_reads_dispatcher, __LINE__, __FILE__);
  {
    l = kbseq_read(b, seq, n);
  }
  xpthread_mutex_unlock(&mutex_reads_dispatcher, __LINE__, __FILE__);
  return l;
}

inline static size_t mt_kpbseq_read(kpbseq_t *b, kseq_t *mate1, kseq_t *mate2, const size_t n)
{
  size_t l = 0;
  // Update the number of active threads
  xpthread_mutex_lock(&mutex_reads_dispatcher, __LINE__, __FILE__);
  {
    l = kpbseq_read(b, mate1, mate2, n);
  }
  xpthread_mutex_unlock(&mutex_reads_dispatcher, __LINE__, __FILE__);
  return l;
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
  size_t b_size;
  kseq_t *seq;
  kseq_t *mate1;
  kseq_t *mate2;
  // Return values
  size_t n_reads;
  size_t n_aligned_reads;
};

// template <typename aligner_t>
// void *mt_align_worker(void *param)
// {
//   mt_param_t<aligner_t> *p = (mt_param_t<aligner_t>*) param;
//   size_t n_reads = 0;
//   size_t n_aligned_reads = 0;

//   FILE *sam_fd;
//   gzFile fp;

//   if ((sam_fd = fopen(p->sam_filename.c_str(), "w")) == nullptr)
//     error("open() file " + p->sam_filename + " failed");

//   if ((fp = gzopen(p->pattern_filename.c_str(), "r")) == Z_NULL)
//     error("open() file " + p->pattern_filename + " failed");

//   gzseek(fp, p->start, SEEK_SET);

//   kseq_t rev;
//   int l;

//   kseq_t *seq = kseq_init(fp);
//   while ((ks_tell(seq) < p->end) && ((l = kseq_read(seq)) >= 0))
//   {

//     bool fwd_align = p->aligner->align(seq, sam_fd, 0);

//     //copy seq
//     copy_kseq_t(&rev, seq);

//     for (size_t i = 0; i < seq->seq.l; ++i)
//       rev.seq.s[i] = complement(seq->seq.s[seq->seq.l - i - 1]);

//     if (rev.seq.m > rev.seq.l)
//       rev.seq.s[rev.seq.l] = 0;

//     bool rev_align = p->aligner->align(&rev, sam_fd, 1);

//     if (fwd_align or rev_align)
//       n_aligned_reads++;
//     n_reads++;

//     free(rev.name.s);
//     free(rev.comment.s);
//     free(rev.seq.s);
//     free(rev.qual.s);
//   }

//   verbose("Number of aligned reads block ", p->wk_id, " : ", n_aligned_reads, "/", n_reads);
//   p->n_reads = n_reads;
//   p->n_aligned_reads = n_aligned_reads;
//   kseq_destroy(seq);
//   gzclose(fp);
//   fclose(sam_fd);

//   return NULL;
// }
// template <typename aligner_t>
// size_t mt_align(aligner_t *aligner, std::string pattern_filename, std::string sam_filename, size_t n_threads)
// {
//   pthread_t t[n_threads] = {0};
//   mt_param_t<aligner_t> params[n_threads];
//   std::vector<size_t> starts = split_fastq(pattern_filename, n_threads);
//   for(size_t i = 0; i < n_threads; ++i)
//   {
//     params[i].aligner = aligner;
//     params[i].pattern_filename = pattern_filename;
//     params[i].sam_filename = sam_filename + "_" + std::to_string(i) + ".sam";
//     params[i].start = starts[i];
//     params[i].end = starts[i+1];
//     params[i].wk_id = i;
//     xpthread_create(&t[i], NULL, &mt_align_worker<aligner_t>, &params[i], __LINE__, __FILE__);
//   }

//   size_t tot_reads = 0;
//   size_t tot_aligned_reads = 0;

//   for(size_t i = 0; i < n_threads; ++i)
//   {
//     xpthread_join(t[i],NULL,__LINE__,__FILE__);
//   }

//   sleep(5);
//   for(size_t i = 0; i < n_threads; ++i)
//   {
//     tot_reads += params[i].n_reads;
//     tot_aligned_reads += params[i].n_aligned_reads;
//   }

//   verbose("Number of aligned reads: ", tot_aligned_reads, "/", tot_reads);
//   return tot_aligned_reads;
// }

template <typename aligner_t>
void *mt_align_worker(void *param)
{
  mt_param_t<aligner_t> *p = (mt_param_t<aligner_t> *)param;
  size_t n_reads = 0;
  size_t n_aligned_reads = 0;

  FILE *sam_fd;
  gzFile fp;

  if ((sam_fd = fopen(p->sam_filename.c_str(), "w")) == nullptr)
    error("open() file " + p->sam_filename + " failed");

  size_t b_size = p->b_size;
  if (p->seq != nullptr)
  {
    kseq_t *seq = p->seq;
    kbseq_t *b = kbseq_init();
    int l = 0;
    while ((l = mt_kbseq_read(b, seq, b_size)) > 0)
    {
      for (size_t i = 0; i < l; ++i)
      {
        if (p->aligner->align(&b->buf[i], sam_fd))
          n_aligned_reads++;        
        n_reads++;
      }
      // std::cout << "Block p " << p->wk_id << " end!" << std::endl;
    }
    kbseq_destroy(b);
  }
  else
  {
    kseq_t *mate1 = p->mate1;
    kseq_t *mate2 = p->mate2;
    kpbseq_t *b = kpbseq_init();
    int l = 0;
    while ((l = mt_kpbseq_read(b, mate1, mate2, b_size)) > 0)
    {
      // for (size_t i = 0; i < l; ++i)
      // {
      //   if (p->aligner->align(&b->mate1->buf[i], &b->mate2->buf[i],sam_fd))
      //     n_aligned_reads++;
      //   n_reads++;
      // }
      // std::cout << "Block p " << p->wk_id << " end!" << std::endl;
      n_aligned_reads+=p->aligner->align(b,sam_fd);
      n_reads += l;
    }
    kpbseq_destroy(b);
  }

  verbose("Number of aligned reads block ", p->wk_id, " : ", n_aligned_reads, "/", n_reads);
  p->n_reads = n_reads;
  p->n_aligned_reads = n_aligned_reads;
  fclose(sam_fd);

  return NULL;
}

template <typename aligner_t>
size_t mt_align(aligner_t *aligner, std::string pattern_filename, std::string sam_filename, size_t n_threads, size_t b_size, std::string mate2_filename = "")
{
  xpthread_mutex_init(&mutex_reads_dispatcher, NULL, __LINE__, __FILE__);
  kseq_t *seq = nullptr;
  kseq_t *mate1 = nullptr;
  kseq_t *mate2 = nullptr;

  gzFile fp_mate2 = nullptr;
  gzFile fp = gzopen(pattern_filename.c_str(), "r");
  if (mate2_filename != "")
  {
    fp_mate2 = gzopen(mate2_filename.c_str(), "r");
    mate1 = kseq_init(fp);
    mate2 = kseq_init(fp_mate2);
  }
  else
  {
    seq = kseq_init(fp);
  }

  pthread_t t[n_threads] = {0};
  mt_param_t<aligner_t> params[n_threads];
  for (size_t i = 0; i < n_threads; ++i)
  {
    // Create a new thread
    params[i].aligner = aligner;
    params[i].pattern_filename = pattern_filename;
    params[i].sam_filename = sam_filename + "_" + std::to_string(i) + ".sam";
    params[i].b_size = b_size;
    params[i].seq = seq;
    params[i].mate1 = mate1;
    params[i].mate2 = mate2;
    params[i].wk_id = i;
    xpthread_create(&t[i], NULL, &mt_align_worker<aligner_t>, &params[i], __LINE__, __FILE__);
  }

  size_t tot_reads = 0;
  size_t tot_aligned_reads = 0;

  for (size_t i = 0; i < n_threads; ++i)
  {
    xpthread_join(t[i], NULL, __LINE__, __FILE__);
  }

  if (fp_mate2 != nullptr)
  {
    kseq_destroy(mate1);
    kseq_destroy(mate2);
    gzclose(fp_mate2);
  }
  else
  {
    kseq_destroy(seq);
  }
  gzclose(fp);

  // sleep(5);
  verbose("Merging temporary SAM files");

  FILE *fd;

  if ((fd = fopen(std::string(sam_filename + ".sam").c_str(), "w")) == nullptr)
    error("open() file " + std::string(sam_filename + ".sam") + " failed");

  fprintf(fd, "%s", aligner->to_sam().c_str());

  for (size_t i = 0; i < n_threads; ++i)
  {
    tot_reads += params[i].n_reads;
    tot_aligned_reads += params[i].n_aligned_reads;

    append_file(params[i].sam_filename, fd);
    if (std::remove(params[i].sam_filename.c_str()) != 0)
      error("remove() file " + params[i].sam_filename + " failed");
  }

  xpthread_mutex_destroy(&mutex_reads_dispatcher, __LINE__, __FILE__);

  verbose("Number of aligned reads: ", tot_aligned_reads, "/", tot_reads);
  return tot_aligned_reads;
}

////////////////////////////////////////////////////////////////////////////////
/// Single Thread
////////////////////////////////////////////////////////////////////////////////
template <typename aligner_t>
size_t st_align(aligner_t *aligner, std::string pattern_filename, std::string sam_filename, size_t b_size, std::string mate2_filename = "")
{
  size_t n_reads = 0;
  size_t n_aligned_reads = 0;

  kseq_t *seq = nullptr;
  kseq_t *mate1 = nullptr;
  kseq_t *mate2 = nullptr;

  gzFile fp_mate2 = nullptr;
  gzFile fp = gzopen(pattern_filename.c_str(), "r");
  if (mate2_filename != "")
  {
    fp_mate2 = gzopen(mate2_filename.c_str(), "r");
    mate1 = kseq_init(fp);
    mate2 = kseq_init(fp_mate2);
  }
  else
  {
    seq = kseq_init(fp);
  }

  int l;
  FILE *sam_fd;

  sam_filename += ".sam";

  if ((sam_fd = fopen(sam_filename.c_str(), "w")) == nullptr)
    error("open() file " + sam_filename + " failed");

  fprintf(sam_fd, "%s", aligner->to_sam().c_str());

  if (seq != nullptr)
  {
    kbseq_t *b = kbseq_init();
    int l = 0;
    while ((l = kbseq_read(b, seq, b_size)) > 0)
    {
      for (size_t i = 0; i < l; ++i)
      {
        if (aligner->align(&b->buf[i], sam_fd))
          n_aligned_reads++;
        n_reads++;
      }
    }
    kbseq_destroy(b);
  } 
  else
  {
    kpbseq_t *b = kpbseq_init();
    int l = 0;
    while ((l = kpbseq_read(b, mate1, mate2, b_size)) > 0)
    {
      // for (size_t i = 0; i < l; ++i)
      // {
      //   if (aligner->align(&b->mate1->buf[i], &b->mate2->buf[i], sam_fd))
      //     n_aligned_reads++;
      //   n_reads++;
      // }
      // std::cout << "Block p " << p->wk_id << " end!" << std::endl;
      n_aligned_reads += aligner->align(b,sam_fd);
      n_reads += l;
    }
    kpbseq_destroy(b);
  }
  verbose("Number of aligned reads: ", n_aligned_reads, "/", n_reads);
  kseq_destroy(seq);
  gzclose(fp);
  fclose(sam_fd);

  return n_aligned_reads;
}



#endif /* end of include guard: _READS_DISPATCHER_HH */
