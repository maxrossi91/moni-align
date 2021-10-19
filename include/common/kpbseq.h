/* kpbseq - read a batch of reads
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
   \file kpbseq.cpp
   \brief kpbseq.cpp read a batch of reads.
   \author Massimiliano Rossi
   \date 03/06/2021
*/

#ifndef _KPBSEQ_HH
#define _KPBSEQ_HH

#include <iostream>
#include <vector>

#include "kseq.h"
#include <zlib.h>

// KSEQ_INIT(gzFile, gzread);

////////////////////////////////////////////////////////////////////////////////
/// kseq extra
////////////////////////////////////////////////////////////////////////////////
static inline size_t ks_tell(kseq_t *seq)
{
    return gztell(seq->f->f) - seq->f->end + seq->f->begin;
}

// test if the file is gzipped
static inline bool is_gzipped(std::string filename)
{
    FILE *fp = fopen(filename.c_str(), "rb");
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

static inline std::string to_string(kstring_t &s)
{
    return std::string(s.s, s.l);
}

static inline void ks_print(kseq_t *seq)
{
    std::cout << "Name: " << to_string(seq->name)
              << "\nComment: " << to_string(seq->comment)
              << "\nSeq:  " << to_string(seq->seq)
              << "\nQual: " << to_string(seq->qual)
              << std::endl;
}

void print_kseq_from_to(gzFile fp, size_t start, size_t end)
{
    gzseek(fp, start, SEEK_SET);
    kseq_t *seq = kseq_init(fp);
    int l;
    while ((ks_tell(seq) < end) && ((l = kseq_read(seq)) >= 0))
    {
        ks_print(seq);
    }
}

inline static void copy_kstring_t(kstring_t &l, kstring_t &r)
{
    l.l = r.l;
    l.m = r.m;
    l.s = (char *)malloc(l.m);
    for (size_t i = 0; i < r.m; ++i)
        l.s[i] = r.s[i];
}
inline static void copy_kseq_t(kseq_t *l, kseq_t *r)
{
    copy_kstring_t(l->name, r->name);
    copy_kstring_t(l->comment, r->comment);
    copy_kstring_t(l->seq, r->seq);
    copy_kstring_t(l->qual, r->qual);
    l->last_char = r->last_char;
}

static const unsigned char seq_compl_table[256] = {
      0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
     16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
     32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
     48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
     64,'T', 66,'G', 68, 69, 70,'C', 72, 73, 74, 75, 76, 77, 78, 79,
     80, 81, 82, 83,'A', 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
     96,'T', 98,'G',100,101,102,'C',104,105,106,107,108,109,110,111,
    112,113,114,115,'A',117,118,119,120,121,122,123,124,125,126,127,
    128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
    144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,
    160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
    176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,
    192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,
    208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,
    224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,
    240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,
    };

inline static void r_copy_kstring_t(kstring_t &l, kstring_t &r)
{
    l.l = r.l;
    l.m = r.m;
    l.s = (char *)malloc(l.m);
    for (size_t i = 0; i < r.m; ++i)
        l.s[i] = r.s[r.l - i - 1];
    if (l.m > l.l)
      l.s[l.l] = 0;
}

inline static void rc_copy_kstring_t(kstring_t &l, kstring_t &r)
{
    l.l = r.l;
    l.m = r.m;
    l.s = (char *)malloc(l.m);
    for (size_t i = 0; i < r.m; ++i)
        l.s[i] = seq_compl_table[r.s[r.l - i - 1]];
    if (l.m > l.l)
      l.s[l.l] = 0;
}

inline static void rc_copy_kseq_t(kseq_t *l, kseq_t *r)
{
    copy_kstring_t(l->name, r->name);
    copy_kstring_t(l->comment, r->comment);
    rc_copy_kstring_t(l->seq, r->seq);
    r_copy_kstring_t(l->qual, r->qual);
    l->last_char = r->last_char;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Batch kseq
////////////////////////////////////////////////////////////////////////////////

typedef struct
{
    size_t l, m;
    kseq_t *buf;
} kbseq_t;

inline static kbseq_t *kbseq_init()
{
    kbseq_t *b = (kbseq_t *)malloc(sizeof(kbseq_t));
    b->m = 0;
    b->l = 0;
    b->buf = (kseq_t *)malloc(sizeof(kseq_t));
    return b;
}

inline static void kbseq_realloc(kbseq_t *b)
{
    b->m = b->l + 1;
    kroundup32(b->m);
    b->buf = (kseq_t *)realloc(b->buf, b->m * sizeof(kseq_t));
}

/**
 * \brief append seq to b.
 * 
 * \param b the batch o sequences
 * \param seq the sequence to append.
 */
inline static void kbseq_push_back(kbseq_t *b, kseq_t *seq)
{
    // allocate new memory if necessary
    if (b->l >= b->m)
        kbseq_realloc(b);
    copy_kseq_t(&b->buf[b->l++], seq);
}

/**
 * \brief Reset the batch b
 * 
 * \param b the batch of sequences.
 */
inline static void kbseq_reset(kbseq_t *b)
{
    for (size_t i = 0; i < b->l; ++i)
    {
        free(b->buf[i].name.s);
        free(b->buf[i].comment.s);
        free(b->buf[i].seq.s);
        free(b->buf[i].qual.s);
    }
    b->l = 0;
}

/**
 * \brief Read n sequences from seq and append them to b.
 * 
 * \param b the batch of sequences
 * \param seq the sequens to read from
 * \param n the number of sequences to read.
 * \return size_t the number of sequences read.
 */
inline static size_t kbseq_read(kbseq_t *b, kseq_t *seq, size_t n)
{
    int l = 0;
    size_t i = 0;
    kbseq_reset(b);
    while ((i < n) and ((l = kseq_read(seq)) >= 0))
    {
        i++;
        kbseq_push_back(b, seq);
    }
    return i;
}

/**
 * @brief destroy kbseq b
 * 
 * @param b the kbseq to be destroyed
 */
inline static void kbseq_destroy(kbseq_t *b)
{
    kbseq_reset(b);
    delete[] b->buf;
    delete b;
}

/**
 * @brief Print a kbseq
 * 
 * @param b the kbseq to be print
 */
static inline void kbs_print(kbseq_t *b)
{
    for(size_t i = 0; i < b->l; ++i )
        ks_print(&b->buf[i]);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Paired batch kseq
////////////////////////////////////////////////////////////////////////////////

typedef struct
{
    kbseq_t* mate1;
    kbseq_t* mate2;
} kpbseq_t;

inline static kpbseq_t *kpbseq_init()
{
    kpbseq_t *p = (kpbseq_t *)malloc(sizeof(kpbseq_t));
    p->mate1 = kbseq_init();
    p->mate2 = kbseq_init();
    return p;
}

/**
 * \brief Read n sequences from seq and append them to b.
 * 
 * \param p the paired batch of sequences
 * \param mate_1 the sequens to read from
 * \param mate_2 the sequens to read from
 * \param n the number of sequences to read.
 * \return size_t the number of sequences read.
 */
inline static size_t kpbseq_read(kpbseq_t *p, kseq_t *mate1_s, kseq_t *mate2_s, size_t n)
{
    size_t l1 = kbseq_read(p->mate1, mate1_s, n);
    size_t l2 = kbseq_read(p->mate2, mate2_s, n);
    if(l1 != l2)
        error("The paired-end files does not have the same number of sequences!");
    return l1;
}

inline static void kpbseq_destroy(kpbseq_t *p)
{
    kbseq_destroy(p->mate1);
    kbseq_destroy(p->mate2);
    delete p;
}

////////////////////////////////////////////////////////////////////////////////


#endif /* end of include guard: _KPBSEQ_HH */
