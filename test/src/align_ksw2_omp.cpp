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

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/io.hpp>

#include <ms_pointers.hpp>

#include <ksw2.h>
// #include <malloc_count.h>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>


#include <omp.h>

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

typedef std::pair<std::string, std::vector<uint8_t>> pattern_t;

struct Patterns
{
private:
  std::ifstream input_file;

public:
  using pattern_t = std::pair<std::string, std::vector<uint8_t>>;

  Patterns(std::string input_path) : input_file(input_path) {}

  pattern_t next();
  bool end();

  size_t n_patterns = 0;
};

Patterns::pattern_t Patterns::next()
{
  n_patterns++;
  std::string header, read;
  std::vector<uint8_t> pattern;
  // header
  std::getline(input_file, header);
  // read
  std::getline(input_file, read);

  pattern_t out;
  out.first = header;
  std::copy(read.begin(), read.end(), std::back_inserter(out.second));
  return out;
}

bool Patterns::end()
{
  return input_file.eof();
}

class aligner_t
{
public:
  using SelSd = SelectSdvec<>;
  using DagcSd = DirectAccessibleGammaCode<SelSd>;

  aligner_t(std::string filename, 
            size_t min_len_ = 50,
            bool forward_only_ = true) : 
                  min_len(min_len_),
                  forward_only(forward_only_)
  {
    verbose("Loading the matching statistics index");
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    std::string filename_ms = filename + ".ms";

    ifstream fs_ms(filename_ms);
    ms.load(fs_ms);

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Matching statistics index construction complete");
    // verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    verbose("Loading random access");
    t_insert_start = std::chrono::high_resolution_clock::now();

    std::string filename_slp = filename + ".slp";

    ifstream fs(filename_slp);
    ra.load(fs);

    n = ra.getLen();

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Matching statistics index loading complete");
    // verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    verbose("Initialize the local aligner");
    t_insert_start = std::chrono::high_resolution_clock::now();

#ifdef HAVE_KALLOC
    km = no_kalloc ? 0 : km_init();
#endif

    ksw_gen_simple_mat(5, mat, a, -b);

    t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Local aligner initialization complete");
    // verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
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

      // Convert A,C,G,T,N into 0,1,2,3,4
      for (i = 0; i < (int)read->seq.l; ++i)
        read->seq.s[i] = seq_nt4_table[(int)read->seq.s[i]];

      for (i = 0; i < (int)len; ++i)
        str[i] = seq_nt4_table[(int)str[i]];

      ksw_extz_t ez;
      ksw_reset_extz(&ez);

      global_aln(algo, km, read->seq.s, str, 5, mat, q, e, q2, e2, w, zdrop, flag, &ez);

      print_aln("first", "second", &ez);


      // if(ez.n_cigar)free(ez.cigar);

      // int score;

      // score = ksw_global(read->seq.l, (uint8_t *)read->seq.s, len, (uint8_t *)str, 5, mat, gapo, gape, w, &q[0]);

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
      // alignment.ref_begin += left_occ;
      // alignment.ref_end += left_occ;
      // alignment.ref_end_next_best += left_occ;

      // if(alignment.sw_score >= min_score)
      //   ssw_write_sam(alignment,"human",read,strand,out);

      // aligned_reads++;
      aligned = true;

      delete str;
    }
    return aligned;
  }

  size_t get_aligned_reads()
  {
    return aligned_reads;
  }

  // // Adapted from SSW
  // static void ssw_write_sam(StripedSmithWaterman::Alignment &a,
  //                           const char *ref_seq_name,
  //                           const kseq_t *read,
  //                           int8_t strand,
  //                           FILE* out) // 0: forward aligned ; 1: reverse complement aligned
  // {
  //   // Sam format output
  //   fprintf(out, "%s\t", read->name.s);
  //   if (a.sw_score == 0)
  //     fprintf(out, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
  //   else
  //   {
  //     int32_t c, p;
  //     uint32_t mapq = -4.343 * log(1 - (double)abs(a.sw_score - a.sw_score_next_best) / (double)a.sw_score);
  //     mapq = (uint32_t)(mapq + 4.99);
  //     mapq = mapq < 254 ? mapq : 254;
  //     if (strand)
  //       fprintf(out, "16\t");
  //     else
  //       fprintf(out, "0\t");
  //     // TODO: Find the correct reference name.
  //     fprintf(out, "%s\t%d\t%d\t", ref_seq_name, a.ref_begin + 1, mapq);
  //     // mismatch = mark_mismatch(a.ref_begin, a.query_begin, a.query_end, ref_num, read_num, read->seq.l, &a->cigar, &a->cigarLen);
  //     // for (c = 0; c < a->cigarLen; ++c)
  //     // {
  //     //   char letter = cigar_int_to_op(a->cigar[c]);
  //     //   uint32_t length = cigar_int_to_len(a->cigar[c]);
  //     //   fprintf(out, "%lu%c", (unsigned long)length, letter);
  //     // }
  //     fprintf(out, "%s", a.cigar_string.c_str());
  //     fprintf(out, "\t*\t0\t0\t");
  //     fprintf(out, "%s", read->seq.s);
  //     fprintf(out, "\t");
  //     if (read->qual.s && strand)
  //     {
  //       for (p = read->qual.l - 1; p >= 0; --p)
  //         fprintf(out, "%c", read->qual.s[p]);
  //     }
  //     else if (read->qual.s)
  //       fprintf(out, "%s", read->qual.s);
  //     else
  //       fprintf(out, "*");
  //     fprintf(out, "\tAS:i:%d", a.sw_score);
  //     fprintf(out, "\tNM:i:%d\t", a.mismatches);
  //     if (a.sw_score_next_best > 0)
  //       fprintf(out, "ZS:i:%d\n", a.sw_score_next_best);
  //     else
  //       fprintf(out, "\n");
  //   }
  // }

  // *************** FROM Ksw2 cli.c *******************************************
  static void print_aln(const char *tname, const char *qname, ksw_extz_t *ez)
  {
    printf("%s\t%s\t%d", tname, qname, ez->score);
    printf("\t%d\t%d\t%d", ez->max, ez->max_t, ez->max_q);
    if (ez->n_cigar > 0)
    {
      int i;
      putchar('\t');
      for (i = 0; i < ez->n_cigar; ++i)
        printf("%d%c", ez->cigar[i] >> 4, "MID"[ez->cigar[i] & 0xf]);
    }
    putchar('\n');
  }

  static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
  {
    int i, j;
    a = a < 0 ? -a : a;
    b = b > 0 ? -b : b;
    for (i = 0; i < m - 1; ++i)
    {
      for (j = 0; j < m - 1; ++j)
        mat[i * m + j] = i == j ? a : b;
      mat[i * m + m - 1] = 0;
    }
    for (j = 0; j < m; ++j)
      mat[(m - 1) * m + j] = 0;
  }

  static void global_aln(const char *algo, void *km, const char *qseq_, const char *tseq_, int8_t m, const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2,
                         int w, int zdrop, int flag, ksw_extz_t *ez)
  {
    int i, qlen, tlen;
    uint8_t *qseq, *tseq;
    ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
    ez->max = 0, ez->mqe = ez->mte = KSW_NEG_INF;
    ez->n_cigar = 0;
    ez->m_cigar = 0;
    qlen = strlen(qseq_);
    tlen = strlen(tseq_);
    qseq = (uint8_t *)calloc(qlen + 33, 1); // 32 for gaba
    tseq = (uint8_t *)calloc(tlen + 33, 1);
    for (i = 0; i < qlen; ++i)
      qseq[i] = seq_nt4_table[(uint8_t)qseq_[i]];
    for (i = 0; i < tlen; ++i)
      tseq[i] = seq_nt4_table[(uint8_t)tseq_[i]];
    if (strcmp(algo, "gg") == 0)
    {
      if (flag & KSW_EZ_SCORE_ONLY)
        ez->score = ksw_gg(km, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, m, mat, q, e, w, 0, 0, 0);
      else
        ez->score = ksw_gg(km, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, m, mat, q, e, w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    }
    else if (strcmp(algo, "gg2") == 0)
    {
      if (flag & KSW_EZ_SCORE_ONLY)
        ez->score = ksw_gg2(km, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, m, mat, q, e, w, 0, 0, 0);
      else
        ez->score = ksw_gg2(km, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, m, mat, q, e, w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    }
    else if (strcmp(algo, "gg2_sse") == 0)
      ez->score = ksw_gg2_sse(km, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, m, mat, q, e, w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    else if (strcmp(algo, "extz") == 0)
      ksw_extz(km, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, m, mat, q, e, w, zdrop, flag, ez);
    else if (strcmp(algo, "extz2_sse") == 0)
      ksw_extz2_sse(km, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, m, mat, q, e, w, zdrop, 0, flag, ez);
    else if (strcmp(algo, "extd") == 0)
      ksw_extd(km, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, m, mat, q, e, q2, e2, w, zdrop, flag, ez);
    else if (strcmp(algo, "extd2_sse") == 0)
      ksw_extd2_sse(km, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, m, mat, q, e, q2, e2, w, zdrop, 0, flag, ez);
    else if (strcmp(algo, "extf2_sse") == 0)
      ksw_extf2_sse(km, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, mat[0], mat[1], e, w, zdrop, ez);
    else if (strcmp(algo, "exts2_sse") == 0)
    {
      int8_t mat[25];
      ksw_gen_simple_mat(5, mat, 1, 2);
      ksw_exts2_sse(km, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, m, mat, 2, 1, 32, 4, zdrop, flag | KSW_EZ_SPLICE_FOR, ez);
    }
    else if (strcmp(algo, "test") == 0)
      ksw_extd2_sse(km, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, m, mat, 4, 2, 24, 1, 751, 400, 0, 8, ez);
#ifdef HAVE_GABA
    else if (strcmp(algo, "gaba") == 0)
    { // libgaba. Note that gaba may not align to the end
      int buf_len = 0x10000;
      char *cigar = (char *)calloc(buf_len, 1);
      for (i = 0; i < qlen; ++i)
        qseq[i] = qseq[i] < 4 ? 1 << qseq[i] : 15;
      for (i = 0; i < tlen; ++i)
        tseq[i] = tseq[i] < 4 ? 1 << tseq[i] : 15;
      struct gaba_section_s qs = gaba_build_section(0, qseq, qlen);
      struct gaba_section_s ts = gaba_build_section(2, tseq, tlen);

      gaba_t *ctx = gaba_init(GABA_PARAMS(.xdrop = (zdrop < 120 ? zdrop : 120), GABA_SCORE_SIMPLE(mat[0], abs(mat[1]), q, e)));
      void const *lim = (void const *)0x800000000000;
      gaba_dp_t *dp = gaba_dp_init(ctx, lim, lim);
      struct gaba_fill_s *f = gaba_dp_fill_root(dp, &qs, 0, &ts, 0);
      f = gaba_dp_fill(dp, f, &qs, &ts);
      gaba_alignment_t *r = gaba_dp_trace(dp, 0, f, 0);
      ez->score = r->score;
      gaba_dp_dump_cigar_forward(cigar, buf_len, r->path->array, 0, r->path->len);
      printf("gaba cigar: %s\n", cigar);
      gaba_dp_clean(dp);
      gaba_clean(ctx);
      free(cigar);
    }
#endif
#ifdef HAVE_PARASAIL
    else if (strncmp(algo, "ps_", 3) == 0)
    {
      parasail_matrix_t *ps_mat;
      parasail_function_t *func;
      parasail_result_t *res;
      ps_mat = parasail_matrix_create("ACGTN", mat[0], mat[1]); // TODO: call parasail_matrix_set_value() to change N<->A/C/G/T scores
      func = parasail_lookup_function(algo + 3);
      if (func == 0)
      {
        fprintf(stderr, "ERROR: can't find parasail function '%s'\n", algo + 3);
        exit(1);
      }
      res = func(qseq_, qlen, tseq_, tlen, q + e, e, ps_mat);
      ez->score = res->score;
      parasail_matrix_free(ps_mat);
      parasail_result_free(res);
    }
#endif
    else
    {
      fprintf(stderr, "ERROR: can't find algorithm '%s'\n", algo);
      exit(1);
    }
    free(qseq);
    free(tseq);
  }
  // *************** END FROM Ksw2 cli.c ***************************************

protected:
  ms_pointers<> ms;
  SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd> ra;

  size_t min_len = 0;
  size_t aligned_reads = 0;
  size_t n = 0;

  void *km = 0;
  int8_t a = 2, b = 4, q = 4, e = 2, q2 = 13, e2 = 1;
  int c, i, pair = 1, w = -1, flag = 0, rep = 1, zdrop = -1, no_kalloc = 0;
  char *algo = "gg", *s;
  int8_t mat[25];

  bool forward_only = true;
};

//*********************** Argument options ***************************************
// struct containing command line parameters and other globals
struct Args
{
  std::string filename = "";
  size_t w = 10;             // sliding window size and its default
  bool store = false;        // store the data structure in the file
  bool memo = false;         // print the memory usage
  bool csv = false;          // print stats on stderr in csv format
  bool rle = false;          // outpt RLBWT
  std::string patterns = ""; // path to patterns file
  size_t l = 25;             // minumum MEM length
  size_t th = 1;             // number of threads
  bool is_fasta = false;     // read a fasta file
};

void parseArgs(int argc, char *const argv[], Args &arg)
{
  int c;
  extern char *optarg;
  extern int optind;

  std::string usage("usage: " + std::string(argv[0]) + " infile [-s store] [-m memo] [-c csv] [-p patterns] [-f fasta] [-r rle] [-t threads] [-l len]\n\n" +
                    "Computes the pfp data structures of infile, provided that infile.parse, infile.dict, and infile.occ exists.\n" +
                    "  wsize: [integer] - sliding window size (def. 10)\n" +
                    "  store: [boolean] - store the data structure in infile.pfp.ds. (def. false)\n" +
                    "   memo: [boolean] - print the data structure memory usage. (def. false)\n" +
                    "  fasta: [boolean] - the input file is a fasta file. (def. false)\n" +
                    "    rle: [boolean] - output run length encoded BWT. (def. false)\n" +
                    "pattens: [string]  - path to patterns file.\n" +
                    "    len: [integer] - minimum MEM lengt (def. 25)\n" +
                    " thread: [integer] - number of threads (def. 1)\n" +
                    "    csv: [boolean] - print the stats in csv form on strerr. (def. false)\n");

  std::string sarg;
  while ((c = getopt(argc, argv, "w:smcfl:rhp:t:")) != -1)
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

int
main(int argc, char *const argv[])
{

  Args args;
  parseArgs(argc, argv, args);

  verbose("Construction of the aligner");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  aligner_t aligner(args.filename, 75);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
  // verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();
  
  size_t n_reads = 0;
  size_t n_aligned_reads = 0;

// #pragma omp parallel
// {
// #pragma omp parallel for schedule(static)
  for(size_t i = 1; i <= 64; ++i)
  {
    size_t n_reads_p = 0;
    size_t n_aligned_reads_p = 0;

    gzFile fp;
    kseq_t *seq;
    int l;

    std::string file_path = args.patterns + "_" + std::to_string(i) + ".fa";

    // Open out file
    std::string filename = file_path + ".sam";
    FILE *sam_fd;

    if ((sam_fd = fopen(filename.c_str(), "w")) == nullptr)
      error("open() file " + filename + " failed");

    fp = gzopen(file_path.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0)
    {
      if(aligner.align(seq,sam_fd,0))
        n_aligned_reads_p ++;
      n_reads_p++;

      // std::cout << "\rSequenced patterns on block " << i << " : "
      //           << n_reads_p << std::flush;
    }

    verbose("Number of aligned reads block ", i, " : ", n_aligned_reads_p, "/", n_reads_p);
    kseq_destroy(seq);
    gzclose(fp);
    fclose(sam_fd);

  }

// #pragma omp critical
//   {
//     // Count the total number of reads.
//     n_reads += n_reads_p;
//     // Count the total number of aligned reads
//     n_aligned_reads += n_aligned_reads_p;
//   }
// }




  // TODO: Merge the SAM files.

  t_insert_end = std::chrono::high_resolution_clock::now();

  // verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  // auto mem_peak = malloc_count_peak();
  // verbose("Memory peak: ", malloc_count_peak());
  verbose("Number of aligned reads: ", n_aligned_reads, "/", n_reads);

  size_t space = 0;
  if (args.memo)
  {
  }

  if (args.store)
  {
  }

  // if (args.csv)
    // std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;
}