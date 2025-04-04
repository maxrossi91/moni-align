/* pfp-ds - prefix free parsing data structures
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
   \file common.hpp
   \brief common.hpp contains common features.
   \author Massimiliano Rossi
   \date 12/03/2020
*/

#ifndef _COMMON_HH
#define _COMMON_HH

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <assert.h>

#include <sys/time.h>

#include <sys/mman.h> // for mmap
#include <unistd.h>
#include <sys/stat.h>
 #include <fcntl.h>

#include <sstream>      // std::stringstream

#include <vector>      // std::vector

#include <chrono>       // high_resolution_clock

#include <sdsl/io.hpp>  // serialize and load
#include <type_traits>  // enable_if_t and is_fundamental

#include <mutex>

#include <stdexcept>
#include <execinfo.h>

//**************************** From  Big-BWT ***********************************
// special symbols used by the construction algorithm:
//   they cannot appear in the input file
//   the 0 symbol is used in the final BWT file as the EOF char

#define Dollar 2     // special char for the parsing algorithm, must be the highest special char
#define EndOfWord 1  // word delimiter for the plain dictionary file
#define EndOfDict 0  // end of dictionary delimiter
//******************************************************************************

#define THRBYTES 5 // The number of bytes for the thresholds
#define SSABYTES 5 // The number of bytes for the thresholds

std::string NowTime();
void _internal_messageInfo(const std::string message);
void _internal_messageWarning( const std::string file, const unsigned int line, const std::string message);
void _internal_messageError( const std::string file, const unsigned int line,const std::string message);


std::string NowTime()
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    char buffer[100];
    tm r;
    strftime(buffer, sizeof(buffer), "%X", localtime_r(&tv.tv_sec, &r));
    char result[100];
    snprintf(result, 100, "%s"/*.%06ld"*/, buffer/*, (long)tv.tv_usec*/);
    return result;
}


template<typename T>
inline void _internal_message_helper(std::stringstream &ss, T const &first) { ss << first; }
template<typename T, typename... Args>
inline void _internal_message_helper(std::stringstream &ss, T const &first, const Args&... args) { ss << first << " "; _internal_message_helper(ss,args...); }
template<typename T, typename... Args>
inline std::string _internal_message(T const &first, const Args&... args) { std::stringstream ss; _internal_message_helper(ss,first,args...); return ss.str(); }


void _internal_messageInfo(const std::string message)
{
  std::cout << "[INFO] " << NowTime() << " - " << "Message: " << message << std::endl;
}

void _internal_messageWarning( const std::string file, const unsigned int line,
  const std::string message)
{
  std::cout << "[WARNING] " << NowTime() << " - "
  << "File: " << file << '\n'
  << "Line: " << line << '\n'
  << "Message: " << message << std::endl;
}

void _internal_messageError( const std::string file, const unsigned int line,
  const std::string message)
{
  std::cerr << "[ERROR] " << NowTime() << " - "
  << "File: " << file << '\n'
  << "Line: " << line << '\n'
  << "Message: " << message << std::endl;
  assert( false );
  exit( 1 );
}



#define info( args... ) \
    _internal_messageInfo( _internal_message(args) )

#ifdef VERBOSE
  #define verbose( args... ) \
      _internal_messageInfo( _internal_message(args) )
#else
  #define verbose( args... )
#endif

#define warning( args... ) \
    _internal_messageWarning( __FILE__, __LINE__, _internal_message(args) )

#define error( args... ) \
    _internal_messageError( __FILE__, __LINE__, _internal_message(args) )


// converts elemens in csv format
template <typename T>
inline void csv_helper(std::stringstream &ss, T const &first){ss << first;}
template <typename T, typename... Args>
inline void csv_helper(std::stringstream &ss, T const &first, const Args &... args){ ss << first << ", "; csv_helper(ss, args...);}
template <typename T, typename... Args>
inline std::string csv(T const &first, const Args &... args){std::stringstream ss;csv_helper(ss, first, args...); return ss.str();}

//*********************** File I/O *********************************************
template<typename T>
void map_file(const char *filename, T*& ptr, size_t& length){
    struct stat filestat;
    int fd;

    if ((fd = open(filename, O_RDONLY)) < 0)
        error("open() file " + std::string(filename) + " failed" );

    if (fstat(fd, &filestat) < 0)
        error("stat() file " + std::string(filename) + " failed" );

    if(filestat.st_size % sizeof(T) != 0)
        error("invilid file " + std::string(filename));

    length = filestat.st_size / sizeof(T);

    if ((ptr = mmap(NULL, filestat.st_size, PROT_READ, MAP_SHARED, fd, 0)) == MAP_FAILED)
        error("mmap() file " + std::string(filename) + " failed");
}

template<typename T>
void read_file(const char *filename, T*& ptr, size_t& length){
    struct stat filestat;
    FILE* fd;

    if ((fd = fopen(filename, "r")) == nullptr)
        error("open() file " + std::string(filename) + " failed" );

    int fn = fileno(fd);
    if (fstat(fn, &filestat) < 0)
        error("stat() file " + std::string(filename) + " failed" );

    if(filestat.st_size % sizeof(T) != 0)
        error("invilid file " + std::string(filename));

    length = filestat.st_size / sizeof(T);
    ptr = new T[length];

    if ((fread(ptr, sizeof(T), length, fd)) != length)
        error("fread() file " + std::string(filename) + " failed");

    fclose(fd);
}

template<typename T>
void read_file(const char *filename, std::vector<T>& ptr){
    struct stat filestat;
    FILE* fd;

    if ((fd = fopen(filename, "r")) == nullptr)
        error("open() file " + std::string(filename) + " failed" );

    int fn = fileno(fd);
    if (fstat(fn, &filestat) < 0)
        error("stat() file " + std::string(filename) + " failed" );

    if(filestat.st_size % sizeof(T) != 0)
        error("invilid file " + std::string(filename));

    size_t length = filestat.st_size / sizeof(T);
    ptr.resize(length);

    if ((fread(&ptr[0], sizeof(T), length, fd)) != length)
        error("fread() file " + std::string(filename) + " failed");

    fclose(fd);
}

void read_file(const char *filename, std::string &ptr)
{
  struct stat filestat;
  FILE *fd;

  if ((fd = fopen(filename, "r")) == nullptr)
    error("open() file " + std::string(filename) + " failed");

  int fn = fileno(fd);
  if (fstat(fn, &filestat) < 0)
    error("stat() file " + std::string(filename) + " failed");

  if (filestat.st_size % sizeof(char) != 0)
    error("invilid file " + std::string(filename));

  size_t length = filestat.st_size / sizeof(char);
  ptr.resize(length);

  if ((fread(&ptr[0], sizeof(char), length, fd)) != length)
    error("fread() file " + std::string(filename) + " failed");

  fclose(fd);
}

template<typename T>
void read_fasta_file(const char *filename, std::vector<T>& v){
    FILE* fd;

    if ((fd = fopen(filename, "r")) == nullptr)
        error("open() file " + std::string(filename) + " failed" );

    v.clear();

    char c;
    while (fread( &c, sizeof(char), 1,fd) == 1) {
      if(c == '>'){
        while(fread( &c, sizeof(char), 1,fd) == 1 && c != '\n');
      }else{
        v.push_back(c);
        while(fread( &c, sizeof(char), 1,fd) == 1 && c!= '\n') v.push_back(c);
      }
  	}
  	fclose(fd);
}

template <typename T>
void write_file(const char *filename, std::vector<T> &ptr)
{
  struct stat filestat;
  FILE *fd;

  if ((fd = fopen(filename, "w")) == nullptr)
    error("open() file " + std::string(filename) + " failed");

  size_t length = ptr.size(); 
  if ((fwrite(&ptr[0], sizeof(T), length, fd)) != length)
    error("fwrite() file " + std::string(filename) + " failed");

  fclose(fd);
}

inline bool file_exists(const std::string &name)
{
  struct stat buffer;
  return (stat(name.c_str(), &buffer) == 0);
}
//*********************** Time resources ***************************************

/*!
 * op the operation that we want to measure
 */
#define _elapsed_time(op)                                                                                               \
  ({                                                                                                                    \
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();          \
    op;                                                                                                                 \
    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();            \
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count()); \
    std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count();                                \
  })

//*********************** Kasai et al. LCP construction algorithm ***************************************
template<typename T, typename S, typename lcp_t>
void LCP_array(S* s, const std::vector<T>& isa, const std::vector<T>& sa, size_t n, std::vector<lcp_t>& lcp){
  lcp[0]  = 0;

  T l = 0;
  for (size_t i = 0; i < n; ++i){
    // if i is the last character LCP is not defined
    T k = isa[i];
    if(k > 0){
      T j = sa[k-1];
      // I find the longest common prefix of the i-th suffix and the j-th suffix.
      while(s[i+l] == s[j+l]) l++;
      // l stores the length of the longest common prefix between the i-th suffix and the j-th suffix
      lcp[k] = l;
      if(l>0) l--;
    }
  }
}

//*********************** Kasai et al. LCP construction algorithm rec text ***************************************
template<typename T, typename S, typename lcp_t>
void LCP_array_cyclic_text(S* s, const std::vector<T>& isa, const std::vector<T>& sa, size_t n, std::vector<lcp_t>& lcp){
  lcp[0] = 0;

  T l = 0;
  for (size_t i = 0; i < n; ++i){
    // if i is the last character LCP is not defined
    T k = isa[i];
    if(k > 0){
      T j = sa[k-1];
      // I find the longest common prefix of the i-th suffix and the j-th suffix.
      while(l <= n && s[(i+l) % n] == s[(j+l) % n]) l++;
      // l stores the length of the longest common prefix between the i-th suffix and the j-th suffix
      lcp[k] = l;
      if(l>0) l--;
    }
  }
}



//********** begin my serialize edit from sdsl ********************
// Those are wrapper around most of the serialization functions of sdsl

template <class T, typename size_type>
uint64_t
my_serialize_array(const T* p, const size_type size, std::ostream &out, typename std::enable_if<std::is_fundamental<T>::value>::type * = 0)
{
  size_t written_bytes = 0;
  if (size > 0)
  {

    size_type idx = 0;
    while (idx + sdsl::conf::SDSL_BLOCK_SIZE < (size))
    {
      out.write((char *)p, sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T));
      written_bytes += sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T);
      p += sdsl::conf::SDSL_BLOCK_SIZE;
      idx += sdsl::conf::SDSL_BLOCK_SIZE;
    }
    out.write((char *)p, ((size) - idx) * sizeof(T));
    written_bytes += ((size) - idx) * sizeof(T);

  }
  return written_bytes;
}

//! Serialize each element of an std::vector
/*!
 * \param vec The vector which should be serialized.
 * \param out Output stream to which should be written.
 * \param v   Structure tree node. Note: If all elements have the same
 *            structure, then it is tried to combine all elements (i.e.
 *            make one node w with size set to the cumulative sum of all
 *           sizes of the children)
 */
// specialization for fundamental types
template <class T>
uint64_t
my_serialize_vector(const std::vector<T> &vec, std::ostream &out, sdsl::structure_tree_node *v, std::string name, typename std::enable_if<std::is_fundamental<T>::value>::type * = 0)
{
  if (vec.size() > 0)
  {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, "std::vector<" + sdsl::util::class_name(vec[0]) + ">");
    // size_t written_bytes = 0;

    // const T *p = &vec[0];
    // typename std::vector<T>::size_type idx = 0;
    // while (idx + sdsl::conf::SDSL_BLOCK_SIZE < (vec.size()))
    // {
    //   out.write((char *)p, sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T));
    //   written_bytes += sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T);
    //   p += sdsl::conf::SDSL_BLOCK_SIZE;
    //   idx += sdsl::conf::SDSL_BLOCK_SIZE;
    // }
    // out.write((char *)p, ((vec.size()) - idx) * sizeof(T));
    // written_bytes += ((vec.size()) - idx) * sizeof(T);

    size_t written_bytes = my_serialize_array<T, typename std::vector<T>::size_type>(&vec[0], vec.size(), out);

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }
  else
  {
    return 0;
  }
}

template <typename X>
uint64_t
my_serialize(const std::vector<X> &x,
             std::ostream &out, sdsl::structure_tree_node *v = nullptr,
             std::string name = "", typename std::enable_if<std::is_fundamental<X>::value>::type * = 0)
{
  return sdsl::serialize(x.size(), out, v, name) + my_serialize_vector(x, out, v, name);
}

/**
 * @brief Load an array of size elements into p. p should be preallocated.
 * 
 * \tparam T 
 * \tparam size_type 
 * @param p 
 * @param size 
 * @param in 
 */
template <class T, typename size_type>
void my_load_array(T *p, const size_type size, std::istream &in, typename std::enable_if<std::is_fundamental<T>::value>::type * = 0)
{
  size_type idx = 0;
  while (idx + sdsl::conf::SDSL_BLOCK_SIZE < (size))
  {
    in.read((char *)p, sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T));
    p += sdsl::conf::SDSL_BLOCK_SIZE;
    idx += sdsl::conf::SDSL_BLOCK_SIZE;
  }
  in.read((char *)p, ((size) - idx) * sizeof(T));
}
//! Load all elements of a vector from a input stream
/*! \param vec  Vector whose elements should be loaded.
 *  \param in   Input stream.
 *  \par Note
 *   The vector has to be resized prior the loading
 *   of its elements.
 */
template <class T>
void my_load_vector(std::vector<T> &vec, std::istream &in, typename std::enable_if<std::is_fundamental<T>::value>::type * = 0)
{
  // T *p = &vec[0];
  // typename std::vector<T>::size_type idx = 0;
  // while (idx + sdsl::conf::SDSL_BLOCK_SIZE < (vec.size()))
  // {
  //   in.read((char *)p, sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T));
  //   p += sdsl::conf::SDSL_BLOCK_SIZE;
  //   idx += sdsl::conf::SDSL_BLOCK_SIZE;
  // }
  // in.read((char *)p, ((vec.size()) - idx) * sizeof(T));
  my_load_array<T, typename std::vector<T>::size_type>(&vec[0], vec.size(), in);
}

template <typename X>
void my_load(std::vector<X> &x, std::istream &in, typename std::enable_if<std::is_fundamental<X>::value>::type * = 0)
{
  typename std::vector<X>::size_type size;
  sdsl::load(size, in);
  x.resize(size);
  my_load_vector(x, in);
}

//*********************** Timing *********************************************

#ifdef MTIME

#ifndef _MTIME
#define _MTIME

#define MTIME_TSAFE_INIT(_n)                                                    \
  std::mutex __mtime_mtx;                                                       \
  std::vector<double> __tsafe_durations(_n, 0.0);                               

#define MTIME_TSAFE_NAMES_INIT( args... )                                       \
  std::vector<std::string> __tsafe_names = {args};

#define MTIME_INIT(_n)                                                                                                                \
  std::vector<std::chrono::high_resolution_clock::time_point> __watches(_n);                                                          \
  std::vector<double> __durations(_n, 0.0);

#define MTIME_START(_i) \
  __watches[_i] = std::chrono::high_resolution_clock::now();

#define MTIME_END(_i) \
  __durations[_i] += std::chrono::duration<double, std::ratio<1>>(std::chrono::high_resolution_clock::now() - __watches[_i]).count();

#define MTIME_REPORT(_i) \
  verbose("Timing variable: ", _i, " ", __durations[_i]);

#define MTIME_REPORT_ALL                          \
  for (size_t i = 0; i < __durations.size(); ++i) \
  MTIME_REPORT(i)

#define MTIME_TSAFE_REPORT(_i) \
  verbose("Timing variable", _i, ":",std::setw(9), __tsafe_durations[_i], "\t(", __tsafe_names[i], ")");

#define MTIME_TSAFE_REPORT_ALL                          \
  for (size_t i = 0; i < __tsafe_durations.size(); ++i) \
  MTIME_TSAFE_REPORT(i)

#define MTIME_TSAFE_MERGE                         \
  __mtime_mtx.lock();                             \
  for (size_t i = 0; i < __durations.size(); ++i) \
    __tsafe_durations[i] += __durations[i];       \
  __mtime_mtx.unlock();



#endif /* _MTIME */

#else
  #define MTIME_TSAFE_INIT(_n)
  #define MTIME_TSAFE_NAMES_INIT( args... )
  #define MTIME_INIT(_n)
  #define MTIME_START(_i)
  #define MTIME_END(_i)
  #define MTIME_REPORT(_i)
  #define MTIME_REPORT_ALL
  #define MTIME_TSAFE_REPORT(_i)
  #define MTIME_TSAFE_REPORT_ALL
  #define MTIME_TSAFE_MERGE
#endif




//***********************  Utils ***********************************************

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

#define maxl(a,b) (a) = std::max((a),(b))
#define minl(a,b) (a) = std::min((a),(b))
#define dist(a,b) ( (a) > (b)? ((a) - (b)) : ((b) - (a)))

////////////////////////////////////////////////////////////////////////////////
/// helper functions
////////////////////////////////////////////////////////////////////////////////

static inline char complement(const char n)
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


//*********************** Print Utils ******************************************


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

#endif /* end of include guard: _COMMON_HH */
