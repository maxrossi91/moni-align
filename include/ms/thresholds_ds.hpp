/* thresholds_ds - Stores the thresholds in compressed and plain ways 
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
   \file thresholds_ds.hpp
   \brief thresholds_ds.hpp Stores the thresholds in compressed and plain ways.
   \author Massimiliano Rossi
   \date 09/07/2020
*/

#ifndef _MS_THRESHOLDS_DS_HH
#define _MS_THRESHOLDS_DS_HH

#include <common.hpp>

#include <malloc_count.h>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include <ms_rle_string.hpp>

template <class rle_string_t = ms_rle_string_sd>
class thr_plain
{
public:
    int_vector<> thresholds;
    rle_string_t &bwt;

    typedef size_t size_type;

    thr_plain() {}

    thr_plain(std::string filename, rle_string_t& bwt_):bwt(bwt_)
    {
        int log_n = bitsize(uint64_t(bwt_.size()));

        verbose("Reading thresholds from file");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string tmp_filename = filename + std::string(".thr_pos");

        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
            error("open() file " + tmp_filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + tmp_filename + " failed");

        if (filestat.st_size % THRBYTES != 0)
            error("invilid file " + tmp_filename);

        size_t length = filestat.st_size / THRBYTES;
        size_t threshold = 0;

        thresholds = int_vector<>(length, 0, log_n);

        for (size_t i = 0; i < length; ++i)
        {
            size_t threshold = 0;
            if ((fread(&threshold, THRBYTES, 1, fd)) != 1)
                error("fread() file " + tmp_filename + " failed");
            thresholds[i] = threshold;
        }

        fclose(fd);

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    }

    size_t operator[] (size_t& i)
    {
        assert( i < thresholds.size());
        return thresholds[i];
    }
};

template <class rle_string_t = ms_rle_string_sd>
class thr_compressed
{
public:
    int_vector<> thresholds;
    rle_string_t &bwt;
    long long min_off;

    typedef size_t size_type;

    thr_compressed() {}

    thr_compressed(std::string filename, rle_string_t& bwt_):bwt(bwt_)
    {
        int log_n = bitsize(uint64_t(bwt_.size()));
        size_t n = uint64_t(bwt.size());

        verbose("Reading thresholds from file");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string tmp_filename = filename + std::string(".thr_pos");

        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
            error("open() file " + tmp_filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + tmp_filename + " failed");

        if (filestat.st_size % THRBYTES != 0)
            error("invilid file " + tmp_filename);

        size_t length = filestat.st_size / THRBYTES;

        size_t pos = 0;

        long long max_off = 0;
        min_off = n;

        for (size_t i = 0; i < length; ++i)
        {
            size_t threshold = 0;
            if ((fread(&threshold, THRBYTES, 1, fd)) != 1)
                error("fread() file " + tmp_filename + " failed");

            long long off = 0;

            if (threshold > 0)
            {
                uint8_t c = bwt.head_of(i);
                size_t pred = bwt.select(bwt.rank(pos - 1, c) - 1, c);
                size_t mid_int = (pos - pred + 1) >> 1;
                assert(threshold > pred);

                threshold = threshold - pred;

                off = mid_int - threshold;

                max_off = max(max_off, off);
                min_off = min(min_off, off);
            }

            pos += bwt.run_at(i);
        }

        // Rewind the file
        fseek(fd,0, SEEK_SET);
        pos = 0;
        
        int log_off = bitsize((size_t)(max_off - min_off + 1));

        min_off = -min_off; // Shift all the values
        thresholds = int_vector<>(length,0,log_off);
        for (size_t i = 0; i < length; ++i)
        {
            size_t threshold = 0;
            if ((fread(&threshold, THRBYTES, 1, fd)) != 1)
                error("fread() file " + tmp_filename + " failed");

            long long off = 0;

            if (threshold > 0)
            {
                uint8_t c = bwt.head_of(i);
                size_t pred = bwt.select(bwt.rank(pos - 1, c) - 1, c);
                size_t mid_int = (pos - pred + 1) >> 1;
                assert(threshold > pred);

                threshold = threshold - pred;

                off = mid_int - threshold + min_off;
            }

            thresholds[i] = off;
            pos += bwt.run_at(i);
        }

        fclose(fd);



        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    }

    size_t operator[] (size_t& i)
    {
        assert( i < thresholds.size());

        // get mid_interval
        uint8_t c = bwt.head_of(i);
        size_t rank = bwt.head_rank(i, c);
        if(rank == 0)
            return 0;

        size_t pred = bwt.select(rank-1, c);
        size_t pos = bwt.select(rank, c);
        size_t mid_int = (pos - pred + 1) >> 1;

        return mid_int + min_off - thresholds[i];
    }
};

#endif /* end of include guard: _MS_THRESHOLDS_DS_HH */
