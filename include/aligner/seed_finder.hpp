/* seed_finder - MEMs types
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
   \file seed_finder.hpp
   \brief seed_finder.hpp MEMs type
   \author Massimiliano Rossi
   \date 19/10/2021
*/

#ifndef _SEED_FINDER_HH
#define _SEED_FINDER_HH


#define VERBOSE
#define MTIME

#include <vector>
#include <common.hpp>


#include <mems.hpp>
#include <moni.hpp>
#include <moni_lcp.hpp>

#include <malloc_count.h>

#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>
#include <PlainSlp.hpp>
#include <FixedBitLenCode.hpp>

#include <slp_definitions.hpp>

template <typename slp_t,
          typename ms_t>
class seed_finder
{
protected:
    size_t n;
    bool filter_seeds = true;
    size_t n_seeds_thr = 5000;
    const size_t min_len = 0;

public:
    ms_t ms;
    slp_t ra;

    typedef slp_t slp_type;

    seed_finder(std::string filename, size_t min_len_, bool filter_seeds_, size_t n_seeds_thr_):
        min_len(min_len_),
        filter_seeds(filter_seeds_),
        n_seeds_thr(n_seeds_thr_)
    {
        verbose("Loading the matching statistics index");
        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string filename_ms = filename + ms.get_file_extension();

        if (not file_exists(filename_ms))
        error("File not found: ", filename_ms);

        ifstream fs_ms(filename_ms);
        ms.load(fs_ms);
        fs_ms.close();

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Matching statistics index loading complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

        t_insert_start = std::chrono::high_resolution_clock::now();

        std::string filename_slp = filename + get_slp_file_extension<slp_t>();
        verbose("Loading random access file: " + filename_slp);

        if (not file_exists(filename_slp))
        error("File not found: ", filename_slp);

        ifstream fs(filename_slp);
        ra.load(fs);
        fs.close();

        n = ra.getLen();
        t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Random access loading complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

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
        size_t prev_pos_plus_one = n + 1;
        for (size_t i = 0; i < pointers.size(); ++i)
        {
            size_t pos = pointers[i];
            while (pos != prev_pos_plus_one && (i + l) < read->seq.l && (pos + l) < n && read->seq.s[i + l] == ra.charAt(pos + l))
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
                // find_MEM_occs(mems.back());
            }

            // Compute next match length
            pl = l;
            l = (l == 0 ? 0 : (l - 1));

            prev_pos_plus_one = pos + 1;
        }

    }

    // Fill the vector of occs with the matches above curr with LCP <= len
    bool find_MEM_above(size_t curr, size_t len, std::vector<size_t>& occs)
    {
        assert (len > 0);
        auto [prev, lcp] = get_prev_occ_with_lcp(curr, len);
        
        while( lcp >= len )
        {
            occs.push_back(prev);
            std::tie(prev,lcp) = get_prev_occ_with_lcp(prev, len);

            if (filter_seeds and (occs.size() > n_seeds_thr) )
                return false;
        }

        return true;
    }


    // Fill the vector of occs with the matches below curr with LCP <= len
    bool find_MEM_below(size_t curr, size_t len, std::vector<size_t>& occs)
    {
        assert(len > 0);
        auto [next, lcp] = get_next_occ_with_lcp(curr, len);
        
        while( lcp >= len )
        {
            occs.push_back(next);
            std::tie(next,lcp) = get_next_occ_with_lcp(next, len);

            if (filter_seeds and (occs.size() > n_seeds_thr) )
                return false;
        }

        return true;
    }



    // Fill the vector of occurrences of the mem_t data structure
    bool find_MEM_occs(mem_t &mem)
    {
        mem.occs.push_back(mem.pos);

        if (!find_MEM_above(mem.pos, mem.len, mem.occs)) return false;
        if (!find_MEM_below(mem.pos, mem.len, mem.occs)) return false;

        return true;
    }


    // Populate the seeds given a list of MEMs
    bool populate_seed(
        mem_t &mem, std::vector<mem_t> &mems)
    {
        size_t l = mem.len;
        size_t i = mem.idx;
        size_t mate = mem.mate;
        size_t pos = mem.pos;
        size_t r = mem.rpos;

        // size_t r = r_offset + (i + l - 1); // compatible with minimap2 chaining algorithm

        mem.occs.push_back(mem.pos);

        find_MEM_above(mem.pos, mem.len, mem.occs);
        size_t upper_suffix = mem.occs.back();
        find_MEM_below(mem.pos, mem.len, mem.occs);
        size_t lower_suffix = mem.occs.back();

        // Take two halves of the MEM
        if (l >= (min_len << 1))
        {
            size_t ll = l >> 1;
            size_t rl = r - l + ll; // compatible with minimap2 chaining algorithm
            mems.push_back(mem_t(upper_suffix, ll, i, mate, rl));

            mem_t &mem = mems.back();
            if ((not find_MEM_above(upper_suffix, mem.len, mem.occs)) or
                (not find_MEM_below(lower_suffix, mem.len, mem.occs)))
            {
                mems.pop_back();
                return false;
            }

            size_t lr = l - ll;
            size_t rr = r; // compatible with minimap2 chaining algorithm
            mems.push_back(mem_t(pos + ll, lr, i + ll, mate, rr));
            // find_MEM_occs(mems.back()); // TODO: Optimize this
            if ((not find_MEM_occs(mems.back())))
            {
                mems.pop_back();
                return false;
            }
        }

        return true;
    }

    // // Populate the seeds given a list of MEMs
    void populate_seeds(
        std::vector<mem_t> &mems)
    {
        size_t n_MEMs = mems.size();
        for (size_t j = 0; j < n_MEMs; ++j)
            populate_seed(mems[j], mems);
    }


    void find_seeds(
        const kseq_t *read,
        std::vector<mem_t> &mems,
        size_t r_offset = 0,
        size_t mate = 0)
    {
        MMEM_START(8);
        find_mems(read, mems, r_offset, mate);
        MMEM_END(8);
        MMEM_START(9);
        populate_seeds(mems);
        MMEM_END(9);
    }

protected:
    inline std::pair<size_t, size_t> get_next_occ_with_lcp(size_t curr, size_t len) 
    {
        if (curr == ms.get_last_run_sample()) return {ms.get_first_run_sample(), 0};

        size_t next = ms.Phi_inv(curr);
        size_t lcp = 0; // lceToR(ra,curr,next);

        if((n-curr) >= len and (n-next) >= len)
            lcp  = lceToRBounded(ra,curr,next,len);

        return {next, lcp};
    }

    inline std::pair<size_t, size_t> get_prev_occ_with_lcp(size_t curr, size_t len)  
    {
        if (curr == ms.get_first_run_sample()) return {ms.get_last_run_sample(), 0};

        size_t prev = ms.Phi(curr);
        size_t lcp = 0; // lceToR(ra,curr,prev);

        if((n-curr) >= len and (n-prev) >= len)
            lcp  = lceToRBounded(ra,curr,prev,len);

        return {prev, lcp};
    }


};


// Computes the prev and next occ with moni_lcp
template<>
// template <typename slp_t>
inline std::pair<size_t, size_t> seed_finder<plain_slp_t, moni_lcp<>>::get_next_occ_with_lcp(size_t curr, size_t len)
{
    if (curr == ms.get_last_run_sample()) return {ms.get_first_run_sample(), 0};

    return ms.Phi_inv_lcp(curr);
}

template<>
// template <typename slp_t>
inline std::pair<size_t, size_t> seed_finder<plain_slp_t, moni_lcp<>>::get_prev_occ_with_lcp(size_t curr, size_t len)
{
    if (curr == ms.get_first_run_sample()) return {ms.get_last_run_sample(), 0};

    return ms.Phi_lcp(curr);
}



#endif /* end of include guard: _SEED_FINDER_HH */