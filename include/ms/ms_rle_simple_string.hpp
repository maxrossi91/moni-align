/* ms_rle_simple_string - Extension of the r-index rle_string to compute matching statistics
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
   \file ms_rle_simple_string.hpp
   \brief ms_rle_simple_string.hpp Extension of the r-index rle_string to compute matching statistics.
   \author Massimiliano Rossi
   \date 20/07/2022
*/

#ifndef _MS_RLE_SIMPLE_STRING_HH
#define _MS_RLE_SIMPLE_STRING_HH

#include <common.hpp>

#include <definitions.hpp>
#include <huff_string.hpp>
#include <sparse_sd_vector.hpp>

// Semplification of the rle_string.hpp of the r-index

template <
    class sparse_bitvector_t = ri::sparse_sd_vector, //predecessor structure storing run length
    class string_t = ri::huff_string                 //run heads
    >
class ms_rle_simple_string
{
public:
    ms_rle_simple_string()
    {
        //NtD
    }

    /*
     * constructor: build structure on the input string
     * \param input the input string without 0x0 bytes in it.
     * \param B block size. The main sparse bitvector has R/B bits set (R being number of runs)
     *
     */
    ms_rle_simple_string(const std::string &input) 
    {
        build_from_string(std::istringstream(input));
    }

    ms_rle_simple_string(std::ifstream &ifs)
    {
        build_from_string(ifs);
    }

    // Construction from run-length encoded BWT
    ms_rle_simple_string(std::ifstream &heads, std::ifstream &lengths)
    {
        build_from_rl_string(heads,lengths);
    }

    inline ri::uchar operator[](ri::ulint i){
        assert(i<n);
        return run_heads[run_of_position(i)];
    }

    inline ri::ulint size() const {return n;}
    /*
     * position of i-th character c. i starts from 0!
     */
    ri::ulint select(ri::ulint i, ri::uchar c) const{
        assert(i<runs_per_letter[c].size());
        //i-th c is inside j-th c-run (j starts from 0)
        assert(i<runs_per_letter[c].size());
        ri::ulint j = runs_per_letter[c].rank(i);
        //starting position of i-th c inside its run
        assert(j==0 || i >= runs_per_letter[c].select(j-1) + 1);
        ri::ulint before = (j==0 ? i : i - (runs_per_letter[c].select(j-1) + 1));
        //position in run_heads
        ri::ulint r = run_heads.select(j,c);
        ri::ulint k = (r>0)?runs.select(r-1)+1:0;
        return k + before;
    }

    /*
     * number of c before position i
     */
    ri::ulint rank(ri::ulint i, ri::uchar c) const{
        assert(i<=n);
        //letter does not exist in the text
        const size_t n_c = number_of_letter(c);
        if(n_c == 0) return 0;
        if(i == n) return n_c;

        ri::ulint current_run = run_of_position(i);
        assert(current_run<R);
        ri::ulint start = (current_run>0)?runs.select(current_run-1)+1:0;
        ri::ulint dist = i - start;
        assert(dist < n);
        //number of c runs before the current run
        ri::ulint rk = run_heads.rank(current_run,c);
        //number of c before i in the current run
        ri::ulint tail = (run_heads[current_run]==c) ? dist : 0;
        //in this case, either there are no c before position i
        //or the current run is the first of the kind ccc...cc
        if(rk==0) return tail;
        return runs_per_letter[c].select(rk-1)+1+tail;
    }

    /*
    * return inclusive range of j-th run in the string
    */
    pair<ri::ulint,ri::ulint> run_range(ri::ulint j) const{
        assert(j<run_heads.size());
        ri::ulint pos_start = (j>0)?runs.select(j-1)+1:0;
        ri::ulint pos_end = runs.select(j);
        return {pos_start,pos_end};
    }

    //break range: given a range <l',r'> on the string and a character c, this function
    //breaks <l',r'> in maximal sub-ranges containing character c.
    //for simplicity and efficiency, we assume that characters at range extremities are both 'c'
    //thanks to the encoding (run-length), this function is quite efficient: O(|result|) ranks and selects
    vector<ri::range_t> break_range(ri::range_t rn,ri::uchar c) const{
        auto l = rn.first;
        auto r = rn.second;
        assert(l<=r);
        assert(r<size());
        assert(operator[](l)==c);
        assert(operator[](r)==c);
        //retrieve runs that contain positions l and r
        auto run_l = run_of(l);
        auto run_r = run_of(r);
        //in this case rn contains only character c: do not break
        if(run_l.first==run_r.first) return {rn};
        vector<ri::range_t> result;
        //first range: from l to the end of the run containing position l
        result.push_back({l,run_l.second});
        //rank of c's of interest in run_heads
        ri::ulint rank_l = run_heads.rank(run_l.first,c);
        ri::ulint rank_r = run_heads.rank(run_r.first,c);
        //now retrieve run bounds of all c-runs of interest
        for(ri::ulint j = rank_l+1;j<rank_r;++j){
            result.push_back(run_range(run_heads.select(j,c)));
        }
        //now last (possibly incomplete) run
        auto range = run_range(run_heads.select(rank_r,c));
        result.push_back({range.first,r});
        return result;
    }
    /*
     * input: inclusive range rn, character c
     *
     * return the position j that is closest to rn.first,
     * such that character in position j is c and that is
     * adjacent to a position j' inside rn that contains a
     * character != c
     *
     * rn must contain c and at least another character d!=c
     *
     */
    ri::ulint closest_run_break(ri::range_t rn, ri::uchar c){
        /*
         * case 1: range begins with a c-run: return last position of the run
         */
        if(operator[](rn.first)==c){
            ri::ulint i = run_of_position(rn.first);
            ri::ulint j = run_range(i).second;
            //j must be inside rn, i.e. rn must not contain only c
            //j must not be last position of rn: this would imply
            //that rn contain only c
            assert(j<rn.second);
            return j;
        }else{
            //case 2: first c-run starts in the middle of the range
            //rank i of first c in the range
            ri::ulint i = rank(rn.first,c);
            assert(i<rank(size(),c));
            //map from rank space to string position:
            //i now is the first position inside the range that contains c
            i = select(i,c);
            assert(operator[](i)==c);
            assert(i<=rn.second);
            return i;
        }
    }

    inline size_t number_of_runs_of_letter(const uint8_t c) const
    {
        return runs_per_letter[c].number_of_1();
    }

    inline size_t number_of_letter(uint8_t c) const
    {
        return runs_per_letter[c].size();
    }

    // i-th run head
    inline uint8_t head_of(const size_t i) const
    {
        assert(i<R);
        return run_heads[i];
    }

    // rank in chracters of the i-th run head
    // i.e., the number of characters c before the first character of the run.
    inline size_t head_rank(const size_t i, const uint8_t c) const
    {
        assert(i < R);
        size_t j = run_heads.rank(i, c);
        if(j < 1) return j;
        assert(j<=i);
        return runs_per_letter[c].select(j-1) + 1; // j-1 because the select is 0 based
    }
    // number of runs of character c in in position i
    inline size_t run_head_rank(const size_t i, const uint8_t c) const
    {
        assert(i < R);
        return run_heads.rank(i, c);
    }

    inline std::pair<size_t,size_t> run_and_head_rank(const size_t i, const uint8_t c) const
    {
        assert(i < R);
        const size_t j = run_heads.rank(i, c);
        if( j < 1) return std::make_pair(j,j);
        const size_t k = runs_per_letter[c].select(j - 1) + 1; // j-1 because the select is 0 based
        return std::make_pair(j, k);
    }

    // Select the i-th run of c
    inline size_t run_head_select(const size_t i, const uint8_t c) const
    {
        assert(i < R and i > 0);
        return run_heads.select(i - 1, c);
    }

    /*
     * text position i is inside this run
     */
    ri::ulint run_of_position(ri::ulint i) const{
        assert(i<n);
        ri::ulint current_run = runs.rank(i);
        assert(current_run<R);
        return current_run;
    }


    //length of i-th run
    inline ri::ulint run_at(ri::ulint i) {
        assert(i<R);
        ri::uchar c = run_heads[i];
        return runs_per_letter[c].gapAt(run_heads.rank(i,c));
    }

    inline ri::ulint number_of_runs() const {return R;}

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    ri::ulint serialize(std::ostream &out)
    {
        ri::ulint w_bytes = 0;
        out.write((char*)&n,sizeof(n));
        out.write((char*)&R,sizeof(R));
        w_bytes += sizeof(n) + sizeof(R);
        if(n==0) return w_bytes;
        w_bytes += runs.serialize(out);
        for(ri::ulint i=0;i<256;++i)
            w_bytes += runs_per_letter[i].serialize(out);
        w_bytes += run_heads.serialize(out);
        return w_bytes;
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in)
    {
        in.read((char*)&n,sizeof(n));
        in.read((char*)&R,sizeof(R));
        if(n==0) return;
        runs.load(in);
        runs_per_letter = vector<sparse_bitvector_t>(256);
        for(ri::ulint i=0;i<256;++i)
            runs_per_letter[i].load(in);
        run_heads.load(in);
    }

    std::string get_file_extension() const
    {
        return ".rle_simple";
    }
protected:
    void build_from_string(std::istream& is)
    {
        is.clear();
        is.seekg(0);
        // assert(not contains0(input)); // We're hacking the 0 away :)

        R = 0;
        auto runs_per_letter_bv = vector<vector<bool> >(256);
        //runs in main bitvector
        vector<bool> runs_bv;
        string run_heads_s;

        int last_c, c;
        last_c = is.get();
        last_c = last_c>TERMINATOR ? last_c : TERMINATOR; // 0->1 :D
        n = 1;
        while((c = is.get()) != EOF) {
            c = c>TERMINATOR ? c : TERMINATOR; // 0->1 :D
            if(c != last_c){
                run_heads_s.push_back(last_c);
                runs_per_letter_bv[last_c].push_back(true);
                last_c = c;
                runs_bv.push_back(true);
                R++;
            } else{
                runs_bv.push_back(false);
                runs_per_letter_bv[last_c].push_back(false);
            }
            n += 1;
        }
        run_heads_s.push_back(last_c);
        runs_per_letter_bv[last_c].push_back(true);
        runs_bv.push_back(false);
        R++;
        assert(run_heads_s.size()==R);

        //now compact structures
        assert(runs_bv.size()==n);
        ri::ulint t = 0;
        for(ri::ulint i=0;i<256;++i)
            t += runs_per_letter_bv[i].size();
        assert(t==n);
        runs = sparse_bitvector_t(runs_bv);
        //a fast direct array: char -> bitvector.
        runs_per_letter = vector<sparse_bitvector_t>(256);
        for(ri::ulint i=0;i<256;++i)
            runs_per_letter[i] = sparse_bitvector_t(runs_per_letter_bv[i]);
        run_heads = string_t(run_heads_s);
        assert(run_heads.size()==R);
    }

    void build_from_rl_string(std::istream &heads, std::istream &lengths) 
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);

        auto runs_per_letter_bv = vector<vector<bool>>(256);
        //runs in main bitvector
        vector<bool> runs_bv;

        // Reads the run heads
        string run_heads_s;
        heads.seekg(0, heads.end);
        run_heads_s.resize(heads.tellg());
        heads.seekg(0, heads.beg);
        heads.read(&run_heads_s[0], run_heads_s.size());

        size_t pos = 0;
        n = 0;
        R = run_heads_s.size();
        // Compute runs_bv and runs_per_letter_bv
        for (size_t i = 0; i < run_heads_s.size(); ++i)
        {
            size_t length = 0;
            lengths.read((char *)&length, 5);
            if (run_heads_s[i] <= TERMINATOR) // change 0 to 1
                run_heads_s[i] = TERMINATOR;

            std::fill_n(std::back_inserter(runs_bv), length - 1, false);
            runs_bv.push_back(true);

            std::fill_n(std::back_inserter(runs_per_letter_bv[run_heads_s[i]]), length - 1, false);
            runs_per_letter_bv[run_heads_s[i]].push_back(true);

            n += length;
        }

        //now compact structures
        assert(runs_bv.size() == n);
        ri::ulint t = 0;
        for (ri::ulint i = 0; i < 256; ++i)
            t += runs_per_letter_bv[i].size();
        assert(t == n);
        runs = sparse_bitvector_t(runs_bv);
        //a fast direct array: char -> bitvector.
        runs_per_letter = vector<sparse_bitvector_t>(256);
        for (ri::ulint i = 0; i < 256; ++i)
            runs_per_letter[i] = sparse_bitvector_t(runs_per_letter_bv[i]);
        run_heads = string_t(run_heads_s);
        assert(run_heads.size() == R);
    }

    static ri::ulint count_runs(string &s){
        ri::ulint runs=1;
        for(ri::ulint i=1;i<s.size();++i){
            if(s[i]!=s[i-1]) runs++;
        }
        return runs;
    }

    //<j=run of position i, last position of j-th run>
    pair<ri::ulint,ri::ulint> run_of(ri::ulint i){
        ri::ulint current_run = run_of_position(i);
        ri::ulint pos = runs.select(current_run);
        assert(pos>0);
        assert(current_run<R);
        return {current_run,pos};
    }

    bool contains0(string &s){
        for(auto c : s)
            if(c==0) return true;
        return false;
    }

    sparse_bitvector_t runs;
    //for each letter, its runs stored contiguously
    vector<sparse_bitvector_t> runs_per_letter;
    //store run heads in a compressed string supporting access/rank
    string_t run_heads;
    //text length and number of runs
    ri::ulint n=0;
    ri::ulint R=0;
    static const ri::uchar TERMINATOR=1;
};

// Construction from run-length encoded BWT specialization for sparse_sd_vector
template <>
ms_rle_simple_string<ri::sparse_sd_vector, ri::huff_string>::ms_rle_simple_string(std::ifstream &heads, std::ifstream &lengths)
{
    heads.clear();
    heads.seekg(0);
    lengths.clear();
    lengths.seekg(0);

    // Reads the run heads
    string run_heads_s;
    heads.seekg(0, heads.end);
    run_heads_s.resize(heads.tellg());
    heads.seekg(0, heads.beg);
    heads.read(&run_heads_s[0], run_heads_s.size());

    size_t pos = 0;
    n = 0;
    R = run_heads_s.size();

    auto runs_per_letter_bv = vector<vector<size_t>> (256);
    auto runs_per_letter_bv_i = vector<size_t> (256,0);
    //runs in main bitvector
    vector<size_t> runs_bv_onset;
    size_t runs_bv_i = 0;
    // Compute runs_bv and runs_per_letter_bv
    for (size_t i = 0; i < run_heads_s.size(); ++i)
    {
        size_t length = 0;
        lengths.read((char *)&length, 5);
        if (run_heads_s[i] <= TERMINATOR) // change 0 to 1
            run_heads_s[i] = TERMINATOR;

        runs_bv_onset.push_back(n + length - 1);

        assert(length > 0);
        runs_per_letter_bv_i[run_heads_s[i]] += length;
        runs_per_letter_bv[run_heads_s[i]].push_back(runs_per_letter_bv_i[run_heads_s[i]] - 1);

        n += length;
    }

    //now compact structures
    ri::ulint t = 0;
    for (ri::ulint i = 0; i < 256; ++i)
        t += runs_per_letter_bv_i[i];
    assert(t == n);
    runs = ri::sparse_sd_vector(runs_bv_onset, n);
    //a fast direct array: char -> bitvector.
    runs_per_letter = vector<ri::sparse_sd_vector>(256);
    for (ri::ulint i = 0; i < 256; ++i)
        runs_per_letter[i] = ri::sparse_sd_vector(runs_per_letter_bv[i],runs_per_letter_bv_i[i]);
    run_heads = ri::huff_string(run_heads_s);
    assert(run_heads.size() == R);
};

typedef ms_rle_simple_string<ri::sparse_sd_vector> ms_rle_simple_string_sd;
typedef ms_rle_simple_string<ri::sparse_hyb_vector> ms_rle_simple_string_hyb;

#endif /* end of include guard: _MS_RLE_SIMPLE_STRING_HH */
