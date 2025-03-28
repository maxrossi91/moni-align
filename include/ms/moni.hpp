/* ms_pointers - Computes the matching statistics pointers from BWT and Thresholds 
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
   \file ms_pointers.hpp
   \brief ms_pointers.hpp Computes the matching statistics pointers from BWT and Thresholds.
   \author Massimiliano Rossi
   \date 09/07/2020
*/

#ifndef _MS_POINTERS_HH
#define _MS_POINTERS_HH

#include <common.hpp>

#include <malloc_count.h>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include <r_index.hpp>

#include <ms_rle_string.hpp>
#include <thresholds_ds.hpp>

template <class sparse_bv_type = ri::sparse_sd_vector,
          class rle_string_t = ms_rle_string_sd,
          class thresholds_t = thr_bv<rle_string_t> >
class ms_pointers : public ri::r_index<sparse_bv_type, rle_string_t>
{
public:
    thresholds_t thresholds;

    // std::vector<ulint> samples_start;
    sparse_bv_type pred_start;
    int_vector<> samples_start; //text positions corresponding to first characters in BWT runs, in BWT order
    int_vector<> pred_start_to_run; //stores the BWT run (0...R-1) corresponding to each position in pred_start, in text order
    // int_vector<> samples_end;
    // std::vector<ulint> samples_last;

    // static const uchar TERMINATOR = 1;
    // bool sais = true;
    // /*
    //  * sparse RLBWT: r (log sigma + (1+epsilon) * log (n/r)) (1+o(1)) bits
    //  */
    // //F column of the BWT (vector of 256 elements)
    // std::vector<ulint> F;
    // //L column of the BWT, run-length compressed
    // rle_string_t bwt;
    // ulint terminator_position = 0;
    // ulint r = 0; //number of BWT runs

    typedef size_t size_type;

    ms_pointers() {}

    ms_pointers(std::string filename, bool rle = false) : ri::r_index<sparse_bv_type, rle_string_t>()
    {
        verbose("Building the r-index from BWT");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string bwt_fname = filename + ".bwt";

        verbose("RLE encoding BWT and computing SA samples");

        if (rle)
        {
            std::string bwt_heads_fname = bwt_fname + ".heads";
            std::ifstream ifs_heads(bwt_heads_fname);
            std::string bwt_len_fname = bwt_fname + ".len";
            std::ifstream ifs_len(bwt_len_fname);
            this->bwt = rle_string_t(ifs_heads, ifs_len);

            ifs_heads.seekg(0);
            ifs_len.seekg(0);
            this->build_F_(ifs_heads, ifs_len);
        }
        else
        {
            std::ifstream ifs(bwt_fname);
            this->bwt = rle_string_t(ifs);

            ifs.seekg(0);
            this->build_F(ifs);
        }
        // std::string istring;
        // read_file(bwt_fname.c_str(), istring);
        // for(size_t i = 0; i < istring.size(); ++i)
        //     if(istring[i]==0)
        //         istring[i] = TERMINATOR;
        // this->bwt = rle_string_t(istring);

        this->r = this->bwt.number_of_runs();
        ri::ulint n = this->bwt.size();
        int log_r = bitsize(uint64_t(this->r));
        int log_n = bitsize(uint64_t(this->bwt.size()));

        verbose("Number of BWT equal-letter runs: r = ", this->r);
        verbose("Rate n/r = ", double(this->bwt.size()) / this->r);
        verbose("log2(r) = ", log2(double(this->r)));
        verbose("log2(n/r) = ", log2(double(this->bwt.size()) / this->r));

        // this->build_F(istring);
        // istring.clear();
        // istring.shrink_to_fit();

        verbose("Building phi");
        build_phi(filename + ".ssa", this->r, n, this->pred, this->pred_to_run);

        verbose("Building inverse of phi");
        build_phi(filename + ".esa", this->r, n, pred_start, pred_start_to_run);

        read_samples(filename + ".ssa", this->r, n, samples_start);
        read_samples(filename + ".esa", this->r, n, this->samples_last);

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("R-index construction complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

        verbose("Reading thresholds from file");

        t_insert_start = std::chrono::high_resolution_clock::now();

        thresholds = thresholds_t(filename,&this->bwt);

        t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    }

    void read_samples(std::string filename, ulint r, ulint n, int_vector<> &samples)
    {
        int log_n = bitsize(uint64_t(n));
        
        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(filename.c_str(), "r")) == nullptr)
            error("open() file " + filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + filename + " failed");

        if (filestat.st_size % SSABYTES != 0)
            error("invilid file " + filename);

        size_t length = filestat.st_size / (2 * SSABYTES);
        //Check that the length of the file is 2*r elements of 5 bytes
        assert(length == r);

        // Create the vector
        samples = int_vector<>(r, 0, log_n);

        // Read the vector
        uint64_t left = 0;
        uint64_t right = 0;
        size_t i = 0;
        while (fread((char *)&left, SSABYTES, 1, fd) && fread((char *)&right, SSABYTES, 1, fd))
        {
            ulint val = (right ? right - 1 : n - 1);
            assert(bitsize(uint64_t(val)) <= log_n);
            samples[i++] = val;
        }

        fclose(fd);
    }

    void build_phi(std::string filename, ulint r, ulint n, sparse_bv_type& pv, int_vector<> &pv_to_run)
    {
        int log_r = bitsize(uint64_t(r));
        int log_n = bitsize(uint64_t(n));
        
        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(filename.c_str(), "r")) == nullptr)
            error("open() file " + filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + filename + " failed");

        if (filestat.st_size % SSABYTES != 0)
            error("invilid file " + filename);

        size_t length = filestat.st_size / (2 * SSABYTES);
        //Check that the length of the file is 2*r elements of 5 bytes
        assert(length == r);

        // Create the vector
        std::vector<std::pair<ulint,ulint>> samples = std::vector<std::pair<ulint,ulint>>(r);

        // Read the vector
        uint64_t left = 0;
        uint64_t right = 0;
        size_t i = 0;
        while (fread((char *)&left, SSABYTES, 1, fd) && fread((char *)&right, SSABYTES, 1, fd))
        {
            ulint val = (right ? right - 1 : n - 1);
            assert(bitsize(uint64_t(val)) <= log_n);
            samples[i] = std::make_pair(val,i);
            i++;
        }

        // Sort the sa samples
        std::sort(samples.begin(),samples.end());
        //build Elias-Fano predecessor
        {
            size_t i = 0;
            std::vector<size_t> pred_bv_onset(r);
            for(auto p : samples){
                assert(p.first < this->bwt.size());
                pred_bv_onset[i++] = p.first;
            }
            pv = sparse_bv_type(pred_bv_onset, this->bwt.size());
        }
        assert(pv.rank(pv.size()) == r);
        // //last text position must be sampled
        // assert(pv[pv.size()-1]);
        pv_to_run = int_vector<>(r,0,log_r); //stores the BWT run (0...R-1) corresponding to each position in pred, in text order
        for(ulint i=0;i<samples.size();++i){
            assert(bitsize(uint64_t(samples[i].second)) <= log_r);
            pv_to_run[i] = samples[i].second;
        }

        //  sdsl::nullstream ns;

        // verbose("Memory consumption (bytes).");
        // verbose("            pv: ", pv.serialize(ns));
        // verbose("     pv_to_run: ", pv_to_run.serialize(ns));
        
        fclose(fd);
    }

    vector<ulint> build_F_(std::ifstream &heads, std::ifstream &lengths)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);

        this->F = vector<ulint>(256, 0);
        int c;
        ulint i = 0;
        while ((c = heads.get()) != EOF)
        {
            size_t length = 0;
            lengths.read((char *)&length, 5);
            if (c > TERMINATOR)
                this->F[c] += length;
            else
            {
                this->F[TERMINATOR] += length;
                this->terminator_position = i;
            }
            i++;
        }
        for (ulint i = 255; i > 0; --i)
            this->F[i] = this->F[i - 1];
        this->F[0] = 0;
        for (ulint i = 1; i < 256; ++i)
            this->F[i] += this->F[i - 1];
        return this->F;
    }

    // Computes the matching statistics pointers for the given pattern
    std::vector<size_t> query(const std::vector<uint8_t> &pattern)
    {
        size_t m = pattern.size();

        return _query(pattern.data(), m);
    }

    std::vector<size_t> query(const char* pattern, const size_t m)
    {
        return _query(pattern, m);
    }

    void print_stats()
    {
        sdsl::nullstream ns;

        verbose("Memory consumption (bytes).");
        verbose("   terminator_position: ", sizeof(this->terminator_position));
        verbose("                     F: ", my_serialize(this->F, ns));
        verbose("                   bwt: ", this->bwt.serialize(ns));
        verbose("            thresholds: ", thresholds.serialize(ns));
        verbose("                  pred: ", this->pred.serialize(ns));
        verbose("           pred_to_run: ", this->pred_to_run.serialize(ns));
        verbose("          samples_last: ", this->samples_last.serialize(ns));
        verbose("            pred_start: ", pred_start.serialize(ns));
        verbose("     pred_start_to_run: ", pred_start_to_run.serialize(ns));
        verbose("         samples_start: ", samples_start.serialize(ns));
    }

    /*
     * \param i position in the BWT
     * \param c character
     * \return lexicographic rank of cw in bwt
     */
    ulint LF(ri::ulint i, ri::uchar c)
    {
        // //if character does not appear in the text, return empty pair
        // if ((c == 255 and this->F[c] == this->bwt_size()) || this->F[c] >= this->F[c + 1])
        //     return {1, 0};
        //number of c before the interval
        ri::ulint c_before = this->bwt.rank(i, c);
        // number of c inside the interval rn
        ri::ulint l = this->F[c] + c_before;
        return l;
    }

    ulint get_first_run_sample() {
        return (samples_start[0]+1) % this->bwt.size();
    }


    /*
     * Phi function. Phi_inv(SA[n-1]) is undefined
     */
    ulint Phi_inv(ulint i){
        assert(i != this->bwt.size()-1);
        //jr is the rank of the predecessor of i (circular)
        ulint jr = pred_start.predecessor_rank_circular(i);
        assert(jr<=this->r-1);
        //the actual predecessor
        ulint j = pred_start.select(jr);
        assert(jr<this->r-1 or j == this->bwt.size()-1);
        //distance from predecessor
        ulint delta = j<i ? i-j : i+1;
        //cannot fall on first run: this can happen only if I call Phi(SA[n-1])
        assert(pred_start_to_run[jr] < samples_start.size() - 1);
        //sample at the end of previous run
        assert(pred_start_to_run[jr]+1 < samples_start.size());
        ulint prev_sample = samples_start[ pred_start_to_run[jr]+1 ];
        return (prev_sample + delta) % this->bwt.size();
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        out.write((char *)&this->terminator_position, sizeof(this->terminator_position));
        written_bytes += sizeof(this->terminator_position);
        written_bytes += my_serialize(this->F, out, child, "F");
        written_bytes += this->bwt.serialize(out);
        written_bytes += this->pred.serialize(out);
        written_bytes += this->pred_to_run.serialize(out);
        written_bytes += this->samples_last.serialize(out);

        written_bytes += thresholds.serialize(out, child, "thresholds");
        // written_bytes += my_serialize(thresholds, out, child, "thresholds");
        // written_bytes += my_serialize(samples_start, out, child, "samples_start");
        written_bytes += pred_start.serialize(out);
        written_bytes += pred_start_to_run.serialize(out, child, "pred_start_to_run");
        written_bytes += samples_start.serialize(out, child, "samples_start");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    std::string get_file_extension() const
    {
        return thresholds.get_file_extension() + ".full.ms";
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in)
    {

        in.read((char *)&this->terminator_position, sizeof(this->terminator_position));
        my_load(this->F, in);
        this->bwt.load(in);
        this->r = this->bwt.number_of_runs();
        this->pred.load(in);
        this->pred_to_run.load(in);
        this->samples_last.load(in);

        thresholds.load(in,&this->bwt);
        // my_load(thresholds, in);
        pred_start.load(in);
        pred_start_to_run.load(in);
        samples_start.load(in);
        // my_load(samples_start,in);
    }

    // // From r-index
    // ulint get_last_run_sample()
    // {
    //     return (samples_last[r - 1] + 1) % bwt.size();
    // }

    // Test methods
    std::vector<size_t> get_SA_Phi()
    {
        std::vector<size_t> sa_phi;

        sa_phi.push_back(this->get_last_run_sample());
        while(sa_phi.size() < this->bwt_size() - 1)
            sa_phi.push_back(this->Phi(sa_phi.back()));

        return sa_phi;
    }
    std::vector<size_t> get_SA_Phi_inv()
    {
        std::vector<size_t> sa_phi_inv;

        sa_phi_inv.push_back((samples_start[1]+1) % this->bwt.size());
        while(sa_phi_inv.size() < this->bwt_size() - 1)
            sa_phi_inv.push_back(Phi_inv(sa_phi_inv.back()));

        return sa_phi_inv;
    }

protected:
    // Computes the matching statistics pointers for the given pattern
    template<typename string_t>
    std::vector<size_t> _query(const string_t &pattern, const size_t m)
    {

        std::vector<size_t> ms_pointers(m);

        // Start with the empty string
        auto pos = this->bwt_size() - 1;
        auto sample = this->get_last_run_sample();

        for (size_t i = 0; i < m; ++i)
        {
            auto c = pattern[m - i - 1];

            if (this->bwt.number_of_letter(c) == 0)
            {
                sample = 0;
            }
            else if (pos < this->bwt.size() && this->bwt[pos] == c)
            {
                sample--;
            }
            else
            {
                // Get threshold
                ri::ulint rnk = this->bwt.rank(pos, c);
                size_t thr = this->bwt.size() + 1;

                ulint next_pos = pos;

                // if (rnk < (this->F[c] - this->F[c-1]) // I can use F to compute it
                if (rnk < this->bwt.number_of_letter(c))
                {
                    // j is the first position of the next run of c's
                    ri::ulint j = this->bwt.select(rnk, c);
                    ri::ulint run_of_j = this->bwt.run_of_position(j);

                    thr = thresholds[run_of_j]; // If it is the first run thr = 0

                    // Here we should use Phi_inv that is not implemented yet
                    // sample = this->Phi(this->samples_last[run_of_j - 1]) - 1;
                    sample = samples_start[run_of_j];

                    next_pos = j;
                }

                if (pos < thr)
                {

                    rnk--;
                    ri::ulint j = this->bwt.select(rnk, c);
                    ri::ulint run_of_j = this->bwt.run_of_position(j);
                    sample = this->samples_last[run_of_j];

                    next_pos = j;
                }

                pos = next_pos;
            }

            ms_pointers[m - i - 1] = sample;

            // Perform one backward step
            pos = LF(pos, c);
        }

        return ms_pointers;
    }
    // // From r-index
    // vector<ulint> build_F(std::ifstream &ifs)
    // {
    //     ifs.clear();
    //     ifs.seekg(0);
    //     F = vector<ulint>(256, 0);
    //     uchar c;
    //     ulint i = 0;
    //     while (ifs >> c)
    //     {
    //         if (c > TERMINATOR)
    //             F[c]++;
    //         else
    //         {
    //             F[TERMINATOR]++;
    //             terminator_position = i;
    //         }
    //         i++;
    //     }
    //     for (ulint i = 255; i > 0; --i)
    //         F[i] = F[i - 1];
    //     F[0] = 0;
    //     for (ulint i = 1; i < 256; ++i)
    //         F[i] += F[i - 1];
    //     return F;
    // }

    // // From r-index
    // vector<pair<ulint, ulint>> &read_run_starts(std::string fname, ulint n, vector<pair<ulint, ulint>> &ssa)
    // {
    //     ssa.clear();
    //     std::ifstream ifs(fname);
    //     uint64_t x = 0;
    //     uint64_t y = 0;
    //     uint64_t i = 0;
    //     while (ifs.read((char *)&x, 5) && ifs.read((char *)&y, 5))
    //     {
    //         ssa.push_back(pair<ulint, ulint>(y ? y - 1 : n - 1, i));
    //         i++;
    //     }
    //     return ssa;
    // }

    // // From r-index
    // vector<ulint> &read_run_ends(std::string fname, ulint n, vector<ulint> &esa)
    // {
    //     esa.clear();
    //     std::ifstream ifs(fname);
    //     uint64_t x = 0;
    //     uint64_t y = 0;
    //     while (ifs.read((char *)&x, 5) && ifs.read((char *)&y, 5))
    //     {
    //         esa.push_back(y ? y - 1 : n - 1);
    //     }
    //     return esa;
    // }
};

// Computes the matching statistics pointers for the given pattern
template <>
template <typename string_t>
std::vector<size_t> ms_pointers<ri::sparse_sd_vector, ms_rle_string_sd, thr_bv<ms_rle_string_sd>>::_query(const string_t &pattern, const size_t m)
{

    std::vector<size_t> ms_pointers(m);

    // Start with the empty string
    auto pos = this->bwt_size() - 1;
    auto sample = this->get_last_run_sample();

    for (size_t i = 0; i < m; ++i)
    {
        auto c = pattern[m - i - 1];
        const auto n_c = this->bwt.number_of_letter(c);
        if (n_c == 0)
        {
            sample = 0;
            // Perform one backward step
            pos = LF(pos, c);
        }
        else if (pos < this->bwt.size() && this->bwt[pos] == c)
        {
            sample--;
            // Perform one backward step
            pos = LF(pos, c);
        }
        else
        {
            // Get threshold
            ri::ulint run_of_pos = this->bwt.run_of_position(pos);
            auto rnk_c = this->bwt.run_and_head_rank(run_of_pos, c);
            size_t thr_c = thresholds.rank(pos + 1,c); // +1 because the rank count the thresiold in pos

            if(rnk_c.first > thr_c)
            {
                // Jump up
                size_t run_of_j = this->bwt.run_head_select(rnk_c.first,c);
                sample = samples_last[run_of_j];
                // Perform one backward step
                pos = this->F[c] + rnk_c.second - 1;
            }
            else
            {
                // Jump down
                size_t run_of_j = this->bwt.run_head_select(rnk_c.first + 1, c);
                sample = samples_start[run_of_j];
                // Perform one backward step
                pos = this->F[c] + rnk_c.second;
            }
        }
        // Store the sample
        ms_pointers[m - i - 1] = sample;
    }

    return ms_pointers;
}

#endif /* end of include guard: _MS_POINTERS_HH */
