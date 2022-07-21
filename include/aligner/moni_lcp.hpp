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

#ifndef _MONI_LCP_HH
#define _MONI_LCP_HH

#include <common.hpp>

#include <malloc_count.h>

#include <moni.hpp>

template <class sparse_bv_type = ri::sparse_sd_vector,
          class rle_string_t = ms_rle_simple_string_sd,
          class thresholds_t = thr_bv<rle_string_t> >
class moni_lcp : public ms_pointers<sparse_bv_type, rle_string_t, thresholds_t>
{
public:
    int_vector<> slcp;

    typedef size_t size_type;

    moni_lcp() {}

    moni_lcp(std::string filename, bool rle = false) : ms_pointers<sparse_bv_type, rle_string_t, thresholds_t>()
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
        this->build_phi(filename + ".ssa", this->r, n, this->pred, this->pred_to_run);

        verbose("Building inverse of phi");
        this->build_phi(filename + ".esa", this->r, n, this->pred_start, this->pred_start_to_run);

        this->read_samples(filename + ".ssa", this->r, n, this->samples_start);
        this->read_samples(filename + ".esa", this->r, n, this->samples_last);

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("R-index construction complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

        verbose("Reading thresholds from file");

        t_insert_start = std::chrono::high_resolution_clock::now();

        this->thresholds = thresholds_t(filename, &this->bwt);

        std::string tmp_filename = filename + std::string(".slcp");

        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
            error("open() file " + tmp_filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + tmp_filename + " failed");

        if (filestat.st_size % SSABYTES != 0)
            error("invilid file " + tmp_filename);

        size_t length = filestat.st_size / SSABYTES;
        slcp.resize(length);

        for (size_t i = 0; i < length; ++i){
            size_t tmp = 0;
            if ((fread(&tmp, SSABYTES, 1, fd)) != 1)
                error("fread() file " + tmp_filename + " failed");
            slcp[i] = tmp;
        }

        fclose(fd);

        // compress slcp array
        sdsl::util::bit_compress(slcp);

        t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    }

   
    void print_stats()
    {
        sdsl::nullstream ns;

        verbose("Memory consumption (bytes).");
        verbose("   terminator_position: ", sizeof(this->terminator_position));
        verbose("                     F: ", my_serialize(this->F, ns));
        verbose("                   bwt: ", this->bwt.serialize(ns));
        verbose("            thresholds: ", this->thresholds.serialize(ns));
        verbose("                  pred: ", this->pred.serialize(ns));
        verbose("           pred_to_run: ", this->pred_to_run.serialize(ns));
        verbose("          samples_last: ", this->samples_last.serialize(ns));
        verbose("            pred_start: ", this->pred_start.serialize(ns));
        verbose("     pred_start_to_run: ", this->pred_start_to_run.serialize(ns));
        verbose("         samples_start: ", this->samples_start.serialize(ns));
        verbose("                  slcp: ", slcp.serialize(ns));
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

        written_bytes += this->thresholds.serialize(out, child, "thresholds");
        // written_bytes += my_serialize(thresholds, out, child, "thresholds");
        // written_bytes += my_serialize(samples_start, out, child, "samples_start");
        written_bytes += this->pred_start.serialize(out);
        written_bytes += this->pred_start_to_run.serialize(out, child, "pred_start_to_run");
        written_bytes += this->samples_start.serialize(out, child, "samples_start");
        written_bytes += slcp.serialize(out, child, "slcp");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    std::string get_file_extension() const
    {
        return this->thresholds.get_file_extension() + ".full.lcp.ms";
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

        this->thresholds.load(in,&this->bwt);
        // my_load(thresholds, in);
        this->pred_start.load(in);
        this->pred_start_to_run.load(in);
        this->samples_start.load(in);
        slcp.load(in);
    }

    /*
     * Phi function. Phi_inv(SA[n-1]) is undefined
     */
    std::pair<ulint,ulint> Phi_inv_lcp(ulint i)
    {
        assert(i != this->bwt.size() - 1);
        // jr is the rank of the predecessor of i (circular)
        ulint jr = this->pred_start.predecessor_rank_circular(i);
        assert(jr <= this->r - 1);
        // the actual predecessor
        ulint j = this->pred_start.select(jr);
        assert(jr < this->r - 1 or j == this->bwt.size() - 1);
        // distance from predecessor
        ulint delta = j < i ? i - j : i + 1;
        // cannot fall on first run: this can happen only if I call Phi(SA[n-1])
        assert(this->pred_start_to_run[jr] < this->samples_start.size() - 1);
        // sample at the end of previous run
        assert(this->pred_start_to_run[jr] + 1 < this->samples_start.size());
        ulint prev_sample = this->samples_start[this->pred_start_to_run[jr] + 1];
        ulint lcp = slcp[this->pred_start_to_run[jr] + 1];
        return {(prev_sample + delta) % this->bwt.size(),lcp - delta + 1};
    }

    /*
     * Phi function. Phi(SA[0]) is undefined
     */
    std::pair<ulint, ulint> Phi_lcp(ulint i)
    {
        assert(i != this->bwt.size() - 1);
        // jr is the rank of the predecessor of i (circular)
        ulint jr = this->pred.predecessor_rank_circular(i);
        assert(jr <= this->r - 1);
        // the actual predecessor
        ulint j = this->pred.select(jr);
        assert(jr < this->r - 1 or j == this->bwt.size() - 1);
        // distance from predecessor
        ulint delta = j < i ? i - j : i + 1;
        // cannot fall on first run: this can happen only if I call Phi(SA[0])
        assert(this->pred_to_run[jr] > 0);
        // sample at the end of previous run
        assert(this->pred_to_run[jr] - 1 < this->samples_last.size());
        ulint idx = this->pred_to_run[jr] - 1;
        ulint prev_sample = this->samples_last[idx];
        ulint lcp = slcp[idx+1];
        return {(prev_sample + delta) % this->bwt.size(), lcp - delta + 1};
    }

protected:

};

#endif /* end of include guard: _MONI_LCP_HH */
