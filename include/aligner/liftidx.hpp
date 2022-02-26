/* liftidx - lifting functions between references
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
   \file liftidx.cpp
   \brief liftidx.cpp lifting functions between references.
   \author Massimiliano Rossi
   \date 22/11/2021
*/
#ifndef _LIFTIDX_HH
#define _LIFTIDX_HH

// LevioSAM has to be included before common because common define verbose
// which is a variable within LevioSAM
#include <leviosam.hpp> 
#include <common.hpp>
#include <seqidx.hpp>


class liftidx : public  seqidx
{
public:
    liftidx(): seqidx()
    {

    }

    // TODO: Include also the list of chains
    /**
     * @brief Construct a new liftidx object from onset, list of sequence names and total length
     *        I am assuming that the lifts are in the same order of the names.
     *        I am also assuming a lift for the reference sequence that return the same position.
     * 
     * @param onset the popsitions
     * @param names_ the list of names
     * @param lifts_ the list of lifts
     * @param l the total length
     */
    liftidx( const std::vector<size_t>& onset, 
            const std::vector<std::string>& names_, 
            std::vector<std::pair<lift::Lift,size_t>>& lifts_, 
            const size_t l
        ):
        seqidx(onset,names_,l),
        lifts(lifts_)
    {
        // NtD
    }

    /**
     * @brief Construct a new liftidx object from onset, list of sequence names and total length
     *        I am assuming that the lifts are in the same order of the names.
     *        I am also assuming a lift for the reference sequence that return the same position.
     * 
     * @param idx seqidx data
     */
    liftidx( const seqidx& idx
        ):
        seqidx(idx)
    {
        const size_t n_seqs = names.size();
        size_t clen = 0;
        lifts.reserve(n_seqs);
        for (size_t i = 0; i < n_seqs; ++i)
        {
            const size_t len = length(i);
            lifts.push_back(std::make_pair(get_null_lift(len),clen));
            clen += len;
        }
    }

    /**
     * @brief Lift the position to the reference
     * 
     * @param pos the position to be lifted
     * @return size_t the position in the reference
     */
    inline size_t lift(const size_t pos)
    {
        size_t rank = rank1(pos + 1);
        size_t start = pos - select1(rank);
        auto& lift_ = lifts[rank - 1];
        return lift_.second + lift_.first.lift_pos(start);
    }

    /**
     * @brief check if two positions corresponds to the same region in the reference.
     * 
     * @param pos1 the first position.
     * @param pos2 the second position.
     * @return true if pos1 and pos2 corresponds to the same reference.
     * @return false if pos1 and pos2 does not correspond to the same reference. 
     */
    inline bool same_region(const size_t pos1, const size_t pos2)
    {
        size_t rank_1 = this->rank1(pos1 + 1);
        size_t rank_2 = this->rank1(pos2 + 1);
        auto& lift1 = lifts[rank_1 - 1];
        auto& lift2 = lifts[rank_2 - 1];
        if ( lift1.second != lift2.second ) return false;

        size_t start1 = pos1 - this->select1(rank_1);
        size_t start2 = pos2 - this->select1(rank_2);
        return lift1.first.lift_pos(start1) == lift2.first.lift_pos(start2);
    }

    size_t serialize(std::ostream &out)
    {

        size_t w_bytes = seqidx::serialize(out);
        w_bytes += sdsl::serialize(lifts.size(), out);
        for(size_t i = 0; i < lifts.size(); ++i){
            w_bytes += sdsl::serialize(lifts[i].second, out);
            w_bytes += lifts[i].first.serialize(out);
        }
        
        return w_bytes;
    }

    void load(std::istream &in)
    {

        seqidx::load(in);
        size_t lifts_size;
        sdsl::load(lifts_size, in);
        lifts.resize(lifts_size);
        for (size_t i = 0; i < lifts.size(); ++i)
        {
            sdsl::load(lifts[i].second, in);
            lifts[i].first.load(in);
        }
    }

    std::string get_file_extension() const
    {
        return  ".ldx";
    }

    inline lift::Lift get_null_lift(const size_t len)
    {
        sdsl::bit_vector ibv(len);
        sdsl::bit_vector dbv(len);
        sdsl::bit_vector sbv(len);

        return lift::Lift(ibv, dbv, sbv);
    }

    inline void lift_cigar(bam1_t* b, const size_t pos)
    {
        size_t rank = rank1(pos + 1);
        auto &lift_ = lifts[rank - 1];
        lift_.first.lift_cigar(b);
    }

protected:
// Pairs of lift and start position of the reference
    std::vector<std::pair<lift::Lift,size_t>> lifts; 

};

#endif /* end of include guard: _LIFTIDX_HH */
