/* chain - Chain MEMs
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
   \file chain.hpp
   \brief chain.hpp Chain MEMs
   \author Massimiliano Rossi
   \date 19/10/2021
*/

#ifndef _CHAIN_HH
#define _CHAIN_HH

#include <mems.hpp>

typedef struct{
    ll score = 0;
    size_t mate = 2;
    bool paired = false;
    bool reversed = false;
    std::vector<size_t> anchors;

    // Lazy reverse of the chain
    void reverse()
    {
        if(not reversed)
        {
            std::reverse(anchors.begin(), anchors.end());
            reversed = true;
        }
    }

    void reset()
    {
        if (reversed)
        {
            std::reverse(anchors.begin(), anchors.end());
            reversed = false;
        }
    }

} chain_t;

// typedef struct{
//     size_t pos;
//     size_t len;
//     size_t mem_i;
//     size_t occ_i;

//     bool operator<=(const anchor_t &other)
//     {
//       return ((this.pos + this.len -1) < (other.pos + other.len -1));
//     }

//     friend bool operator>(const anchor_t& lhs, const anchor_t& rhs)
//     {
//       return ((this.pos + this.len -1) >= (other.pos + other.len -1));
//     }

// } anchor_t;

typedef struct{
    ll G = LLONG_MAX;
    ll max_dist_x = 500;//LLONG_MAX;
    ll max_dist_y = 100;//LLONG_MAX;
    ll max_iter = 50;
    ll max_pred = 50;
    ll min_chain_score = 40;
    ll min_chain_length = 1;
} chain_config_t;


void populate_anchors(
    std::vector< std::pair< size_t, size_t > >& anchors,
    const std::vector<mem_t>& mems,
    size_t& tot_mem_length
    )
{
    for(size_t i = 0; i < mems.size(); ++i) {
        for(size_t j = 0; j < mems[i].occs.size(); ++j){
            anchors.push_back(make_pair(i,j));
        }
        tot_mem_length +=  mems[i].len * mems[i].occs.size();
    }
}

inline void update_score_and_pred(
    std::vector<ll>& f,
    std::vector<ll>& p,
    std::vector<ll>& msc,
    const ll max_f,
    const ll max_j,
    const size_t i
    )
{
    f[i] = max_f;
    p[i] = max_j;
    if( max_j >= 0 and msc[max_j] > max_f)
        msc[i] = msc[max_j];
    else
        msc[i] = max_f;
}

void find_chain_ends(
    const std::vector< std::pair< size_t, size_t > >& anchors,
    const std::vector<ll>& p,
    std::vector<ll>& t
    )
{
    for(size_t i = 0; i < anchors.size(); ++i){
        if(p[i] >= 0){
            t[p[i]] = 1;
        } 
    }
}

size_t find_n_chains(
    size_t n_chains,
    const std::vector< std::pair< size_t, size_t > >& anchors,
    const std::vector<ll>& t,
    const std::vector<ll>& msc,
    const ll min_chain_score 
    )
{
    for(size_t i = 0; i < anchors.size(); ++i){
        if(t[i] == 0 and msc[i] > min_chain_score){
            n_chains ++;
        } 
    }
    return n_chains;
}

void find_chain_starts(
    std::vector<std::pair<ll, size_t>>& chain_starts,
    size_t& k,
    const std::vector<ll>& t,
    const std::vector<ll>& f,
    const std::vector<ll>& p,
    const std::vector<ll>& msc,
    const std::vector< std::pair< size_t, size_t > >& anchors,
    const ll min_chain_score
    )
{
    for(size_t i = 0; i < anchors.size(); ++i)
    {
        if(t[i] == 0 and  msc[i] > min_chain_score)
        {
            size_t j = i;
            while(j >= 0 and f[j] < msc[j]) j = p[j]; // Find the peak that maximizes f
            if (j < 0) i = 0;
            chain_starts[k++] = make_pair(f[j],j);
        }
    }
}

void backtrack(
    const std::vector<std::pair<ll, size_t>>& chain_starts,
    const size_t n_chains,
    std::vector<ll>& t,
    const std::vector<ll>& f,
    const std::vector<ll>& p,
    std::vector<chain_t>& chains,
    const std::vector<std::pair<size_t,size_t>>& anchors,
    const std::vector<mem_t>& mems,
    const ll min_chain_score,
    const ll min_chain_length
    )
{
    for (size_t i = 0; i < n_chains; ++i) {
        ll j = chain_starts[i].second;
        chain_t chain;
        chain.mate = mems[anchors[j].first].mate;
        chain.score = chain_starts[i].first;
        do {
            chain.paired = chain.paired or (chain.mate != mems[anchors[j].first].mate);
            chain.anchors.push_back(j); // stores the reverse of the chain
            t[j] = 1;
            j = p[j];
        } while (j >= 0 && t[j] == 0);
        if (j < 0) { // l - prev_l is the length of the chain
            if (chain.anchors.size() >= min_chain_length){
                chains.push_back(std::move(chain));
            }     
        } else if (chain_starts[i].first - f[j] >= min_chain_score) { // Two chains share a common prefix
            if (chain.anchors.size() >= min_chain_length){
                chains.push_back(std::move(chain));
            }
        }
    }
}

void clear_and_shrink_vecs(
    std::vector<ll>& p,
    std::vector<ll>& t,
    std::vector<ll>& f,
    std::vector<ll>& msc
    )
{
    p.resize(0);
    p.shrink_to_fit();
    t.resize(0);
    t.shrink_to_fit();
    f.resize(0);
    f.shrink_to_fit();
    msc.resize(0);
    msc.shrink_to_fit();
}


// Given a set of mems, find the chains
bool find_chains(
    const   std::vector<mem_t>& mems,
            std::vector< std::pair< size_t, size_t > >& anchors,
            std::vector<chain_t>& chains,
            const chain_config_t config = chain_config_t()        
    )
{
    MTIME_INIT(7);   
    MTIME_START(3); // Timing helper
    /************* minimap2 dynamic programming for mem chaining ***************/
    /* https://github.com/lh3/minimap2/blob/master/chain.c */

    // Sort anchors
    // Lamda helper to sort the anchors
    auto cmp = [&] (const std::pair<size_t, size_t>& i, const std::pair<size_t, size_t>& j) -> bool {
    return (mems[i.first].occs[i.second] + mems[i.first].len - 1) < (mems[j.first].occs[j.second] + mems[j.first].len - 1);
    };

    // TODO: improve this initialization
    size_t tot_mem_length = 0;

    populate_anchors(anchors, mems, tot_mem_length);
    
    float avg_mem_length = (float)tot_mem_length / anchors.size();

    std::sort(anchors.begin(),anchors.end(),cmp);
    MTIME_END(3); //Timing helper
    MTIME_START(4); //Timing helper
    // MTIME_TSAFE_MERGE;

    // Dynamic programming
    const ll G = config.G;
    const ll max_dist_x = config.max_dist_x;
    const ll max_dist_y = config.max_dist_y;
    const ll max_iter = config.max_iter;
    const ll max_pred = config.max_pred;
    const ll min_chain_score = config.min_chain_score;
    const ll min_chain_length = config.min_chain_length;

    // // TODO: Parameters to be defined
    // const ll G = LLONG_MAX;
    // const ll max_dist_x = 500;//LLONG_MAX;
    // // const ll max_dist_x = 100;//LLONG_MAX;
    // const ll max_dist_y = 100;//LLONG_MAX;
    // const ll max_iter = 50;
    // const ll max_pred = 50;
    // const ll min_chain_score = 1;
    // const ll min_chain_length = 1;


    std::vector<ll> f(anchors.size(),0); // Score ending in position i
    std::vector<ll> p(anchors.size(),0); // Position of the next anchor giving the max when chained with the one in position i
    std::vector<ll> msc(anchors.size(),0); // Max score up to position i
    std::vector<ll> t(anchors.size(),0); // Stores i in position p[j], if i is chained with j. See heuristics in minimap2

    ll lb = 0;
    // For all the anchors
    for(size_t i = 0 ; i < anchors.size(); ++i)
    {
    // Get anchor i
    const auto a_i = anchors[i];
    const mem_t& mem_i = mems[a_i.first];
    // const ll k_i = mem_i.occs[a_i.second];
    const ll x_i = mem_i.occs[a_i.second] + mem_i.len - 1;
    // const ll z_i = mem_i.idx;
    const ll y_i = mem_i.rpos;
    // const ll y_i = mem_i.idx + mem_i.len - 1;
    const ll w_i = mem_i.len;
    const size_t mate_i = mem_i.mate;

    ll max_f = w_i;
    ll max_j = -1;
    size_t n_pred = 0;
    // For all previous anchors
    // Heuristics from minimap2 -> do not try more than 50 anchors
    if(i - lb > max_iter) lb = i - max_iter;
    for(ll j = i-1; j >= lb; --j)
    {
        const auto a_j = anchors[j];
        const mem_t& mem_j = mems[a_j.first];
        const ll x_j = mem_j.occs[a_j.second] + mem_j.len - 1;
        const ll y_j = mem_j.rpos;
        // const ll y_j = mem_j.idx + mem_j.len - 1;
        const size_t mate_j = mem_j.mate;

        // Check if anchors are compatible
        // if mate_i == mate_j they are from the same mate and same orientation
        // if mate_i ^ mate_j == 3 they are from different mate and withdifferent orientation
        // TODO: we can change to same orientation replacing 3 with 1
        if((mate_i != mate_j) and ((mate_i ^ mate_j)!= 3)) continue;

        // If the current anchor is too far, exit.
        if(x_i > x_j + max_dist_x)
        {
            // j = lb - 1;
            lb = j;
            continue;
        }
        // // skip if incompatible
        // if(k_i < x_j or z_i < y_j) continue;

        const ll x_d = x_i - x_j;
        const ll y_d = y_i - y_j;
        const int32_t l = (y_d > x_d ? (y_d - x_d) : (x_d - y_d));
        const uint32_t ilog_l = (l > 0? ilog2_32(l): 0);

        if((mate_i == mate_j and (y_j >= y_i or y_d > max_dist_y)) or max(y_d, x_d) > G)
        continue;

        const ll alpha = min(min(y_d,x_d),w_i);
        // const ll beta = (l > 0? (ll)(.01 * l * avg_mem_length) + ilog_l >> 1 : 0); 
        ll beta = 0;
        if(mate_i != mate_j){
            if (x_d == 0) ++beta; // possibly due to overlapping paired ends; give a minor bonus
            else {
                const int c_lin = (int)(l * .01 * avg_mem_length);
                beta = c_lin < ilog_l ? c_lin : ilog_l;
            }
        }else{
            beta = (l > 0? (ll)(.01 * l * avg_mem_length) + ilog_l >> 1 : 0);
        }



        // No gap scale as in minimap2
        ll score = f[j] + (alpha - beta);

        if( score > max_f )
        {
        max_f = score;
        max_j = j;
        if( n_pred > 0) --n_pred;
        }
        else // minimap2: If i is chained wth j, than chaining i with a predecessor of j does not improve the score
        if (t[j] == i and (++n_pred > max_pred))            
            break;
        
        if(p[j] > 0) t[p[j]] = i;
    }

    update_score_and_pred(f, p, msc, max_f, max_j, i);
    }

    MTIME_END(4);   // Timing helper
    MTIME_START(5); // Timing helper

    // Find the end positions of chains
    memset(t.data(),0,sizeof(size_t) * t.size());
    find_chain_ends(anchors, p, t);
    
    size_t n_chains = 0;
    n_chains = find_n_chains(n_chains, anchors, t, msc, min_chain_score);
    
    // TODO: check if we want to report also non aligned reads
    if(n_chains == 0)
    return false;

    // TODO: replace the vector of pairs with a lambda for the sort
    std::vector<std::pair<ll,size_t>> chain_starts(n_chains); // Stores uniqe chains and their index of uniqe chains in anchors
    size_t k = 0;
    find_chain_starts(chain_starts, k, t, f, p, msc, anchors, min_chain_score);

    chain_starts.resize(k), chain_starts.shrink_to_fit();
    n_chains = chain_starts.size();
    std::sort(chain_starts.begin(), chain_starts.end(), std::greater<std::pair<ll,size_t>>());
    
    // std::vector<std::pair<size_t, std::vector<size_t>>> chains;
    MTIME_END(5);   // Timing helper
    MTIME_START(6); // Timing helper

    // Backtrack
    memset(t.data(),0,sizeof(size_t) * t.size());
    backtrack(chain_starts, n_chains, t, f, p, chains, anchors, mems, min_chain_score, min_chain_length);

    // Lamda helper to sort the anchors
    auto chain_t_cmp = [](const chain_t& i, const chain_t& j) -> bool
    {
    return i.score > j.score;
    };

    // Sort the chains by max scores.
    std::sort(chains.begin(), chains.end(), chain_t_cmp);

    // // Backtrack
    // memset(t.data(),0,sizeof(size_t) * t.size());
    // for(size_t i = 0; i < n_chains; ++i)
    // {
    //   ll j = chain_starts[i].second;
    //   std::vector<size_t> chain;
    //   bool paired = true;
    //   size_t mate = mems[anchors[j].first].mate; 
    //   do {
    //     size_t tmp_mate = mems[anchors[j].first].mate;
    //     paired = paired and (mate == tmp_mate);
    //     chain.push_back(j); // stores th reverse of the chain
    //     t[j] = 1;
    //     j = p[j];
    //   } while (j >= 0 && t[j] == 0);
    //   if (j < 0) { // l - prev_l is the length of the chain
    //     if (chain.size() >= min_chain_length)
    //       chains.push_back(std::make_pair(chain_starts[i].first, chain));     
    //   } else if (chain_starts[i].first - f[j] >= min_chain_score) { // Two chains share a common prefix
    //     if (chain.size() >= min_chain_length)
    //       chains.push_back(std::make_pair((chain_starts[i].first - f[j]), chain));
    //   }
    // }

    // // Sort the chains by max scores.
    // std::sort(chains.begin(), chains.end(), std::greater<std::pair<ll,std::vector<size_t>>>());

    // Clear space
    clear_and_shrink_vecs(p, t, f, msc);

    MTIME_END(6); // Timing helper
    MTIME_TSAFE_MERGE;

    return true;
}


// Given a set of mems, find the chains
bool find_chains_secondary(
    const   std::vector<mem_t>& mems,
            std::vector< std::pair< size_t, size_t > >& anchors,
            std::vector<chain_t>& chains,
            const chain_config_t config = chain_config_t()        
    )
{
    MTIME_INIT(7);   
    MTIME_START(3); // Timing helper
    /************* minimap2 dynamic programming for mem chaining ***************/
    /* https://github.com/lh3/minimap2/blob/master/chain.c */

    // Sort anchors
    // Lamda helper to sort the anchors
    auto cmp = [&] (const std::pair<size_t, size_t>& i, const std::pair<size_t, size_t>& j) -> bool {
    return (mems[i.first].occs[i.second] + mems[i.first].len - 1) < (mems[j.first].occs[j.second] + mems[j.first].len - 1);
    };

    // TODO: improve this initialization
    size_t tot_mem_length = 0;

    populate_anchors(anchors, mems, tot_mem_length);
    
    float avg_mem_length = (float)tot_mem_length / anchors.size();

    std::sort(anchors.begin(),anchors.end(),cmp);
    MTIME_END(3); //Timing helper
    MTIME_START(4); //Timing helper
    // MTIME_TSAFE_MERGE;

    // Dynamic programming
    const ll G = config.G;
    const ll max_dist_x = config.max_dist_x;
    const ll max_dist_y = config.max_dist_y;
    const ll max_iter = config.max_iter;
    const ll max_pred = config.max_pred;
    const ll min_chain_score = config.min_chain_score;
    const ll min_chain_length = config.min_chain_length;

    // // TODO: Parameters to be defined
    // const ll G = LLONG_MAX;
    // const ll max_dist_x = 500;//LLONG_MAX;
    // // const ll max_dist_x = 100;//LLONG_MAX;
    // const ll max_dist_y = 100;//LLONG_MAX;
    // const ll max_iter = 50;
    // const ll max_pred = 50;
    // const ll min_chain_score = 1;
    // const ll min_chain_length = 1;


    std::vector<ll> f(anchors.size(),0); // Score ending in position i
    std::vector<ll> f_sec(anchors.size(),0); //Secondary score ending in position i
    std::vector<ll> p(anchors.size(),0); // Position of the next anchor giving the max when chained with the one in position i
    std::vector<ll> p_sec(anchors.size(),0); // Position of the next anchor giving the secondary max when chained with the one in position i
    std::vector<ll> msc(anchors.size(),0); // Max score up to position i
    std::vector<ll> msc_sec(anchors.size(),0); // Secondary max score up to position i
    std::vector<ll> t(anchors.size(),0); // Stores i in position p[j], if i is chained with j. See heuristics in minimap2
    std::vector<ll> t_sec(anchors.size(),0); // Stores i in position p[j], if i is chained with j secondary. See heuristics in minimap2

    ll lb = 0;
    // For all the anchors
    for(size_t i = 0 ; i < anchors.size(); ++i)
    {
        // Get anchor i
        const auto a_i = anchors[i];
        const mem_t& mem_i = mems[a_i.first];
        // const ll k_i = mem_i.occs[a_i.second];
        const ll x_i = mem_i.occs[a_i.second] + mem_i.len - 1;
        // const ll z_i = mem_i.idx;
        const ll y_i = mem_i.rpos;
        // const ll y_i = mem_i.idx + mem_i.len - 1;
        const ll w_i = mem_i.len;
        const size_t mate_i = mem_i.mate;

        ll max_f = w_i;
        ll max_sec_f = w_i;

        ll max_j = -1;
        ll max_sec_j = -1;

        size_t n_pred = 0;
        // For all previous anchors
        // Heuristics from minimap2 -> do not try more than 50 anchors
        if(i - lb > max_iter) lb = i - max_iter;
        for(ll j = i-1; j >= lb; --j)
        {
            const auto a_j = anchors[j];
            const mem_t& mem_j = mems[a_j.first];
            const ll x_j = mem_j.occs[a_j.second] + mem_j.len - 1;
            const ll y_j = mem_j.rpos;
            // const ll y_j = mem_j.idx + mem_j.len - 1;
            const size_t mate_j = mem_j.mate;

            // Check if anchors are compatible
            // if mate_i == mate_j they are from the same mate and same orientation
            // if mate_i ^ mate_j == 3 they are from different mate and withdifferent orientation
            // TODO: we can change to same orientation replacing 3 with 1
            if((mate_i != mate_j) and ((mate_i ^ mate_j)!= 3)) continue;

            // If the current anchor is too far, exit.
            if(x_i > x_j + max_dist_x)
            {
                // j = lb - 1;
                lb = j;
                continue;
            }
            // // skip if incompatible
            // if(k_i < x_j or z_i < y_j) continue;

            const ll x_d = x_i - x_j;
            const ll y_d = y_i - y_j;
            const int32_t l = (y_d > x_d ? (y_d - x_d) : (x_d - y_d));
            const uint32_t ilog_l = (l > 0? ilog2_32(l): 0);

            if((mate_i == mate_j and (y_j >= y_i or y_d > max_dist_y)) or max(y_d, x_d) > G)
            continue;

            const ll alpha = min(min(y_d,x_d),w_i);
            // const ll beta = (l > 0? (ll)(.01 * l * avg_mem_length) + ilog_l >> 1 : 0); 
            ll beta = 0;
            if(mate_i != mate_j){
                if (x_d == 0) ++beta; // possibly due to overlapping paired ends; give a minor bonus
                else {
                    const int c_lin = (int)(l * .01 * avg_mem_length);
                    beta = c_lin < ilog_l ? c_lin : ilog_l;
                }
            }else{
                beta = (l > 0? (ll)(.01 * l * avg_mem_length) + ilog_l >> 1 : 0);
            }



            // No gap scale as in minimap2
            ll score = f[j] + (alpha - beta);
            ll score_sec = f_sec[j] + (alpha - beta);

            // Primary Chain of anchor i
            if( score > max_f )
            {
                max_f = score;
                max_j = j;
                if( n_pred > 0) --n_pred;
            }
            // Secondary Chain of anchor i
            else if (score_sec > max_sec_f)
            {
                // Check if max_j of primary chain is defined 
                if (max_j >= 0){
                    ll tmp = max_j;
                    bool uniq_chain = true;
                    const size_t& mem_j_pos = mems[anchors[j].first].occs[anchors[j].second];

                    // Check whether potential j anchor for secondary chain is used in primary chain
                    // Secondary chain should not use anchors in primary chain, otherwise secondary chain can potentially be a copy of primary chain
                    while(tmp >= 0)
                    {
                        const size_t& mem_tmp_pos = mems[anchors[tmp].first].occs[anchors[tmp].second];

                        if (mem_j_pos == mem_tmp_pos)
                        {
                            uniq_chain = false;
                            break;
                        }
                        
                        tmp = p[tmp];
                    }
                    if (uniq_chain){ // If anchor is not used in primary chain then set secondary score.
                        max_sec_f = score_sec;
                        max_sec_j = j;
                    }
                }
            }
            else // minimap2: If i is chained wth j, than chaining i with a predecessor of j does not improve the score
            {
                if (t[j] == i and (++n_pred > max_pred)) break;
            }
            
            if(p[j] > 0) t[p[j]] = i;
            if(p_sec[j] > 0) t_sec[p_sec[j]] = i;
        }

        update_score_and_pred(f, p, msc, max_f, max_j, i);
        update_score_and_pred(f_sec, p_sec, msc_sec, max_sec_f, max_sec_j, i);
    }

    MTIME_END(4);   // Timing helper
    MTIME_START(5); // Timing helper

    // Find the end positions of chains
    memset(t.data(),0,sizeof(size_t) * t.size());
    memset(t_sec.data(),0,sizeof(size_t) * t_sec.size());
    find_chain_ends(anchors, p, t);
    find_chain_ends(anchors, p_sec, t_sec);
    
    
    size_t n_chains = 0;
    size_t n_chains_sec = 0;
    n_chains = find_n_chains(n_chains, anchors, t, msc, min_chain_score);
    n_chains_sec = find_n_chains(n_chains_sec, anchors, t_sec, msc_sec, min_chain_score);
    
    // TODO: check if we want to report also non aligned reads
    if(n_chains == 0)
    return false;

    // TODO: replace the vector of pairs with a lambda for the sort
    std::vector<std::pair<ll,size_t>> chain_starts(n_chains); // Stores uniqe chains and their index of uniqe chains in anchors
    std::vector<std::pair<ll,size_t>> chain_starts_sec(n_chains_sec); // Stores uniqe secondary chains and their index of uniqe secondary chains in anchors

    // Find chain starts for primary chains
    size_t k = 0;
    find_chain_starts(chain_starts, k, t, f, p, msc, anchors, min_chain_score);
    // Find chain starts for secondary chains
    size_t k_sec = 0;
    find_chain_starts(chain_starts_sec, k_sec, t_sec, f_sec, p_sec, msc_sec, anchors, min_chain_score);

    chain_starts.resize(k), chain_starts.shrink_to_fit();
    chain_starts_sec.resize(k_sec), chain_starts_sec.shrink_to_fit();

    n_chains = chain_starts.size();
    n_chains_sec = chain_starts_sec.size();

    // Lambda for sorting chain starts 
    auto chain_start_cmp = [](const std::pair<ll, size_t>& lhs, const std::pair<ll, size_t>& rhs){
        return lhs.first > rhs.first;
    };

    std::sort(chain_starts.begin(), chain_starts.end(), chain_start_cmp);
    std::sort(chain_starts_sec.begin(), chain_starts_sec.end(), chain_start_cmp);
    
    // std::vector<std::pair<size_t, std::vector<size_t>>> chains;
    MTIME_END(5);   // Timing helper
    MTIME_START(6); // Timing helper

    // Backtrack for primary chains
    memset(t.data(),0,sizeof(size_t) * t.size());
    backtrack(chain_starts, n_chains, t, f, p, chains, anchors, mems, min_chain_score, min_chain_length);

    // Backtrack for secondary chains
    memset(t_sec.data(),0,sizeof(size_t) * t_sec.size());
    backtrack(chain_starts_sec, n_chains_sec, t_sec, f_sec, p_sec, chains, anchors, mems, min_chain_score, min_chain_length);

    // Lamda helper to sort the anchors
    auto chain_t_cmp = [](const chain_t& i, const chain_t& j) -> bool
    {
        return i.score > j.score;
    };

    // Sort the chains by max scores.
    std::sort(chains.begin(), chains.end(), chain_t_cmp);

    // // Backtrack
    // memset(t.data(),0,sizeof(size_t) * t.size());
    // for(size_t i = 0; i < n_chains; ++i)
    // {
    //   ll j = chain_starts[i].second;
    //   std::vector<size_t> chain;
    //   bool paired = true;
    //   size_t mate = mems[anchors[j].first].mate; 
    //   do {
    //     size_t tmp_mate = mems[anchors[j].first].mate;
    //     paired = paired and (mate == tmp_mate);
    //     chain.push_back(j); // stores th reverse of the chain
    //     t[j] = 1;
    //     j = p[j];
    //   } while (j >= 0 && t[j] == 0);
    //   if (j < 0) { // l - prev_l is the length of the chain
    //     if (chain.size() >= min_chain_length)
    //       chains.push_back(std::make_pair(chain_starts[i].first, chain));     
    //   } else if (chain_starts[i].first - f[j] >= min_chain_score) { // Two chains share a common prefix
    //     if (chain.size() >= min_chain_length)
    //       chains.push_back(std::make_pair((chain_starts[i].first - f[j]), chain));
    //   }
    // }

    // // Sort the chains by max scores.
    // std::sort(chains.begin(), chains.end(), std::greater<std::pair<ll,std::vector<size_t>>>());

    // Clear space
    clear_and_shrink_vecs(p, t, f, msc); 
    clear_and_shrink_vecs(p_sec, t_sec, f_sec, msc_sec);

    MTIME_END(6); // Timing helper
    MTIME_TSAFE_MERGE;

    return true;
}



#endif /* end of include guard: _CHAIN_HH */