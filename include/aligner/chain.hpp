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

typedef long long int ll;

#define MATE_1 0 // The MEM cames from mate_1
#define MATE_2 1 // The MEM cames from mate_2
#define MATE_F 0 // The MEM cames from the forward strand
#define MATE_RC 2 // The MEM cames from the reverse-complement strand

typedef struct mem_t{
    size_t pos = 0;  // Position in the reference
    size_t len = 0;  // Length
    size_t idx = 0;  // Position in the pattern
    size_t mate = 0; // Left mate (0) or Right mate (1)
    size_t rpos = 0; // Position in the read for chaining
                        // With a Forward-Reverse library
                        // If the mem is in the FWD strand it is the position of the last character in the read
                        // If the mem is in the REV strand it is the position of the first character in the read
    std::vector<size_t> occs; // List of occurrences of the MEM

    mem_t(size_t p, size_t l, size_t i)
    {
        pos = p;  // Position in the reference
        len = l;  // Length of the MEM
        idx = i;  // Position in the read
    }
    mem_t(size_t p, size_t l, size_t i, size_t m, size_t r)
    {
        pos = p;  // Position in the reference
        len = l;  // Length of the MEM
        idx = i;  // Position in the read
        mate = m; // Left mate (0) or Right mate (1)
        rpos = r; // Position in the read for chaining
    }
} mem_t;

typedef struct{
    ll score = 0;
    size_t mate = 2;
    bool paired = false;
    std::vector<size_t> anchors;
} chain_t;

typedef struct{
    ll G = LLONG_MAX;
    ll max_dist_x = 500;//LLONG_MAX;
    ll max_dist_y = 100;//LLONG_MAX;
    ll max_iter = 50;
    ll max_pred = 50;
    ll min_chain_score = 1;
    ll min_chain_length = 1;
} chain_config_t;

// Given a set of mems, find the chains
bool find_chains(
    const   std::vector<mem_t>& mems,
            std::vector< std::pair< size_t, size_t > >& anchors,
            std::vector<chain_t>& chains,
            const chain_config_t config = chain_config_t()        
    )
{
    /************* minimap2 dynamic programming for mem chaining ***************/
    /* https://github.com/lh3/minimap2/blob/master/chain.c */

    // Sort anchors
    // Lamda helper to sort the anchors
    auto cmp = [&] (std::pair<size_t, size_t> i, std::pair<size_t, size_t> j) -> bool {
    return (mems[i.first].occs[i.second] + mems[i.first].len - 1) < (mems[j.first].occs[j.second] + mems[j.first].len - 1);
    };

    // TODO: improve this initialization
    size_t tot_mem_length = 0;

    // std::vector< std::pair< size_t, size_t > > anchors;
    for(size_t i = 0; i < mems.size(); ++i)
    {
    for(size_t j = 0; j < mems[i].occs.size(); ++j)
        anchors.push_back(make_pair(i,j));
    tot_mem_length +=  mems[i].len * mems[i].occs.size();
    }

    float avg_mem_length = (float)tot_mem_length / anchors.size();

    std::sort(anchors.begin(),anchors.end(),cmp);

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
    const mem_t mem_i = mems[a_i.first];
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
        const mem_t mem_j = mems[a_j.first];
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

        if((mate_i == mate_j and y_j >= y_i) or max(y_d, x_d) > G)
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

    f[i] = max_f;
    p[i] = max_j;
    if( max_j >= 0 and msc[max_j] > max_f)
        msc[i] = msc[max_j];
    else
        msc[i] = max_f;
    }

    // Find the end positions of chains
    memset(t.data(),0,sizeof(size_t) * t.size());
    for(size_t i = 0; i < anchors.size(); ++i)
    if(p[i] >= 0) t[p[i]] = 1;
    
    size_t n_chains = 0;
    for(size_t i = 0; i < anchors.size(); ++i)
    if(t[i] == 0 and msc[i] > min_chain_score) n_chains ++;
    
    // TODO: check if we want to report also non aligned reads
    if(n_chains == 0)
    return false;

    // TODO: replace the vector of pairs with a lambda for the sort
    std::vector<std::pair<ll,size_t>> chain_starts(n_chains); // Stores uniqe chains and their index of uniqe chains in anchors
    size_t k = 0;
    for(size_t i = 0; i < anchors.size(); ++i)
    {
        if(t[i] == 0 and  msc[i] > min_chain_score)
        {
            size_t j = i;
            while(j >= 0 and f[j] < msc[j]) j = p[j]; // Find teh peak that maximizes f
            if (j < 0) i = 0;
            chain_starts[k++] = make_pair(f[j],j);
        }
    }

    chain_starts.resize(k), chain_starts.shrink_to_fit();
    n_chains = chain_starts.size();
    std::sort(chain_starts.begin(), chain_starts.end(), std::greater<std::pair<ll,size_t>>());
    
    // std::vector<std::pair<size_t, std::vector<size_t>>> chains;

    // Backtrack
    memset(t.data(),0,sizeof(size_t) * t.size());
    for(size_t i = 0; i < n_chains; ++i)
    {
        ll j = chain_starts[i].second;
        chain_t chain;
        chain.mate = mems[anchors[j].first].mate;
        chain.score = chain_starts[i].first;
        do {
            size_t tmp_mate = mems[anchors[j].first].mate;
            chain.paired = chain.paired or (chain.mate != tmp_mate);
            chain.anchors.push_back(j); // stores th reverse of the chain
            t[j] = 1;
            j = p[j];
        } while (j >= 0 && t[j] == 0);
        if (j < 0) { // l - prev_l is the length of the chain
            if (chain.anchors.size() >= min_chain_length)
            chains.push_back(chain);     
        } else if (chain_starts[i].first - f[j] >= min_chain_score) { // Two chains share a common prefix
            if (chain.anchors.size() >= min_chain_length)
            chains.push_back(chain);
        }
    }

    // Lamda helper to sort the anchors
    auto chain_t_cmp = [](chain_t i, chain_t j) -> bool
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
    p.resize(0), p.shrink_to_fit();
    t.resize(0), t.shrink_to_fit();
    f.resize(0), f.shrink_to_fit();
    msc.resize(0), msc.shrink_to_fit();

    return true;
}



#endif /* end of include guard: _CHAIN_HH */