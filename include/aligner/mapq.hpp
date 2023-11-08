/* mapq - Compute the mapping quality of the alignment.
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
   \file mapq.hpp
   \brief mapq.hpp Compute the mapping quality of the alignment.
   \author Massimiliano Rossi
   \date 19/10/2021
*/

#ifndef _MAPQ_HH
#define _MAPQ_HH


// From https://github.com/BenLangmead/bowtie2/blob/4512b199768e562e8627ffdfd9253affc96f6fc6/unique.cpp
// There is no valid second-best alignment and the best alignment has a
// perfect score.
const uint32_t unp_nosec_perf = 44;

// There is no valid second-best alignment.  We stratify the alignment
// score of the best alignment into 10 bins.
const uint32_t unp_nosec[11] = {
43, 42, 41, 36, 32, 27, 20, 11, 4, 1, 0
};

// The best alignment has a perfect score, and we stratify the distance
// between best and second-best alignment scores into 10 bins.
const uint32_t unp_sec_perf[11] = {
2, 16, 23, 30, 31, 32, 34, 36, 38, 40, 42
};

// The best alignment has a non-perfect score, and we stratify both by best
// alignment score (specifically, the maximum score minus the best "best")
// and by the distance between the best and second-best alignment scores
// ("difference").  Each is stratified into 10 bins.  Each row is a
// difference (smaller elts = smaller differences) and each column is a
// best score (smaller elts = higher best alignment scores).
const uint32_t unp_sec[11][11] = {
{  2,  2,  2,  1,  1, 0, 0, 0, 0, 0, 0},
{ 20, 14,  7,  3,  2, 1, 0, 0, 0, 0, 0},
{ 20, 16, 10,  6,  3, 1, 0, 0, 0, 0, 0},
{ 20, 17, 13,  9,  3, 1, 1, 0, 0, 0, 0},
{ 21, 19, 15,  9,  5, 2, 2, 0, 0, 0, 0},
{ 22, 21, 16, 11, 10, 5, 0, 0, 0, 0, 0},
{ 23, 22, 19, 16, 11, 0, 0, 0, 0, 0, 0},
{ 24, 25, 21, 30,  0, 0, 0, 0, 0, 0, 0},
{ 30, 26, 29,  0,  0, 0, 0, 0, 0, 0, 0},
{ 30, 27,  0,  0,  0, 0, 0, 0, 0, 0, 0},
{ 30,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0},
};

//
// Paired mapping quality:
//

// There is no valid second-best alignment and the best alignment has a
// perfect score.
const uint32_t pair_nosec_perf = 44;


/*!
    Compute the mapping quality of the alignment
    Inspired from https://github.com/BenLangmead/bowtie2/blob/4512b199768e562e8627ffdfd9253affc96f6fc6/unique.h
  */
size_t compute_mapq(
    const int32_t score,      // Best alignment score
    const int32_t score2,     // Second best alignemt score
    const int32_t min_score,  // Minimum alignemt score
    const int32_t max_score   // Maximum alignment score
    )
    {
    if(score2 > score)
        return 0;
    assert(score2 <= score);
    int32_t best = max_score - score;
    size_t best_bin = (size_t)((double)best * (10.0 / (double)(max_score - min_score)) + 0.5);
    if(score2 >= min_score)
    {
        int32_t diff = score - score2;
        size_t diff_bin = (size_t)((double)diff * (10.0 / (double)(max_score - min_score)) + 0.5);
        if(best == max_score)
        return unp_sec_perf[best_bin];
        else
        return unp_sec[diff_bin][best_bin];
    }
    else
    {
        if(best == max_score)
        return unp_nosec_perf;
        else
        return unp_nosec[best_bin];
    }
}

// /*!
//     Compute the mapping quality of the alignment
//     Inspired from https://github.com/BenLangmead/bowtie2/blob/4512b199768e562e8627ffdfd9253affc96f6fc6/unique.h
//   */
// size_t compute_mapq(
//     const int32_t score,      // Best alignment score
//     const int32_t score2,     // Second best alignemt score
//     const int32_t min_score,  // Minimum alignemt score
//     const int32_t max_score   // Maximum alignment score
//     )
//     {
//     assert(score2 <= score);
//     int32_t best = max_score - score;
//     size_t best_bin = (size_t)((double)best * (10.0 / (double)(max_score - min_score)) + 0.5);
//     if(score2 >= min_score)
//     {
//         int32_t diff = score - score2;
//         size_t diff_bin = (size_t)((double)diff * (10.0 / (double)(max_score - min_score)) + 0.5);
//         if(best == max_score)
//         return unp_sec_perf[best_bin];
//         else
//         return unp_sec[diff_bin][best_bin];
//     }
//     else
//     {
//         if(best == max_score)
//         return unp_nosec_perf;
//         else
//         return unp_nosec[best_bin];
//     }
// }


// o->mapQ_coef_len = 50; o->mapQ_coef_fac = log(o->mapQ_coef_len);
/*!
    Compute the mapping quality of the alignment
    Inspired from https://github.com/lh3/bwa/blob/34374c56139b7a08b3ef7dae38aca36f43f3cdd1/bwamem.c#L981
  */
#define MEM_MAPQ_COEF 30.0
#define raw_mapq(diff, a) ((int)(6.02 * (diff) / (a) + .499))

size_t compute_mapq_se_bwa(
    const int32_t score,        // Best alignment score
    const int32_t score2,       // Second best alignemt score
    const int32_t rlen,         // Length of the reference in the alignment
    const int32_t qlen,             // Length of the read in the alignment
    const int32_t min_seed_length,  // Minimum seed lenght
    const int32_t match_score,      // Match score
    const int32_t mismatch_score,   // Mismatch score
    const double mapq_coeff_len,    // MAPQ coefficient length
    const int32_t mapq_coeff_fac,   // MAPQ coefficient factor
    const int32_t sub_n,            // Number of sub-optimal alignments
    const int32_t seed_cov,         // Length of the region covered by seeds (https://github.com/lh3/bwa/blob/0747fcc09d96ff44ce555f0c258d0f9762c20611/bwamem.c#L800)
    const double frac_rep          // Length of the region covered by seeds (https://github.com/lh3/bwa/blob/0747fcc09d96ff44ce555f0c258d0f9762c20611/bwamem.c#L291)
    )
    {

    int32_t mapq = 0; // Mapping quality 
    int32_t l = std::max(rlen, qlen); // Length of the longest segment
    int32_t sub = score2 ? score2: min_seed_length * match_score;
    if( sub >= score) return mapq;

    double identity = 1. - (double)(l * match_score - score) / (match_score + mismatch_score) / l;
	if ( score == 0 ) {
		mapq = 0;
	} else if (mapq_coeff_len > 0) {
		double tmp;
		tmp = l < mapq_coeff_len? 1. : mapq_coeff_fac / log(l);
		tmp *= identity * identity;
		mapq = (int)(6.02 * (score - sub) / match_score * tmp * tmp + .499);
	} else {
		mapq = (int)(MEM_MAPQ_COEF * (1. - (double)sub / score) * log(seed_cov) + .499);
		mapq = identity < 0.95? (int)(mapq * identity * identity + .499) : mapq;
	}
	if (sub_n > 0) mapq -= (int)(4.343 * log(sub_n+1) + .499);
	if (mapq > 60) mapq = 60;
	if (mapq < 0) mapq = 0;
	mapq = (int)(mapq * (1. - frac_rep) + .499);
	return mapq;
}

size_t compute_mapq_pe_bwa(
    const int32_t score,       // Best alignment score
    const int32_t score2,      // Second best alignemt score
    const int32_t score_un,    // Score if unpaired
    const int32_t match_score, // Match score
    const int32_t sub_n,       // Number of sub-optimal alignments
    const double frac_rep_m1,  // Length of the region covered by seeds of mate 1(https://github.com/lh3/bwa/blob/0747fcc09d96ff44ce555f0c258d0f9762c20611/bwamem.c#L291)
    const double frac_rep_m2,  // Length of the region covered by seeds of mate 2(https://github.com/lh3/bwa/blob/0747fcc09d96ff44ce555f0c258d0f9762c20611/bwamem.c#L291)
    const int32_t score_m1,    // Best score m1
    const int32_t score_m2,    // Best score m2
    const int32_t score2_m1,   // Second best score m1
    const int32_t score2_m2,   // Second best score m2
    size_t &mapq_m1,           // MAPQ mate1
    size_t &mapq_m2            // MAPQ mate2
)
{
    assert(score2 <= score);
    int32_t mapq = 0; // Mapping quality 

    int32_t sub = std::max(score2, score_un);
    mapq = raw_mapq(score - sub, match_score);
    if (sub_n > 0) mapq -= (int)(4.343 * log(sub_n+1) + .499);
    if (mapq < 0) mapq = 0;
    if (mapq > 60) mapq = 60;
    mapq = (int)(mapq * (1. - .5 * (frac_rep_m1 + frac_rep_m2)) + .499);
    
    if (score > score_un) { // paired alignment


        mapq_m1 = mapq_m1 > mapq? mapq_m1 : mapq < mapq_m1 + 40? mapq : mapq_m1 + 40;
        mapq_m2 = mapq_m2 > mapq? mapq_m2 : mapq < mapq_m2 + 40? mapq : mapq_m2 + 40;
        // // cap at the tandem repeat score
        mapq_m1 = mapq_m1 < raw_mapq(score_m1 - score2_m1, match_score) ? mapq_m1 : raw_mapq(score_m1 - score2_m1, match_score);
        mapq_m2 = mapq_m2 < raw_mapq(score_m2 - score2_m2, match_score)? mapq_m2 : raw_mapq(score_m2 - score2_m2, match_score);
    }

	return mapq;
}


#endif /* end of include guard: _MAPQ_HH */