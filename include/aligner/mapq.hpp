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

#endif /* end of include guard: _MAPQ_HH */