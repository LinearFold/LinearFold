/*
 *utility.h*
 provides feature functions.

 author: Kai Zhao, Dezhong Deng
 edited by: 02/2018
*/

#ifndef FASTCKY_UTILITY_H
#define FASTCKY_UTILITY_H

#include <algorithm>
#include <cstring>
#include <assert.h>

#include "shared.h" // lhuang
#include "feature_weight.h"

//#define HELIX_STACKING_OLD(x, y, z, w) (_helix_stacking[GET_ACGU_NUM(x)][GET_ACGU_NUM(y)][GET_ACGU_NUM(z)][GET_ACGU_NUM(w)])

bool _helix_stacking[NOTON][NOTON][NOTON][NOTON];
double cache_single[SINGLE_MAX_LEN+1][SINGLE_MAX_LEN+1];

void initialize_cachesingle()
{
    memset(cache_single, 0, sizeof(cache_single));
    for (int l1 = SINGLE_MIN_LEN; l1 <= SINGLE_MAX_LEN; l1 ++)
        for (int l2 = SINGLE_MIN_LEN; l2 <= SINGLE_MAX_LEN; l2 ++)
        {
            if (l1 == 0 && l2 == 0)
                continue;

                // bulge
            else if (l1 == 0)
                cache_single[l1][l2] += bulge_length[l2];
            else if (l2 == 0)
                cache_single[l1][l2] += bulge_length[l1];
            else
            {

                // internal
                cache_single[l1][l2] += internal_length[std::min(l1+l2, INTERNAL_MAX_LEN)];

                // internal explicit
                if (l1 <= EXPLICIT_MAX_LEN && l2 <= EXPLICIT_MAX_LEN)
                    cache_single[l1][l2] +=
                            internal_explicit[l1<=l2 ? l1*EXPLICIT_MAX_LEN+l2 : l2*EXPLICIT_MAX_LEN+l1];

                // internal symmetry
                if (l1 == l2)
                    cache_single[l1][l2] += internal_symmetric_length[std::min(l1, SYMMETRIC_MAX_LEN)];

                else {  // internal asymmetry
                    int diff = l1 - l2; if (diff < 0) diff = -diff;
                    cache_single[l1][l2] += internal_asymmetry[std::min(diff, ASYMMETRY_MAX_LEN)];
                }
            }
        }
    return;
}

// ------------- nucs based scores -------------

// parameters: nucs[i], nucs[j]
inline double base_pair_score(int nuci, int nucj) {
  return base_pair[nucj*NOTON + nuci];
}

// parameters: nucs[i], nucs[i+1], nucs[j-1], nucs[j]
inline double helix_stacking_score(int nuci, int nuci1, int nucj_1, int nucj) {
  return helix_stacking[nuci*NOTONT + nucj*NOTOND + nuci1*NOTON + nucj_1];
}

// parameters: nucs[i], nucs[j]
inline double helix_closing_score(int nuci, int nucj) {
    return helix_closing[nuci*NOTON + nucj];
}

// parameters: nucs[i], nucs[i+1], nucs[j-1], nucs[j]
inline double terminal_mismatch_score(int nuci, int nuci1, int nucj_1, int nucj) {
    return terminal_mismatch[nuci*NOTONT+nucj*NOTOND + nuci1*NOTON + nucj_1];
}



// parameter: nucs[i]
inline double bulge_nuc_score(int nuci) {
    return bulge_0x1_nucleotides[nuci];
}

// parameters: nucs[i], nucs[j]
inline double internal_nuc_score(int nuci, int nucj) {
  return internal_1x1_nucleotides[nuci*NOTON + nucj];
}

// parameters: nucs[i], nucs[i+1], nucs[j]
inline double dangle_left_score(int nuci, int nuci1, int nucj) {
    return dangle_left[nuci*NOTOND + nucj*NOTON + nuci1];
}

// parameters: nucs[i], nucs[j-1], nucs[j]
inline double dangle_right_score(int nuci, int nucj_1, int nucj) {
    return dangle_right[nuci*NOTOND + nucj*NOTON + nucj_1];
}



// ------------- length based scores -------------

inline double hairpin_score(int i, int j) {
    return hairpin_length[std::min(j-i-1, HAIRPIN_MAX_LEN)];
}

inline double internal_length_score(int l) {
    return internal_length[std::min(l, INTERNAL_MAX_LEN)];
}

inline double internal_explicit_score(int l1, int l2){
    int l1_ = std::min(l1, EXPLICIT_MAX_LEN);
    int l2_ = std::min(l2, EXPLICIT_MAX_LEN);
    return internal_explicit[l1_<=l2_ ? l1_*NOTON+l2_ : l2_*NOTON+l1_];
}

inline double internal_sym_score(int l) {
    return internal_symmetric_length[std::min(l, SYMMETRIC_MAX_LEN)];
}

inline double internal_asym_score(int l1, int l2)
{
    int diff = l1 - l2; if (diff < 0) diff = -diff;
    return internal_asymmetry[std::min(diff, ASYMMETRY_MAX_LEN)];
}

inline double bulge_length_score(int l){
    return bulge_length[std::min(l, BULGE_MAX_LEN)];
}

inline double hairpin_at_least_score(int l) {
    return hairpin_length_at_least[std::min(l, HAIRPIN_MAX_LEN)];
}

inline double buldge_length_at_least_score(int l) {
    return bulge_length_at_least[std::min(l, BULGE_MAX_LEN)];
}

inline double internal_length_at_least_score(int l) {
    return internal_length_at_least[std::min(l, INTERNAL_MAX_LEN)];
}


//-----------------------------------------------------
inline double score_junction_A(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int len) {
    return helix_closing_score(nuci, nucj) +
            (i < len - 1 ? dangle_left_score(nuci, nuci1, nucj) : 0) +
            (j > 0 ? dangle_right_score(nuci, nucj_1, nucj) : 0);
}

inline double score_junction_B(int i, int j, int nuci, int nuci1, int nucj_1, int nucj) {
    return helix_closing_score(nuci, nucj) + terminal_mismatch_score(nuci, nuci1, nucj_1, nucj);
}

inline double score_hairpin_length(int len) {
  return hairpin_length[std::min(len, HAIRPIN_MAX_LEN)];
}

inline double score_hairpin(int i, int j, int nuci, int nuci1, int nucj_1, int nucj) {
    return hairpin_length[std::min(j-i-1, HAIRPIN_MAX_LEN)] +
            score_junction_B(i, j, nuci, nuci1, nucj_1, nucj);
}

inline double score_helix(int nuci, int nuci1, int nucj_1, int nucj) {
    return helix_stacking_score(nuci, nuci1, nucj_1, nucj) + base_pair_score(nuci1, nucj_1);
}

inline double score_single_nuc(int i, int j, int p, int q, int nucp_1, int nucq1) {
    int l1 = p-i-1, l2=j-q-1;
    if (l1==0 && l2==1) return bulge_nuc_score(nucq1);
    if (l1==1 && l2==0) return bulge_nuc_score(nucp_1);
    if (l1==1 && l2==1) return internal_nuc_score(nucp_1, nucq1);
    return 0;
}

inline double score_single(int i, int j, int p, int q, int len,
                           int nuci, int nuci1, int nucj_1, int nucj,
                           int nucp_1, int nucp, int nucq, int nucq1) {
    int l1 = p-i-1, l2=j-q-1;
    return cache_single[l1][l2] +
           base_pair_score(nucp, nucq) +
           score_junction_B(i, j, nuci, nuci1, nucj_1, nucj) +
           score_junction_B(q, p, nucq, nucq1, nucp_1, nucp) +
           score_single_nuc(i, j, p, q, nucp_1, nucq1);
}

// score_single without socre_junction_B
inline double score_single_without_junctionB(int i, int j, int p, int q,
                           int nucp_1, int nucp, int nucq, int nucq1) {
    int l1 = p-i-1, l2=j-q-1;
    return cache_single[l1][l2] +
           base_pair_score(nucp, nucq) +
           score_single_nuc(i, j, p, q, nucp_1, nucq1);
}

inline double score_multi(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int len) {
    return score_junction_A(i, j, nuci, nuci1, nucj_1, nucj, len) +
           multi_paired + multi_base;
}

inline double score_multi_unpaired(int i, int j) {
    return (j-i+1) * multi_unpaired;
}

inline double score_M1(int i, int j, int k, int nuci_1, int nuci, int nuck, int nuck1, int len) {
    return score_junction_A(k, i, nuck, nuck1, nuci_1, nuci, len) +
           score_multi_unpaired(k+1, j) + base_pair_score(nuci, nuck) + multi_paired;
}

inline double score_external_paired(int i, int j, int nuci_1, int nuci, int nucj, int nucj1, int len) {
    return score_junction_A(j, i, nucj, nucj1, nuci_1, nuci, len) +
           external_paired + base_pair_score(nuci, nucj);
}

inline double score_external_unpaired(int i, int j) {
    return (j-i+1) * external_unpaired;
}

#endif //FASTCKY_UTILITY_H
