/*
 *LinearFold.cpp*
 The main code for LinearFold: Linear-Time RNA Structure Prediction.

 author: Kai Zhao, Dezhong Deng, He Zhang
 edited by: 11/2018
*/

#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>
#include <set>

#include "LinearFold.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"
#include "LinearFoldEval.cpp" // adding eval mode

#define SPECIAL_HP

using namespace std;


#ifdef lv
    bool comparefunc(std::pair<int,State> a, std::pair<int,State> b) {
        return a.first > b.first;
    }

    void BeamCKYParser::sort_keys(std::unordered_map<int, State> &map, std::vector<std::pair<int,State>> &sorted_keys) {
        sorted_keys.clear();
        for(auto &kv : map) {
            sorted_keys.push_back(kv);
        }
        sort(sorted_keys.begin(), sorted_keys.end(), comparefunc);    
    }
#endif

void BeamCKYParser::get_parentheses(char* result, string& seq) {
    memset(result, '.', seq_length);
    result[seq_length] = 0;

    stack<tuple<int, int, State>> stk;
    stk.push(make_tuple(0, seq_length-1, bestC[seq_length-1]));

    if(is_verbose){
            printf(">verbose\n");
    }
    // verbose stuff
    vector<pair<int,int>> multi_todo;
    unordered_map<int,int> mbp; // multi bp
    double total_energy = .0;
    double external_energy = .0;

    while ( !stk.empty() ) {
        tuple<int, int, State> top = stk.top();
        int i = get<0>(top), j = get<1>(top);
        State& state = get<2>(top);
        stk.pop();

        int k, p, q;

        switch (state.manner) {
            case MANNER_H:
                // this state should not be traced
                break;
            case MANNER_HAIRPIN:
                {
                    result[i] = '(';
                    result[j] = ')';
                    if(is_verbose) {
                        int tetra_hex_tri = -1;
                        if (j-i-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[i];
                        else if (j-i-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[i];
                        else if (j-i-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[i];
                        int nuci = nucs[i], nucj = nucs[j];
                        int nuci1 = (i + 1) < seq_length ? nucs[i + 1] : -1;
                        int nucj_1 = (j - 1) > -1 ? nucs[j - 1] : -1;

                        value_type newscore = - v_score_hairpin(i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri);
                        printf("Hairpin loop ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], newscore / -100.0);
                        total_energy += newscore;
                    }
                }
                break;
            case MANNER_SINGLE:
                {
                    result[i] = '(';
                    result[j] = ')';
                    p = i + state.trace.paddings.l1;
                    q = j - state.trace.paddings.l2;
                    stk.push(make_tuple(p, q, bestP[q][p]));
                    if(is_verbose) {
                        int nuci = nucs[i], nuci1 = nucs[i+1], nucj_1 = nucs[j-1], nucj = nucs[j];
                        int nucp_1 = nucs[p-1], nucp = nucs[p], nucq = nucs[q], nucq1 = nucs[q+1];

                        value_type newscore = -v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj,
                                                          nucp_1, nucp, nucq, nucq1);
                        printf("Interior loop ( %d, %d) %c%c; ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], p+1, q+1, seq[p],seq[q], newscore / -100.0);
                        total_energy += newscore;
                    }
                }
                break;
            case MANNER_HELIX:
                {
                    result[i] = '(';
                    result[j] = ')';
                    stk.push(make_tuple(i+1, j-1, bestP[j-1][i+1]));
                    if(is_verbose){
                        p = i + 1;
                        q = j - 1;
                        int nuci = nucs[i], nuci1 = nucs[i+1], nucj_1 = nucs[j-1], nucj = nucs[j];
                        int nucp_1 = nucs[p-1], nucp = nucs[p], nucq = nucs[q], nucq1 = nucs[q+1];

                        value_type newscore = -v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj,
                                                          nucp_1, nucp, nucq, nucq1);
                        printf("Interior loop ( %d, %d) %c%c; ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], p+1, q+1, seq[p],seq[q], newscore / -100.0);
                        total_energy += newscore;
                    }
                }
                break;
            case MANNER_MULTI: 
                p = i + state.trace.paddings.l1;
                q = j - state.trace.paddings.l2;
                stk.push(make_tuple(p, q, bestM2[q][p]));
                break;
            case MANNER_MULTI_eq_MULTI_plus_U:
                p = i + state.trace.paddings.l1;
                q = j - state.trace.paddings.l2;
                stk.push(make_tuple(p, q, bestM2[q][p]));
                break;
            case MANNER_P_eq_MULTI:
                result[i] = '(';
                result[j] = ')';
                stk.push(make_tuple(i, j, bestMulti[j][i]));
                if(is_verbose) {
                    multi_todo.push_back(make_pair(i,j));
                }
                break;
            case MANNER_M2_eq_M_plus_P:
                k = state.trace.split;
                stk.push(make_tuple(i, k, bestM[k][i]));
                stk.push(make_tuple(k+1, j, bestP[j][k+1]));
                if(is_verbose)
                    mbp[k+1] = j;
                break;
            case MANNER_M_eq_M2:
                stk.push(make_tuple(i, j, bestM2[j][i]));
                break;
            case MANNER_M_eq_M_plus_U:
                stk.push(make_tuple(i, j-1, bestM[j-1][i]));
                break;
            case MANNER_M_eq_P:
                stk.push(make_tuple(i, j, bestP[j][i]));
                if(is_verbose)
                    mbp[i] = j;
                break;
            case MANNER_C_eq_C_plus_U:
                k = j - 1;
                if (k != -1)
                    stk.push(make_tuple(0, k, bestC[k]));
                if (is_verbose) 
                    external_energy += - v_score_external_unpaired(0, 0); // zero at this moment
                break;
            case MANNER_C_eq_C_plus_P:
                {
                    k = state.trace.split;
                    if (k != -1) {
                        stk.push(make_tuple(0, k, bestC[k]));
                        stk.push(make_tuple(k+1, j, bestP[j][k+1]));
                    }
                    else {
                        stk.push(make_tuple(i, j, bestP[j][i]));
                    }
                    if (is_verbose) {
                        int nuck = k > -1 ? nucs[k] : -1;
                        int nuck1 = nucs[k+1], nucj = nucs[j];
                        int nucj1 = (j + 1) < seq_length ? nucs[j + 1] : -1;
                        external_energy +=  - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                    nucj, nucj1, seq_length);
                    }
                }
                break;
            default:  // MANNER_NONE or other cases
                if (use_constraints){
                    printf("We can't find a valid structure for this sequence and constraint.\n");
                    printf("There are two minor restrictions in our real system:\n");
                    printf("the length of an interior loop is bounded by 30nt \n");
                    printf("(a standard limit found in most existing RNA folding software such as CONTRAfold)\n");
                    printf("so is the leftmost (50-end) unpaired segment of a multiloop (new constraint).\n");
                    exit(1);
                } 
                printf("wrong manner at %d, %d: manner %d\n", i, j, state.manner); fflush(stdout);
                assert(false);
                
        }
    }

    if(is_verbose) {
        for (auto item : multi_todo) {
            int i = item.first;
            int j = item.second;
            int nuci = nucs[i], nuci1 = nucs[i+1], nucj_1 = nucs[j-1], nucj = nucs[j];
            value_type multi_energy = - v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, seq_length);
            int num_unpaired = 0;
            for (int k=i+1; k<j; ++k) {
                if (result[k] == '.')
                    num_unpaired += 1;
                else if (result[k] == '(') {
                    int p = k, q = mbp[k];
                    int nucp_1 = nucs[p-1], nucp = nucs[p], nucq = nucs[q], nucq1 = nucs[q+1];

                    multi_energy += - v_score_M1(p, q, q, nucp_1, nucp, nucq, nucq1, seq_length);
                    k = q;
                }
            }
            multi_energy += - v_score_multi_unpaired(1, num_unpaired);

            printf("Multi loop ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], multi_energy / -100.0);
            total_energy += multi_energy;
        }

        printf("External loop : %.2f\n", external_energy / -100.0);
        total_energy += external_energy;

#ifndef lv
            printf("Energy(kcal/mol): %.2f\n", total_energy / -100.0);
#endif
    }

    return;
}

unsigned long quickselect_partition(vector<pair<value_type, int>>& scores, unsigned long lower, unsigned long upper) {
    value_type pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);

    }
    return upper;
}

// in-place quick-select
value_type quickselect(vector<pair<value_type, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}


value_type BeamCKYParser::beam_prune(std::unordered_map<int, State> &beamstep) {
    scores.clear();
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;
        value_type newscore;
        // lisiz: for _V, avoid -inf-int=+inf
        if ((k >= 0) && (bestC[k].score == VALUE_MIN)) newscore = VALUE_MIN;
        else newscore = (k >= 0 ? bestC[k].score : 0) + cand.score;
        scores.push_back(make_pair(newscore, i));
    }
    if (scores.size() <= beam) return VALUE_MIN;
    value_type threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}

void BeamCKYParser::sortM(value_type threshold,
                          std::unordered_map<int, State> &beamstep,
                          std::vector<std::pair<value_type, int>> &sorted_stepM) {
    sorted_stepM.clear();
    if (threshold == VALUE_MIN) {
        // no beam pruning before, so scores vector not usable
        for (auto &item : beamstep) {
            int i = item.first;
            State &cand = item.second;
            int k = i - 1;
            value_type newscore;
            // lisiz: constraints may cause all VALUE_MIN, sorting has no use
            if ((use_constraints) && (k >= 0) && (bestC[k].score == VALUE_MIN)) newscore = cand.score;
            else newscore = (k >= 0 ? bestC[k].score : 0) + cand.score;
            sorted_stepM.push_back(make_pair(newscore, i));
        }
    } else {
        for (auto &p : scores) {
            if (p.first >= threshold) sorted_stepM.push_back(p);
        }
    }

    sort(sorted_stepM.begin(), sorted_stepM.end(), std::greater<pair<value_type, int>>());
}


// void BeamCKYParser::prepare(unsigned len) {
void BeamCKYParser::prepare(unsigned len) {
    seq_length = len;

    bestH.clear();
    bestH.resize(seq_length);
    bestP.clear();
    bestP.resize(seq_length);
    bestM2.clear();
    bestM2.resize(seq_length);
    bestM.clear();
    bestM.resize(seq_length);
    bestC.clear();
    bestC.resize(seq_length);
    bestMulti.clear();
    bestMulti.resize(seq_length);

#ifdef is_cube_pruning
        sorted_bestM.clear();
        sorted_bestM.resize(seq_length);
#endif

    nucs.clear();
    nucs.resize(seq_length);

    scores.reserve(seq_length);

    if (use_constraints){
        allow_unpaired_position.clear();
        allow_unpaired_position.resize(seq_length);

        allow_unpaired_range.clear();
        allow_unpaired_range.resize(seq_length);
    }
}

// lisiz, constraints
bool BeamCKYParser::allow_paired(int i, int j, vector<int>* cons, char nuci, char nucj) {
    return ((*cons)[i] == -1 || (*cons)[i] == j) && ((*cons)[j] == -1 || (*cons)[j] == i) && _allowed_pairs[nuci][nucj];
}

// BeamCKYParser::DecoderResult BeamCKYParser::parse(string& seq) {

BeamCKYParser::DecoderResult BeamCKYParser::parse(string& seq, vector<int>* cons) {

    struct timeval parse_starttime, parse_endtime;

    // number of states
    unsigned long nos_H = 0, nos_P = 0, nos_M2 = 0,
            nos_M = 0, nos_C = 0, nos_Multi = 0;

    gettimeofday(&parse_starttime, NULL);

    // prepare(static_cast<unsigned>(seq.length()));
    prepare(static_cast<unsigned>(seq.length()));

    for (int i = 0; i < seq_length; ++i)
        nucs[i] = GET_ACGU_NUM(seq[i]);

    // lisiz, constraints
    if (use_constraints) {
        for (int i=0; i<seq_length; i++){
            int cons_idx = (*cons)[i];
            allow_unpaired_position[i] = cons_idx == -1 || cons_idx == -2;
            if (cons_idx > -1){
                if (!_allowed_pairs[nucs[i]][nucs[cons_idx]]){
                    printf("Constrains on non-classical base pairs (non AU, CG, GU pairs)\n");
                    exit(1);
                }
            }
        }
        int firstpair = seq_length;
        for (int i=seq_length-1; i>-1; i--){
            allow_unpaired_range[i] = firstpair;
            if ((*cons)[i] >= 0)
                firstpair = i;
        }
    }

    vector<int> next_pair[NOTON];
    {
        if (use_constraints){
            for (int nuci = 0; nuci < NOTON; ++nuci) {
                next_pair[nuci].resize(seq_length, -1);
                int next = -1;
                for (int j = seq_length-1; j >=0; --j) {
                    next_pair[nuci][j] = next;
                    if ((*cons)[j] > -2 && _allowed_pairs[nuci][nucs[j]]) next = j;
                }
            }
        } else {
            for (int nuci = 0; nuci < NOTON; ++nuci) {
                next_pair[nuci].resize(seq_length, -1);
                int next = -1;
                for (int j = seq_length-1; j >=0; --j) {
                    next_pair[nuci][j] = next;
                    if (_allowed_pairs[nuci][nucs[j]]) next = j;
                }
            }
        }
    }

#ifdef SPECIAL_HP
#ifdef lv
    v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
#else
    if (is_verbose)
        v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
#endif
#endif

    // start CKY decoding
#ifdef lv
        if(seq_length > 0) bestC[0].set(- v_score_external_unpaired(0, 0), MANNER_C_eq_C_plus_U);
        if(seq_length > 1) bestC[1].set(- v_score_external_unpaired(0, 1), MANNER_C_eq_C_plus_U);
#else
        if(seq_length > 0) bestC[0].set(score_external_unpaired(0, 0), MANNER_C_eq_C_plus_U);
        if(seq_length > 1) bestC[1].set(score_external_unpaired(0, 1), MANNER_C_eq_C_plus_U);
#endif
    ++nos_C;

    // from left to right
    for(int j = 0; j < seq_length; ++j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of H
        {
            if (beam > 0 && beamstepH.size() > beam) beam_prune(beamstepH);

            {
                // for nucj put H(j, j_next) into H[j_next]
                int jnext = next_pair[nucj][j];
                if (no_sharp_turn) while (jnext - j < 4 && jnext != -1) jnext = next_pair[nucj][jnext];

                // lisiz, constriants
                if (use_constraints){
                    if (!allow_unpaired_position[j]){
                        jnext = (*cons)[j] > j ? (*cons)[j] : -1; // lisiz: j must be left bracket, jump to the constrainted pair (j, j') directly
                    }
                    if (jnext != -1){
                        int nucjnext = nucs[jnext];
                        if (jnext > allow_unpaired_range[j] || !allow_paired(j, jnext, cons, nucj, nucjnext))  // lisiz: avoid cross constrainted brackets or unallowed pairs
                            jnext = -1;
                    }
                }

                if (jnext != -1) {
                    int nucjnext = nucs[jnext];
                    int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                    value_type newscore;

#ifdef lv
                        int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                        if (jnext-j-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[j];
                        else if (jnext-j-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[j];
                        else if (jnext-j-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[j];
#endif
                        newscore = - v_score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext, tetra_hex_tri);
#else
                        newscore = score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext);
#endif
                    // this candidate must be the best one at [j, jnext]
                    // so no need to check the score
                    update_if_better(bestH[jnext][j], newscore, MANNER_H);
                    ++ nos_H;
                }
            }

            {
                // for every state h in H[j]
                //   1. extend h(i, j) to h(i, jnext)
                //   2. generate p(i, j)
#ifdef lv
                sort_keys(beamstepH, keys);
                for (auto &item : keys) {
#else
                for (auto &item : beamstepH) {
#endif
                    int i = item.first;
                    // printf("%d\n", i);
                    State &state = item.second;
                    int nuci = nucs[i];
                    int jnext = next_pair[nuci][j];

                    // 2. generate p(i, j)
                    // lisiz, change the order because of the constriants
                    {
                        update_if_better(beamstepP[i], state.score, MANNER_HAIRPIN);
                        ++ nos_P;
                    }

                    // lisiz, constraints
                    if (jnext != -1 && use_constraints){
                        int nucjnext = nucs[jnext];
                        if (jnext > allow_unpaired_range[i] || !allow_paired(i, jnext, cons, nuci, nucjnext))
                            continue;
                    }

                    if (jnext != -1) {
                        int nuci1 = (i + 1) < seq_length ? nucs[i + 1] : -1;
                        int nucjnext = nucs[jnext];
                        int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                        // 1. extend h(i, j) to h(i, jnext)
                        value_type newscore;

#ifdef lv
                            int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                            if (jnext-i-1 == 4) // 6:tetra
                                tetra_hex_tri = if_tetraloops[i];
                            else if (jnext-i-1 == 6) // 8:hexa
                                tetra_hex_tri = if_hexaloops[i];
                            else if (jnext-i-1 == 3) // 5:tri
                                tetra_hex_tri = if_triloops[i];
#endif
                            newscore = - v_score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);
#else
                            newscore = score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext);
#endif
                        // this candidate must be the best one at [i, jnext]
                        // so no need to check the score
                        update_if_better(bestH[jnext][i], newscore, MANNER_H);
                        ++nos_H;
                    }
                }
            }
        }
        if (j == 0) continue;

        // beam of Multi
        {
            if (beam > 0 && beamstepMulti.size() > beam) beam_prune(beamstepMulti);

            // for every state in Multi[j]
            //   1. extend (i, j) to (i, jnext)
            //   2. generate P (i, j)
#ifdef lv
            sort_keys(beamstepMulti, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepMulti) {
#endif
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 2. generate P (i, j)
                // lisiz, change the order because of the constraits
                {
                    value_type newscore;
#ifdef lv
                        newscore = state.score - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
#else
                        newscore = state.score + score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
#endif
                    update_if_better(beamstepP[i], newscore, MANNER_P_eq_MULTI);
                    ++ nos_P;
                }

                // lisiz cnstriants
                if (jnext != -1 && use_constraints){
                    int nucjnext = nucs[jnext];
                    if (jnext > allow_unpaired_range[j] || !allow_paired(i, jnext, cons, nuci, nucjnext))
                        continue;
                }

                // 1. extend (i, j) to (i, jnext)
                {
                    char new_l1 = state.trace.paddings.l1;
                    int new_l2 = state.trace.paddings.l2 + jnext - j;
                    // if (jnext != -1 && new_l1 + new_l2 <= SINGLE_MAX_LEN) {
                    if (jnext != -1) {
                        // 1. extend (i, j) to (i, jnext)
                        value_type newscore;
#ifdef lv
                            newscore = state.score - v_score_multi_unpaired(j, jnext - 1);
#else
                            newscore = state.score + score_multi_unpaired(j, jnext - 1);
#endif
                        // this candidate must be the best one at [i, jnext]
                        // so no need to check the score
                        update_if_better(bestMulti[jnext][i], newscore, MANNER_MULTI_eq_MULTI_plus_U,
                                         new_l1,
                                         new_l2
                        );
                        ++nos_Multi;
                    }
                }
            }
        }

        // beam of P
        {
            if (beam > 0 && beamstepP.size() > beam) beam_prune(beamstepP);

            // for every state in P[j]
            //   1. generate new helix/bulge
            //   2. M = P
            //   3. M2 = M + P
            //   4. C = C + P
#ifdef is_cube_pruning
            bool use_cube_pruning = beam > MIN_CUBE_PRUNING_SIZE
                                    && beamstepP.size() > MIN_CUBE_PRUNING_SIZE;
#else
            bool use_cube_pruning = false;
#endif               

#ifdef lv
            sort_keys(beamstepP, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepP) {
#endif
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                // 2. M = P
                if(i > 0 && j < seq_length-1){
                    value_type newscore;
#ifdef lv
                        newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#else
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#endif
                    update_if_better(beamstepM[i], newscore, MANNER_M_eq_P);
                    ++ nos_M;
                }
                //printf(" M = P at %d\n", j); fflush(stdout);

                // 3. M2 = M + P
                if(!use_cube_pruning) {
                    int k = i - 1;
                    if ( k > 0 && !bestM[k].empty()) {
                        value_type M1_score;
#ifdef lv
                            M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#else
                            M1_score = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#endif
                        // candidate list
                        auto bestM2_iter = beamstepM2.find(i);
#ifndef is_candidate_list
                        for (auto &m : bestM[k]) {
                                int newi = m.first;
                                // eq. to first convert P to M1, then M2/M = M + M1
                                value_type newscore = M1_score + m.second.score;
                                update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                                //update_if_better(bestM[j][newi], newscore, MANNER_M_eq_M_plus_P, k);
                                ++nos_M2;
                                //++nos_M;
                            }
#else
                        if (bestM2_iter==beamstepM2.end() || M1_score > bestM2_iter->second.score) {
                            for (auto &m : bestM[k]) {
                                int newi = m.first;
                                // eq. to first convert P to M1, then M2/M = M + M1
                                value_type newscore = M1_score + m.second.score;
                                update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                                //update_if_better(bestM[j][newi], newscore, MANNER_M_eq_M_plus_P, k);
                                ++nos_M2;
                                //++nos_M;
                            }
                        }
#endif
                    }
                }
                //printf(" M/M2 = M + P at %d\n", j); fflush(stdout);

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                      State& prefix_C = bestC[k];
                      if (prefix_C.manner != MANNER_NONE) {
                        int nuck = nuci_1;
                        int nuck1 = nuci;
                        value_type newscore;
#ifdef lv
                            newscore = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                 nucj, nucj1, seq_length) +
                                prefix_C.score + state.score;
#else
                            newscore = score_external_paired(k+1, j, nuck, nuck1,
                                                             nucj, nucj1, seq_length) +
                                prefix_C.score + state.score;
#endif
                        update_if_better(beamstepC, newscore, MANNER_C_eq_C_plus_P, k);
                        ++ nos_C;
                      }
                    } else {
                        value_type newscore;
#ifdef lv
                            newscore = - v_score_external_paired(0, j, -1, nucs[0],
                                                                 nucj, nucj1, seq_length) +
                                state.score;
#else
                            newscore = score_external_paired(0, j, -1, nucs[0],
                                                             nucj, nucj1, seq_length) +
                                state.score;
#endif
                        update_if_better(beamstepC, newscore, MANNER_C_eq_C_plus_P, -1);
                        ++ nos_C;
                    }
                }
                //printf(" C = C + P at %d\n", j); fflush(stdout);

                // 1. generate new helix / single_branch
                // new state is of shape p..i..j..q
                if (i >0 && j<seq_length-1) {
                    value_type precomputed;
#ifdef lv
                        precomputed = 0;
#else
                        precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci);
#endif
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; // hzhang: move here
                        int q = next_pair[nucp][j];

                        // lisiz constraints
                        if (use_constraints){
                            if (p < i-1 && !allow_unpaired_position[p+1]) // lisiz: if p+1 must be paired, break
                                break;
                            if (!allow_unpaired_position[p]){             // lisiz: if p must be paired, p must be left bracket
                                q = (*cons)[p];
                                if (q < p) break;
                            }
                        }

                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];

                            // lisiz constraints
                            if (use_constraints){
                                if (q>j+1 && q > allow_unpaired_range[j])  // lisiz: if q-1 must be paired, break
                                    break;
                                if (!allow_paired(p, q, cons, nucp, nucq)) // lisiz: if p q are )(, break
                                    break;
                            }

                            // int nucp_1 = nucs[p - 1]; // hzhang: no need
                            // int nucp1 = nucs[p + 1]; // hzhang: move to outside of while loop
                            int nucq_1 = nucs[q - 1];
                            // int nucq1 = nucs[q + 1]; // hzhang: no need
                            // int nuci1 = nucs[i + 1]; // hzhang: no need
                            // int nucj_1 = nucs[j - 1]; // hzhang: no need
                            // int nucq = nucs[q];
                            // int nucp1 = nucs[p + 1];
                            // int nucq_1 = nucs[q - 1];

                            if (p == i - 1 && q == j + 1) {
                                // helix
                                value_type newscore;
#ifdef lv
                                    newscore = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                             nuci_1, nuci, nucj, nucj1)
                                        + state.score;
#else
                                    newscore = score_helix(nucp, nucp1, nucq_1, nucq) + state.score;
#endif
                                update_if_better(bestP[q][p], newscore, MANNER_HELIX);
                                ++nos_P;
                            } else {
                                // single branch
                                value_type newscore;
#ifdef lv
                                    newscore = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                   nuci_1, nuci, nucj, nucj1)
                                        + state.score;
#else
                                    newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
                                        precomputed +
                                        score_single_without_junctionB(p, q, i, j,
                                                                       nuci_1, nuci, nucj, nucj1) +
                                        state.score;
#endif
                                update_if_better(bestP[q][p], newscore, MANNER_SINGLE,
                                                 static_cast<char>(i - p),
                                                 q - j);
                                ++nos_P;
                            }
                            q = next_pair[nucp][q];
                        }
                    }
                }
                //printf(" helix / single at %d\n", j); fflush(stdout);
            }

            if (use_cube_pruning) {
                // 3. M2 = M + P with cube pruning
                vector<int> valid_Ps;
                vector<value_type> M1_scores;
#ifdef lv
            sort_keys(beamstepP, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepP) {
#endif
                    int i = item.first;
                    State &state = item.second;
                    int nuci = nucs[i];
                    int nuci_1 = (i - 1 > -1) ? nucs[i - 1] : -1;
                    int k = i - 1;

                    // group candidate Ps
                    if (k > 0 && !bestM[k].empty()) {
                        assert(bestM[k].size() == sorted_bestM[k].size());
                        value_type M1_score;
#ifdef lv
                            M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#else
                            M1_score = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#endif
                        auto bestM2_iter = beamstepM2.find(i);
#ifndef is_candidate_list
                        valid_Ps.push_back(i);
                        M1_scores.push_back(M1_score);
#else
                        if (bestM2_iter == beamstepM2.end() || M1_score > bestM2_iter->second.score) {
                            valid_Ps.push_back(i);
                            M1_scores.push_back(M1_score);
                        }
#endif
                    }
                }

                // build max heap
                // heap is of form (heuristic score, (index of i in valid_Ps, index of M in bestM[i-1]))
                vector<pair<value_type, pair<int, int>>> heap;
                for (int p = 0; p < valid_Ps.size(); ++p) {
                    int i = valid_Ps[p];
                    int k = i - 1;
                    heap.push_back(make_pair(M1_scores[p] + sorted_bestM[k][0].first,
                                             make_pair(p, 0)
                    ));
                    push_heap(heap.begin(), heap.end());
                }

                // start cube pruning
                // stop after beam size M2 states being filled
                int filled = 0;
                // exit when filled >= beam and current score < prev score
                value_type prev_score = VALUE_MIN;
                value_type current_score = VALUE_MIN;
                while ((filled < beam || current_score == prev_score) && !heap.empty()) {
                    auto &top = heap.front();
                    prev_score = current_score;
                    current_score = top.first;
                    int index_P = top.second.first;
                    int index_M = top.second.second;
                    int i = valid_Ps[top.second.first];
                    int k = i - 1;
                    int newi = sorted_bestM[k][index_M].second;
                    value_type newscore = M1_scores[index_P] + bestM[k][newi].score;
                    pop_heap(heap.begin(), heap.end());
                    heap.pop_back();

                    if (beamstepM2[newi].manner == MANNER_NONE) {
                        ++filled;
                        update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                        ++nos_M2;
                    } else {
                        assert(beamstepM2[newi].score > newscore - 1e-8);
                    }

                    ++index_M;
                    while (index_M < sorted_bestM[k].size()) {
                        // candidate_score is a heuristic score
                        value_type candidate_score = M1_scores[index_P] + sorted_bestM[k][index_M].first;
                        int candidate_newi = sorted_bestM[k][index_M].second;
                        if (beamstepM2.find(candidate_newi) == beamstepM2.end()) {
                            heap.push_back(make_pair(candidate_score,
                                                     make_pair(index_P, index_M)));
                            push_heap(heap.begin(), heap.end());
                            break;
                        } else {
                            // based on the property of cube pruning, the new score must be worse
                            // than the state already inserted
                            // so we keep iterate through the candidate list to find the next
                            // candidate
                            ++index_M;
                            assert(beamstepM2[candidate_newi].score >
                                   M1_scores[index_P] + bestM[k][candidate_newi].score - 1e-8);
                        }
                    }
                }
            }
        }
        //printf("P at %d\n", j); fflush(stdout);

        // beam of M2
        {
            if (beam > 0 && beamstepM2.size() > beam) beam_prune(beamstepM2);

            // for every state in M2[j]
            //   1. multi-loop  (by extending M2 on the left)
            //   2. M = M2
#ifdef lv
            sort_keys(beamstepM2, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepM2) {
#endif
                int i = item.first;
                State& state = item.second;

                // 2. M = M2
                {
                    update_if_better(beamstepM[i], state.score, MANNER_M_eq_M2);
                    ++ nos_M;
                }

                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp][j];

                        if (use_constraints){
                            if (p < i - 1 && !allow_unpaired_position[p+1])
                                break;
                            if (!allow_unpaired_position[p]){
                                q = (*cons)[p];
                                if (q < p) break;
                            }
                            if (q > j+1 && q > allow_unpaired_range[j])
                                continue;
                            int nucq = nucs[q];
                            if (!allow_paired(p, q, cons, nucp, nucq))
                                continue;
                        }

                        if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {
                            // the current shape is p..i M2 j ..q

                            value_type newscore;
#ifdef lv
                                newscore = - v_score_multi_unpaired(p+1, i-1) -
                                    v_score_multi_unpaired(j+1, q-1) + state.score;
#else
                                newscore = score_multi_unpaired(p+1, i-1) +
                                    score_multi_unpaired(j+1, q-1) + state.score;
#endif
                            update_if_better(bestMulti[q][p], newscore, MANNER_MULTI,
                                             static_cast<char>(i - p),
                                             q - j);
                            ++ nos_Multi;

                            //q = next_pair[nucp][q];
                        }
                    }
                }
            }
        }
        //printf("M2 at %d\n", j); fflush(stdout);

        // beam of M
        {
            value_type threshold = VALUE_MIN;
            if (beam > 0 && beamstepM.size() > beam) threshold = beam_prune(beamstepM);

#ifdef is_cube_pruning
                sortM(threshold, beamstepM, sorted_bestM[j]);
            // }
#endif

            // for every state in M[j]
            //   1. M = M + unpaired
#ifdef lv
            sort_keys(beamstepM, keys);
            for (auto &item : keys) {
#else
            for(auto& item : beamstepM) {
#endif
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
                    if (use_constraints && !allow_unpaired_position[j+1]) // if j+1 must be paired
                        continue;
                    value_type newscore;
                    // if (use_vienna)
#ifdef lv
                        newscore = - v_score_multi_unpaired(j + 1, j + 1) + state.score;
                    // else
#else
                        newscore = score_multi_unpaired(j + 1, j + 1) + state.score;
#endif
                    update_if_better(bestM[j+1][i], newscore, MANNER_M_eq_M_plus_U);
                    ++ nos_M;
                }
            }
        }
        //printf("M at %d\n", j); fflush(stdout);

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
                if (use_constraints && !allow_unpaired_position[j+1])
                        continue;
                value_type newscore;
#ifdef lv
                    newscore = -v_score_external_unpaired(j+1, j+1) + beamstepC.score;
#else
                    newscore = score_external_unpaired(j+1, j+1) + beamstepC.score;
#endif
                update_if_better(bestC[j+1], newscore, MANNER_C_eq_C_plus_U);
                ++ nos_C;
            }
        }
        //printf("C at %d\n", j); fflush(stdout);

    }  // end of for-loo j

    State& viterbi = bestC[seq_length-1];

    char result[seq_length+1];
    get_parentheses(result, seq);

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

    unsigned long nos_tot = nos_H + nos_P + nos_M2 + nos_Multi + nos_M + nos_C;
    if (is_verbose) {
        printf("Parse Time: %f len: %d score %f #states %lu H %lu P %lu M2 %lu Multi %lu M %lu C %lu\n",
               parse_elapsed_time, seq_length, double(viterbi.score), nos_tot,
               nos_H, nos_P, nos_M2, nos_Multi, nos_M, nos_C);
    }

    fflush(stdout);

    return {string(result), viterbi.score, nos_tot, parse_elapsed_time};
}

BeamCKYParser::BeamCKYParser(int beam_size,
                             bool nosharpturn,
                             bool verbose,
                             bool constraints)
    : beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      is_verbose(verbose),
      use_constraints(constraints){
#ifdef lv
        initialize();
#else
        initialize();
        initialize_cachesingle();
#endif
}


// -------------------------------------------------------------

int main(int argc, char** argv){

    int beamsize = 100;
    bool sharpturn = false;
    bool is_verbose = false;
    bool is_eval = false;
    bool is_constraints = false; // lisiz, add constraints

    if (argc > 1) {
        beamsize = atoi(argv[1]);
        sharpturn = atoi(argv[2]) == 1;
        is_verbose = atoi(argv[3]) == 1;
        is_eval = atoi(argv[4]) == 1; // adding eval mode
        is_constraints = atoi(argv[5]) == 1; // lisiz, add constraints
    }

    // variables for decoding
    int num=0, total_len = 0;
    unsigned long long total_states = 0;
    double total_score = .0;
    double total_time = .0;

    if (is_eval) {
        int lineIndex = 0;
        string seq, ref;
        bool seqflag = false;
        for (string input; getline(cin, input);) {
            // printf("Input: %s\n", input.c_str());
            if (lineIndex % 2 == 0) { // seq
                seq = input;
                if (seq.length() == 0) {
                    // printf("empty sequence!\n");
                    // seqflag = true;
                    // lineIndex ++;
                    continue;
                }

                else if (seq[0] == ';' || seq[0] == '>' || (!isalpha(seq[0]))) {
                    printf("Unrecognized sequence: %s\n", seq.c_str());
                    seqflag = true;
                    lineIndex ++;
                    continue;
                }
            }

            else { // structure
                ref = input;
                // printf("%c\n", ref[0]);
                if (ref.length() == 0) {
                    // printf("empty structure!\n");
                    // return 0;
                    continue;
                }

                else if ((ref[0] != '.') && (ref[0] != '(') && (ref[0] != ')')) {
                    printf("Unrecognized structure: %s\n", ref.c_str());
                    return 0;
                }

                if (seqflag) {
                    printf("Reference with wrong sequence!\n");
                    lineIndex ++;
                    seqflag = false;
                    continue;
                }

                if (seq.length() != ref.length()) {
                    printf("sequence length is not equal to structure length!\n");
                    lineIndex ++;
                    continue;
                } 

                // remove peudoknots
                char r;
                map<char, char> rs = { {'[', '.'}, {']', '.'}, {'{', '.'}, {'}', '.'}, {'<', '.'}, {'>', '.'} };
                replace_if(ref.begin(), ref.end(), [&](char c){ return r = rs[c]; }, r);

                // print input
                // printf("%s\n", ref.c_str());

                // call eval function;
                double MFE_energy = eval(seq, ref, is_verbose) / -100.0;
                // printf("structure energy: %.2f\n", MFE_energy/ -100.0);
                printf("%s\n", seq.c_str());
                printf("%s (%.2f)\n", input.c_str(), MFE_energy);

            }
            lineIndex ++;
        }
    }

    else {
        if (is_constraints){
            int lineIndex = 0;
            string seq, constr;
            bool seqflag, consflag;
            set<char> consSet {'?', '.', '(', ')'};
            for (string input; getline(cin, input);){
                if (input[0] == ';' || input[0] == '>') {
                    printf("%s\n", input.c_str());
                    continue;
                }
                if (lineIndex % 2 == 0) { // seq
                    seq = input;
                    
                    // check the seq
                    if (seq.length() == 0)
                    continue;
                    
                    if (!isalpha(seq[0])){
                        printf("Unrecognized sequence: %s\n", seq.c_str());
                        continue;
                    }
                    
                    // valid seq
                    printf("%s\n", seq.c_str());
                    
                    // convert to uppercase
                    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

                    // convert T to U
                    replace(seq.begin(), seq.end(), 'T', 'U');
                    
                    seqflag = true;
                    lineIndex ++;
                } else { // constraint
                    constr = input;
                    if (seq.length() != constr.length())
                        printf("The lengths don't match between sequence and constraints: %s, %s\n", seq.c_str(), constr.c_str());
                        // return 0;
                    int n = seq.length();
                    vector<int> cons(n);
                    stack<int> leftBrackets;
                    consflag = true;
                    for (int i=0; i < n; i++){
                        char coni = constr[i];
                        if (consSet.count(coni) == 0){
                            printf("Unrecognized constraint character, should be ? . ( or )\n");
                            consflag = false;
                            break;
                        }
                        switch(coni){
                            case '.':
                                cons[i] = -2;
                                break;
                            case '?':
                                cons[i] = -1;
                                break;
                            case '(':
                                leftBrackets.push(i);
                                break;
                            case ')':
                                int leftIndex = leftBrackets.top();
                                leftBrackets.pop();
                                cons[leftIndex] = i;
                                cons[i] = leftIndex;
                                break;
                        }
                    }
                    
                    seqflag = false;
                    lineIndex ++;

                    if (consflag) {
                        printf("%s\n", constr.c_str());
                        
                        // lhuang: moved inside loop, fixing an obscure but crucial bug in initialization
                        BeamCKYParser parser(beamsize, !sharpturn, is_verbose, is_constraints);

                        BeamCKYParser::DecoderResult result = parser.parse(seq, &cons);

                #ifdef lv
                        double printscore = (result.score / -100.0);
                #else
                        double printscore = result.score;
                #endif
                        printf("%s (%.2f)\n", result.structure.c_str(), printscore);

                        ++num;
                        total_len += seq.length();
                        total_score += result.score;
                        total_states += result.num_states;
                        total_time += result.time;
                    }
                }
            }
        } else {
            for (string seq; getline(cin, seq);) {
                if (seq.length() == 0)
                    continue;

                if (seq[0] == ';' || seq[0] == '>') {
                    printf("%s\n", seq.c_str());
                    continue;
                }

                if (!isalpha(seq[0])){
                    printf("Unrecognized sequence: %s\n", seq.c_str());
                    continue;
                }

                printf("%s\n", seq.c_str());
                
                // convert to uppercase
                transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

                // convert T to U
                replace(seq.begin(), seq.end(), 'T', 'U');

                // lhuang: moved inside loop, fixing an obscure but crucial bug in initialization
                BeamCKYParser parser(beamsize, !sharpturn, is_verbose);

                BeamCKYParser::DecoderResult result = parser.parse(seq, NULL);

        #ifdef lv
                double printscore = (result.score / -100.0);
        #else
                double printscore = result.score;
        #endif
                printf("%s (%.2f)\n", result.structure.c_str(), printscore);

                ++num;
                total_len += seq.length();
                total_score += result.score;
                total_states += result.num_states;
                total_time += result.time;
            }
        }
        // lhuang: TODO add --time switch
        // printf("beam %d\tlen %d\ttime %.5f\tscore %.2f\n", beamsize, seq.length(), result.time, printscore); 

        // ++num;
        // total_len += seq.length();
        // total_score += result.score;
        // total_states += result.num_states;
        // total_time += result.time;
    }

    return 0;
}
