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

#include "LinearFold.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"
#include "LinearFoldEval.cpp" // adding eval mode

#define SPECIAL_HP

using namespace std;

#ifdef lv
    bool comparefunc(pair<int,State> a, pair<int,State> b) {
        return a.first > b.first;
    }

    void BeamCKYParser::sort_keys(unordered_map<int, State> &map, vector<pair<int,State>> &sorted_keys) {
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
                }
                break;
            case MANNER_SINGLE:
                {
                    result[i] = '(';
                    result[j] = ')';
                    p = i + state.trace.paddings.l1;
                    q = j - state.trace.paddings.l2;
                    stk.push(make_tuple(p, q, bestP[q][p]));
                }
                break;
            case MANNER_HELIX:
                {
                    result[i] = '(';
                    result[j] = ')';
                    stk.push(make_tuple(i+1, j-1, bestP[j-1][i+1]));
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
                break;
            case MANNER_M2_eq_M_plus_P:
                k = state.trace.split;
                stk.push(make_tuple(i, k, bestM[k][i]));
                stk.push(make_tuple(k+1, j, bestP[j][k+1]));
                break;
            case MANNER_M_eq_M2:
                stk.push(make_tuple(i, j, bestM2[j][i]));
                break;
            case MANNER_M_eq_M_plus_U:
                stk.push(make_tuple(i, j-1, bestM[j-1][i]));
                break;
            case MANNER_M_eq_P:
                stk.push(make_tuple(i, j, bestP[j][i]));
                break;
            case MANNER_C_eq_C_plus_U:
                k = j - 1;
                if (k != -1)
                    stk.push(make_tuple(0, k, bestC[k]));
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
                }
                break;
            default:  // MANNER_NONE or other cases
                printf("wrong manner at %d, %d: manner %d\n", i, j, state.manner); fflush(stdout);
                assert(false);
        }
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
        value_type newscore = (k >= 0 ? bestC[k].score : 0) + cand.score;
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
  //sorted_stepM.clear();
    if (threshold == VALUE_MIN) {
        // no beam pruning before, so scores vector not usable
        for (auto &item : beamstep) {
            int i = item.first;
            State &cand = item.second;
            int k = i - 1;
            value_type newscore = (k >= 0 ? bestC[k].score : 0) + cand.score;
            sorted_stepM.push_back(make_pair(newscore, i));
        }
    } else {
        for (auto &p : scores) {
            if (p.first >= threshold) sorted_stepM.push_back(p);
        }
    }

    sort(sorted_stepM.begin(), sorted_stepM.end(), std::greater<pair<value_type, int>>());
}


void BeamCKYParser::prepare(unsigned len) {
    seq_length = len;
    nucs = new int[seq_length];
    bestC = new State[seq_length];
    bestH = new unordered_map<int, State>[seq_length];
    bestP = new unordered_map<int, State>[seq_length];
    bestM = new unordered_map<int, State>[seq_length];
    bestM2 = new unordered_map<int, State>[seq_length];
    bestMulti = new unordered_map<int, State>[seq_length];
     
#ifdef is_cube_pruning
    sorted_bestM = new vector<std::pair<value_type, int>>[seq_length];
#endif

    scores.reserve(seq_length);
}

void BeamCKYParser::postprocess() {

    delete[] bestC;  
    delete[] bestH;  
    delete[] bestP;  
    delete[] bestM;  
    delete[] bestM2;  
    delete[] bestMulti;  

    delete[] nucs;  
    delete[] sorted_bestM;
}



BeamCKYParser::DecoderResult BeamCKYParser::parse(string& seq) {

    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

    prepare(static_cast<unsigned>(seq.length()));

    for (int i = 0; i < seq_length; ++i)
        nucs[i] = GET_ACGU_NUM(seq[i]);

    int *next_pair = new int[NOTON * seq_length]{-1};
    {
      for (int nuci = 0; nuci < NOTON; ++nuci) {
	int next = -1;
	for (int j = seq_length-1; j >=0; --j) {
	  next_pair[nuci * seq_length + j] = next;
	  if (_allowed_pairs[nuci][nucs[j]]) next = j;
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
        if(seq_length > 0) bestC[0].set(0, MANNER_C_eq_C_plus_U);
        if(seq_length > 1) bestC[1].set(0, MANNER_C_eq_C_plus_U);
#else
        if(seq_length > 0) bestC[0].set(score_external_unpaired(0, 0), MANNER_C_eq_C_plus_U);
        if(seq_length > 1) bestC[1].set(score_external_unpaired(0, 1), MANNER_C_eq_C_plus_U);
#endif

    // from left to right
	value_type newscore;
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
                int jnext = next_pair[nucj * seq_length + j];
                if (no_sharp_turn) while (jnext - j < 4 && jnext != -1) jnext = next_pair[nucj * seq_length + jnext];
                if (jnext != -1) {
                    int nucjnext = nucs[jnext];
                    int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                    
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
                    // ++ nos_H;
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
                    State &state = item.second;

                    int nuci = nucs[i];
                    int jnext = next_pair[nuci*seq_length+j];

                    if (jnext != -1) {
                        int nuci1 = (i + 1) < seq_length ? nucs[i + 1] : -1;
                        int nucjnext = nucs[jnext];
                        int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                        // 1. extend h(i, j) to h(i, jnext)
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
                        // ++nos_H;
                    }

                    // 2. generate p(i, j)
                    {
                        update_if_better(beamstepP[i], state.score, MANNER_HAIRPIN);
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
                int jnext = next_pair[nuci*seq_length+j];

                // 1. extend (i, j) to (i, jnext)
                {
                    char new_l1 = state.trace.paddings.l1;
                    int new_l2 = state.trace.paddings.l2 + jnext - j;
                    if (jnext != -1) {
                        // 1. extend (i, j) to (i, jnext)

#ifdef lv
                            newscore = state.score;
#else
                            newscore = state.score + score_multi_unpaired(j, jnext - 1);
#endif
                        // this candidate must be the best one at [i, jnext]
                        // so no need to check the score
                        update_if_better(bestMulti[jnext][i], newscore, MANNER_MULTI_eq_MULTI_plus_U,
                                         new_l1,
                                         new_l2
                        );
                    }
                }

                // 2. generate P (i, j)
                {

#ifdef lv
                        newscore = state.score -
                            v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
#else
                        newscore = state.score +
                            score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
#endif
                    update_if_better(beamstepP[i], newscore, MANNER_P_eq_MULTI);
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

                // 1. generate new helix / single_branch
                // new state is of shape p..i..j..q
                if (i >0 && j<seq_length-1) {
                    value_type precomputed;
#ifndef lv
                    precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci);
#endif
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; // hzhang: move here
                        int q = next_pair[nucp*seq_length+j];
                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];
                            int nucq_1 = nucs[q - 1];

#ifdef lv
			     newscore = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
							 nuci_1, nuci, nucj, nucj1)
				+ state.score;

			     update_if_better(bestP[q][p], newscore, MANNER_SINGLE, static_cast<char>(i - p),
					      q - j);
      
#else
                            if (p == i - 1 && q == j + 1) {
                                    newscore = score_helix(nucp, nucp1, nucq_1, nucq) + state.score;
                                     update_if_better(bestP[q][p], newscore, MANNER_HELIX);
                            } else {
                                    newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
                                        precomputed +
                                        score_single_without_junctionB(p, q, i, j,
                                                                       nuci_1, nuci, nucj, nucj1) +
                                        state.score;
                                update_if_better(bestP[q][p], newscore, MANNER_SINGLE,
                                                 static_cast<char>(i - p),
                                                 q - j);
			    }
#endif
                            

                            q = next_pair[nucp*seq_length+q];
                        }
                    }
                }
                //printf(" helix / single at %d\n", j); fflush(stdout);

                // 2. M = P
                if(i > 0 && j < seq_length-1){

#ifdef lv
                        newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#else
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score;
#endif
                    update_if_better(beamstepM[i], newscore, MANNER_M_eq_P);
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
                                newscore = M1_score + m.second.score;
                                update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                            }
#else
                        if (bestM2_iter==beamstepM2.end() || M1_score > bestM2_iter->second.score) {
                            for (auto &m : bestM[k]) {
                                int newi = m.first;
                                // eq. to first convert P to M1, then M2/M = M + M1
                                newscore = M1_score + m.second.score;
                                update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
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
                      }
                    } else {
                       
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
                    }
                }
                //printf(" C = C + P at %d\n", j); fflush(stdout);
            }

            if (use_cube_pruning) {
                // 3. M2 = M + P with cube pruning
                vector<int> valid_Ps;
                vector<value_type> M1_scores;
#ifdef lv
		//sort_keys(beamstepP, keys);
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
		      			//assert(bestM[k].size() == sorted_bestM[k].size());
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
                     newscore = M1_scores[index_P] + bestM[k][newi].score;
                    pop_heap(heap.begin(), heap.end());
                    heap.pop_back();

                    if (beamstepM2[newi].manner == MANNER_NONE) {
                        ++filled;
                        update_if_better(beamstepM2[newi], newscore, MANNER_M2_eq_M_plus_P, k);
                        // ++nos_M2;
                    } 
		    		//else assert(beamstepM2[newi].score > newscore - 1e-8);
                    

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
                            //assert(beamstepM2[candidate_newi].score >
                            //       M1_scores[index_P] + bestM[k][candidate_newi].score - 1e-8);
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

                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp*seq_length+j];
                        if (q != -1 && (i - p - 1 <= SINGLE_MAX_LEN)) {
                            // the current shape is p..i M2 j ..q

#ifdef lv
                                newscore = state.score;
#else
                                newscore = score_multi_unpaired(p+1, i-1) +
                                    score_multi_unpaired(j+1, q-1) + state.score;
#endif
                            update_if_better(bestMulti[q][p], newscore, MANNER_MULTI,
                                             static_cast<char>(i - p),
                                             q - j);
                        }
                    }
                }

                // 2. M = M2
                {
                    update_if_better(beamstepM[i], state.score, MANNER_M_eq_M2);
                    // ++ nos_M;
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
                  // if (use_vienna)
#ifdef lv
                        newscore = state.score;
                    // else
#else
                        newscore = score_multi_unpaired(j + 1, j + 1) + state.score;
#endif
                    update_if_better(bestM[j+1][i], newscore, MANNER_M_eq_M_plus_U);
                    // ++ nos_M;
                }
            }
        }
        //printf("M at %d\n", j); fflush(stdout);

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {

#ifdef lv
                    newscore = beamstepC.score;
#else
                    newscore = score_external_unpaired(j+1, j+1) + beamstepC.score;
#endif
                update_if_better(bestC[j+1], newscore, MANNER_C_eq_C_plus_U);
                // ++ nos_C;
            }
        }
        //printf("C at %d\n", j); fflush(stdout);

    }  // end of for-loo j

    State& viterbi = bestC[seq_length-1];
    value_type total_score = viterbi.score; 
    char result[seq_length+1];
    get_parentheses(result, seq);

	postprocess();

    // return {string(result), viterbi.score, nos_tot, parse_elapsed_time};
    return {string(result), total_score};
}

BeamCKYParser::BeamCKYParser(int beam_size,
                             bool nosharpturn,
                             bool verbose)
    : beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      is_verbose(verbose) {
#ifdef lv
        initialize();
#else
        initialize();
        initialize_cachesingle();
#endif
}


// -------------------------------------------------------------

int main(int argc, char** argv){

        struct timeval total_starttime, total_endtime;
    gettimeofday(&total_starttime, NULL);

    int beamsize = 100;
    bool sharpturn = false;
    bool is_verbose = false;
    bool is_eval = false;

    if (argc > 1) {
        beamsize = atoi(argv[1]);
        sharpturn = atoi(argv[2]) == 1;
        is_verbose = atoi(argv[3]) == 1;
        is_eval = atoi(argv[4]) == 1; // adding eval mode
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
                if (ref.length() == 0) {
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

                double MFE_energy = eval(seq, ref, is_verbose) / -100.0;
                // printf("structure energy: %.2f\n", MFE_energy/ -100.0);
                printf("%s\n", seq.c_str());
                printf("%s (%.2f)\n", input.c_str(), MFE_energy);

            }
            lineIndex ++;
        }
    }

    else {
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

	    BeamCKYParser::DecoderResult result = parser.parse(seq);

    #ifdef lv
             double printscore = (result.score / -100.0);
    #else
             double printscore = result.score;
#endif
	     gettimeofday(&total_endtime, NULL);
	     double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;
	     printf("%s (%.2f)\n", result.structure.c_str(), printscore);
        }
    }

    return 0;
}
