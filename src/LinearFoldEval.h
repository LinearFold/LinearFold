/*
 *LinearFoldEval.cpp*
 Evaluate the energy of a given RNA structure.

 author: He Zhang
 edited by: 12/2018
*/

#ifndef LINEARFOLDEVAL_H
#define LINEARFOLDEVAL_H

#include <stack>
#include <string>
#include <vector>

#include "LinearFold.h"

#include "Utils/utility_v.h"

using namespace std;

long eval(string seq, string ref, bool is_verbose, int dangle_model) {

    int seq_length = seq.length();

    vector<int> if_tetraloops;
    vector<int> if_hexaloops;
    vector<int> if_triloops;

    v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops); // calculate if_tetraloops, if_hexaloops, if_triloops

    vector<int> eval_nucs;
    eval_nucs.clear();
    eval_nucs.resize(seq_length);
    for (int i = 0; i < seq_length; ++i) {
      eval_nucs[i] = GET_ACGU_NUM_V(seq[i]); // lhuang: explicitly use Vienna coding (not very nice)
    }

    long total_energy = 0;
    long external_energy = 0;
    long M1_energy[seq_length];
    long multi_number_unpaired[seq_length];
    // int external_number_unpaired = 0;

    stack<pair<int, int>> stk; // tuple of (index, page)
    tuple<int, int> inner_loop;

    for (int j=0; j<seq_length; j++) {
        M1_energy[j] = 0; // init multi of position j
        multi_number_unpaired[j] = 0;

        if (ref[j] == '.') {
            if (!stk.empty())
                multi_number_unpaired[stk.top().first] += 1;
        }

        else if (ref[j] == '(') {
            if (!stk.empty()) { // +1 for outer loop page
                stk.top().second ++;
            }
            stk.push(make_pair(j, 0)); // init page=0
        }

        else if (ref[j] == ')') {
            assert(!stk.empty());
            tuple<int, int> top = stk.top();
            int i = get<0>(top), page = get<1>(top);
            stk.pop();

            int nuci = eval_nucs[i];
            int nucj = eval_nucs[j];
            int nuci1 = (i + 1) < seq_length ? eval_nucs[i + 1] : -1;
            int nucj_1 = (j - 1) > -1 ? eval_nucs[j - 1] : -1;
            int nuci_1 = (i-1>-1) ? eval_nucs[i-1] : -1; // only for calculating v_score_M1
            int nucj1 = (j+1) < seq_length ? eval_nucs[j+1] : -1; // only for calculating v_score_M1

            if (page == 0) { // hairpin
                int tetra_hex_tri = -1;
                if (j-i-1 == 4) // 6:tetra
                    tetra_hex_tri = if_tetraloops[i];
                else if (j-i-1 == 6) // 8:hexa
                    tetra_hex_tri = if_hexaloops[i];
                else if (j-i-1 == 3) // 5:tri
                    tetra_hex_tri = if_triloops[i];
                
                int newscore = - v_score_hairpin(i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri);
                if (is_verbose)
                    printf("Hairpin loop ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], newscore / -100.0);
                total_energy += newscore;
            }

            else if (page == 1) { //single
                int p = get<0>(inner_loop), q = get<1>(inner_loop);

                int nucp_1 = eval_nucs[p-1], nucp = eval_nucs[p], nucq = eval_nucs[q], nucq1 = eval_nucs[q+1];

                int newscore = - v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj,
                                                  nucp_1, nucp, nucq, nucq1);
                if (is_verbose)
                    printf("Interior loop ( %d, %d) %c%c; ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], p+1, q+1, seq[p],seq[q], newscore / -100.0);
                total_energy += newscore;
            }

            else { //multi
                int multi_score = 0;
                multi_score += M1_energy[i];
                multi_score += - v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, seq_length, dangle_model);
                multi_score += - v_score_multi_unpaired(i+1, i + multi_number_unpaired[i]); // current model is 0
                if (is_verbose)
                    printf("Multi loop ( %d, %d) %c%c : %.2f\n", i+1, j+1, seq[i], seq[j], multi_score / -100.0);
                total_energy += multi_score;
            }

            //update inner_loop
            inner_loop = make_tuple(i, j);

            // possible M
            if (!stk.empty())
                M1_energy[stk.top().first] += - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length, dangle_model);

            // check if adding external energy
            if (stk.empty()) {
                int k = i - 1;
                int nuck = k > -1 ? eval_nucs[k] : -1;
                int nuck1 = eval_nucs[k+1];
                external_energy +=  - v_score_external_paired(k+1, j, nuck, nuck1,
                                                            nucj, nucj1, seq_length, dangle_model);
                // external_energy += 0; currently external unpaired is 0
            }
        }
    }

    if (is_verbose)
        printf("External loop : %.2f\n", external_energy / -100.0);
    total_energy += external_energy;
    return total_energy;
}

#endif // LINEARFOLDEVAL_H
