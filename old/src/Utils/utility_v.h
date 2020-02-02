/*
 *utility_v.h*
 provides feature functions for vienna model.

 author: Kai Zhao, Dezhong Deng
 edited by: 02/2018
*/

#ifndef FASTCKY_UTILITY_V_H
#define FASTCKY_UTILITY_V_H

// pairs: 0:NP 1:CG 2:GC 3:GU 4:UG 5:AU 6:UA 7:NN
// nucleotides: CONTRAfold: 0:A 1:C 2:G 3:U 4:N ; Vienna: 0:N 1:A 2:C 3:G 4:U
// TODO: unify
#define NUM_TO_NUC(x) (x==-1?-1:((x==4?0:(x+1))))

// lhuang: NB: before NUM_TO_NUC
#define NUM_TO_PAIR(x,y) (x==0? (y==3?5:0) : (x==1? (y==2?1:0) : (x==2 ? (y==1?2:(y==3?3:0)) : (x==3 ? (y==2?4:(y==0?6:0)) : 0))))
#define NUC_TO_PAIR(x,y) (x==1? (y==4?5:0) : (x==2? (y==3?1:0) : (x==3 ? (y==2?2:(y==4?3:0)) : (x==4 ? (y==3?4:(y==1?6:0)) : 0))))

// bool _allowed_pairs[NOTON][NOTON];

#include <string.h>
#include <cmath>

#include "energy_parameter.h" // energy_parameter stuff
#include "intl11.h"
#include "intl21.h"
#include "intl22.h"

#define MAXLOOP 30

inline int MIN2(int a, int b) {if (a <= b)return a;else return b;}
inline int MAX2(int a, int b) {if (a >= b)return a;else return b;}

/* void v_initialize() */
/* { */
/*     _allowed_pairs[GET_ACGU_NUM('A')][GET_ACGU_NUM('U')] = true; */
/*     _allowed_pairs[GET_ACGU_NUM('U')][GET_ACGU_NUM('A')] = true; */
/*     _allowed_pairs[GET_ACGU_NUM('C')][GET_ACGU_NUM('G')] = true; */
/*     _allowed_pairs[GET_ACGU_NUM('G')][GET_ACGU_NUM('C')] = true; */
/*     _allowed_pairs[GET_ACGU_NUM('G')][GET_ACGU_NUM('U')] = true; */
/*     _allowed_pairs[GET_ACGU_NUM('U')][GET_ACGU_NUM('G')] = true; */

/* } */

inline void v_init_tetra_hex_tri(std::string& seq, int seq_length, std::vector<int>& if_tetraloops, std::vector<int>& if_hexaloops, std::vector<int>& if_triloops) {

    // TetraLoops
    if_tetraloops.resize(seq_length-5<0?0:seq_length-5, -1);
    for (int i = 0; i < seq_length-5; ++i) {
        if (!(seq[i] == 'C' && seq[i+5] == 'G'))
            continue;
        char *ts;
        const char* tl = seq.substr(i,6).c_str();
                if ((ts=strstr(Tetraloops, tl)))
                    if_tetraloops[i] = (ts - Tetraloops)/7;
    }

    // Triloops
    if_triloops.resize(seq_length-4<0?0:seq_length-4, -1);
    for (int i = 0; i < seq_length-4; ++i) {
        if (!((seq[i] == 'C' && seq[i+4] == 'G') || (seq[i] == 'G' && seq[i+4] == 'C')))
            continue;
        char *ts;
        const char* tl = seq.substr(i,5).c_str();
        if ((ts=strstr(Triloops, tl)))
            if_triloops[i] = (ts - Triloops)/6;
    }

    // Hexaloops
    if_hexaloops.resize(seq_length-7<0?0:seq_length-7, -1);
    for (int i = 0; i < seq_length-7; ++i) {
        if (!(seq[i] == 'A' && seq[i+7] == 'U'))
            continue;
        char *ts;
        const char* tl = seq.substr(i,8).c_str();
                if ((ts=strstr(Hexaloops, tl)))
                    if_hexaloops[i] = (ts - Hexaloops)/9;
    }
    return;
}

inline int v_score_hairpin(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int tetra_hex_tri_index = -1) {
    int size = j-i-1;
    int type = NUM_TO_PAIR(nuci, nucj);
    int si1 = NUM_TO_NUC(nuci1);
    int sj1 = NUM_TO_NUC(nucj_1);

    int energy;

    if(size <= 30)
        energy = hairpin37[size];
    else
        energy = hairpin37[30] + (int)(lxc37*log((size)/30.));

    if(size < 3) return energy; /* should only be the case when folding alignments */
#ifdef SPECIAL_HP
    // if(special_hp){
        if (size == 4 && tetra_hex_tri_index > -1)
            return Tetraloop37[tetra_hex_tri_index];
        else if (size == 6 && tetra_hex_tri_index > -1)
            return Hexaloop37[tetra_hex_tri_index];
        else if (size == 3) {
            if (tetra_hex_tri_index > -1)
                return Triloop37[tetra_hex_tri_index];
            return (energy + (type>2 ? TerminalAU37 : 0));
        }
    // }
#endif

    energy += mismatchH37[type][si1][sj1];

    return energy;
}

inline int v_score_single(int i, int j, int p, int q,
                        int nuci, int nuci1, int nucj_1, int nucj,
                        int nucp_1, int nucp, int nucq, int nucq1){
    int si1 = NUM_TO_NUC(nuci1);
    int sj1 = NUM_TO_NUC(nucj_1);
    int sp1 = NUM_TO_NUC(nucp_1);
    int sq1 = NUM_TO_NUC(nucq1);
    int type = NUM_TO_PAIR(nuci, nucj);
    int type_2 = NUM_TO_PAIR(nucq, nucp);
    int n1 = p-i-1;
    int n2 = j-q-1;
    int nl, ns, u, energy;
    energy = 0;

    if (n1>n2) { nl=n1; ns=n2;}
    else {nl=n2; ns=n1;}

    if (nl == 0)
        return stack37[type][type_2];  /* stack */

    if (ns==0) {                      /* bulge */
        energy = (nl<=MAXLOOP)?bulge37[nl]:
      (bulge37[30]+(int)(lxc37*log(nl/30.)));
    if (nl==1) energy += stack37[type][type_2];
    else {
      if (type>2) energy += TerminalAU37;
      if (type_2>2) energy += TerminalAU37;
    }
    return energy;
  }
  else {                            /* interior loop */
    if (ns==1) {
      if (nl==1)                    /* 1x1 loop */
        return int11_37[type][type_2][si1][sj1];
      if (nl==2) {                  /* 2x1 loop */
        if (n1==1)
          energy = int21_37[type][type_2][si1][sq1][sj1];
        else
          energy = int21_37[type_2][type][sq1][si1][sp1];
        return energy;
      }
      else {  /* 1xn loop */
        energy = (nl+1<=MAXLOOP)?(internal_loop37[nl+1]) : (internal_loop37[30]+(int)(lxc37*log((nl+1)/30.)));
        energy += MIN2(MAX_NINIO, (nl-ns)*ninio37);
        energy += mismatch1nI37[type][si1][sj1] + mismatch1nI37[type_2][sq1][sp1];
        return energy;
      }
    }
    else if (ns==2) {
      if(nl==2)      {              /* 2x2 loop */
        return int22_37[type][type_2][si1][sp1][sq1][sj1];}
      else if (nl==3){              /* 2x3 loop */
        energy = internal_loop37[5]+ninio37;
        energy += mismatch23I37[type][si1][sj1] + mismatch23I37[type_2][sq1][sp1];
        return energy;
      }

    }
    { /* generic interior loop (no else here!)*/
      u = nl + ns;
      energy = (u <= MAXLOOP) ? (internal_loop37[u]) : (internal_loop37[30]+(int)(lxc37*log((u)/30.)));

      energy += MIN2(MAX_NINIO, (nl-ns)*ninio37);

      energy += mismatchI37[type][si1][sj1] + mismatchI37[type_2][sq1][sp1];
    }
  }
  return energy;
}

// multi_loop
inline int E_MLstem(int type, int si1, int sj1) {
    int energy = 0;

    if(si1 >= 0 && sj1 >= 0){
        energy += mismatchM37[type][si1][sj1];
    }
    else if (si1 >= 0){
        energy += dangle5_37[type][si1];
    }
    else if (sj1 >= 0){
        energy += dangle3_37[type][sj1];
    }

    if(type > 2) {
        energy += TerminalAU37;
    }

    energy += ML_intern37;

    return energy;
}

inline int v_score_M1(int i, int j, int k, int nuci_1, int nuci, int nuck, int nuck1, int len) {
    int p = i;
    int q = k;
    int tt = NUM_TO_PAIR(nuci, nuck);
    int sp1 = NUM_TO_NUC(nuci_1);
    int sq1 = NUM_TO_NUC(nuck1);

    return E_MLstem(tt, sp1, sq1);

}

inline int v_score_multi_unpaired(int i, int j) {
    return 0;
}

inline int v_score_multi(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int len) {
    int tt = NUM_TO_PAIR(nucj, nuci);
    int si1 = NUM_TO_NUC(nuci1);
    int sj1 = NUM_TO_NUC(nucj_1);

    return E_MLstem(tt, sj1, si1) + ML_closing37;
}

// exterior_loop
inline int v_score_external_paired(int i, int j, int nuci_1, int nuci, int nucj, int nucj1, int len) {
    int type = NUM_TO_PAIR(nuci, nucj);
    int si1 = NUM_TO_NUC(nuci_1);
    int sj1 = NUM_TO_NUC(nucj1);
    int energy = 0;

    if(si1 >= 0 && sj1 >= 0){
        energy += mismatchExt37[type][si1][sj1];
    }
    else if (si1 >= 0){
        energy += dangle5_37[type][si1];
    }
    else if (sj1 >= 0){
        energy += dangle3_37[type][sj1];
    }

    if(type > 2)
        energy += TerminalAU37;
  return energy;
}

inline int v_score_external_unpaired(int i, int j) {
    return 0;
}

#endif //FASTCKY_UTILITY_V_H
