/*
 *LinearFold.h*
 header file for LinearFold.cpp.

 author: Kai Zhao, Dezhong Deng, He Zhang
 edited by: 02/2018
*/

#ifndef FASTCKY_BEAMCKYPAR_H
#define FASTCKY_BEAMCKYPAR_H

#include <string>
#include <limits>
#include <vector>
#include <unordered_map>
//#include <stdint.h>

#define MIN_CUBE_PRUNING_SIZE 20

// #define lv

#ifdef lv
  typedef int value_type;
#define VALUE_MIN std::numeric_limits<int>::lowest()
#else
  typedef double value_type;
  #define VALUE_MIN std::numeric_limits<double>::lowest()
#endif

enum Manner {
  MANNER_NONE = 0,              // 0: empty
  MANNER_H,                     // 1: hairpin candidate
  MANNER_HAIRPIN,               // 2: hairpin
  MANNER_SINGLE,                // 3: single
  MANNER_HELIX,                 // 4: helix
  MANNER_MULTI,                 // 5: multi = ..M2. [30 restriction on the left and jump on the right]
  MANNER_MULTI_eq_MULTI_plus_U, // 6: multi = multi + U
  MANNER_P_eq_MULTI,            // 7: P = (multi)
  MANNER_M2_eq_M_plus_P,        // 8: M2 = M + P
  MANNER_M_eq_M2,               // 9: M = M2
  MANNER_M_eq_M_plus_U,         // 10: M = M + U
  MANNER_M_eq_P,                // 11: M = P
  /* MANNER_C_eq_U, */
  /* MANNER_C_eq_P, */
  MANNER_C_eq_C_plus_U,     // 12: C = C + U
  MANNER_C_eq_C_plus_P,     // 13: C = C + P
};

struct State {
    // double score;
    value_type score;
    Manner manner;

    union TraceInfo {
        int split;
        struct {
            char l1;
            int l2;
        } paddings;
    };

    TraceInfo trace;

    State(): manner(MANNER_NONE), score(VALUE_MIN) {};
    State(value_type s, Manner m): score(s), manner(m) {};

    void set(value_type score_, Manner manner_) {
        score = score_; manner = manner_;
    }

    void set(value_type score_, Manner manner_, int split_) {
        score = score_; manner = manner_; trace.split = split_;
    }

    void set(value_type score_, Manner manner_, char l1_, int l2_) {
        score = score_; manner = manner_;
        trace.paddings.l1 = l1_; trace.paddings.l2 = l2_;
    }
};


class BeamCKYParser {
public:
    int beam;

    // int nos_set_update = 0;

    // bool use_vienna;
    // bool is_candidate_list;
    // bool is_cube_pruning;
    bool no_sharp_turn;
    bool is_verbose;
    bool use_constraints; // lisiz, add constraints

    struct DecoderResult {
        std::string structure;
        value_type score;
        unsigned long num_states;
        double time;
    };

    BeamCKYParser(int beam_size=100,
                  // bool vienna=false,
                  // bool candidate_list=true,
                  bool nosharpturn=true,
                  // bool cube_pruning=true,
                  bool is_verbose=false,
                  bool is_constraints=false); // lisiz, add constraints

    // DecoderResult parse(std::string& seq);
    DecoderResult parse(std::string& seq, std::vector<int>* cons);

private:
    void get_parentheses(char* result, std::string& seq);

    unsigned seq_length;

    std::vector<std::unordered_map<int, State>> bestH, bestP, bestM2, bestMulti, bestM;

    std::vector<int> if_tetraloops;
    std::vector<int> if_hexaloops;
    std::vector<int> if_triloops;

    // same as bestM, but ordered
    std::vector<std::vector<std::pair<value_type, int>>> sorted_bestM;

    // hzhang: sort keys in each beam to avoid randomness
    std::vector<std::pair<int, State>> keys;

    // hzhang: sort keys in each beam to avoid randomness
    void sort_keys(std::unordered_map<int, State> &map, std::vector<std::pair<int,State>> &sorted_keys);

    void sortM(value_type threshold,
               std::unordered_map<int, State> &beamstep,
               std::vector<std::pair<value_type, int>>& sorted_stepM);

    std::vector<State> bestC;

    std::vector<int> nucs;

    // void prepare(unsigned len);

    // lisiz: constraints
    std::vector<int> allow_unpaired_position;
    std::vector<int> allow_unpaired_range;
    bool allow_paired(int i, int j, std::vector<int>* cons, char nuci, char nucj);

    void prepare(unsigned len);

    void update_if_better(State &state, value_type newscore, Manner manner) {
      if (state.score < newscore)
            // ++ nos_set_update;
            state.set(newscore, manner);
    };

    void update_if_better(State &state, value_type newscore, Manner manner, int split) {
        if (state.score < newscore || state.manner == MANNER_NONE)
            // ++ nos_set_update;
            state.set(newscore, manner, split);
    };

    void update_if_better(State &state, value_type newscore, Manner manner, char l1, int l2) {
        if (state.score < newscore || state.manner == MANNER_NONE)
            // ++ nos_set_update;
            state.set(newscore, manner, l1, l2);
    };

    value_type beam_prune(std::unordered_map<int, State>& beamstep);

    // vector to store the scores at each beam temporarily for beam pruning
    std::vector<std::pair<value_type, int>> scores;
};

#endif //FASTCKY_BEAMCKYPAR_H
