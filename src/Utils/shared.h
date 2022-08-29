#ifndef SHARED_H
#define SHARED_H

#define INF 1000000007

#define NOTON 5 // NUM_OF_TYPE_OF_NUCS
#define NOTOND 25
#define NOTONT 125

#define EXPLICIT_MAX_LEN 4
#define SINGLE_MIN_LEN 0
#define SINGLE_MAX_LEN 30  // NOTE: *must* <= sizeof(char), otherwise modify State::TraceInfo accordingly

#define HAIRPIN_MAX_LEN 30
#define BULGE_MAX_LEN SINGLE_MAX_LEN
#define INTERNAL_MAX_LEN SINGLE_MAX_LEN
#define SYMMETRIC_MAX_LEN 15
#define ASYMMETRY_MAX_LEN 28

bool _allowed_pairs[NOTON][NOTON];

// Vienna & CONTRAfold encoding
#define GET_ACGU_NUM_V(x) ((x=='A'? 1 : (x=='C'? 2 : (x=='G'? 3 : (x=='U'?4: 0)))))
#define GET_ACGU_NUM_C(x) ((x=='A'? 0 : (x=='C'? 1 : (x=='G'? 2 : (x=='U'?3: 4)))))

#ifdef lv // lhuang: Vienna: 0N 1A 2C 3G 4U
#define GET_ACGU_NUM(x)   GET_ACGU_NUM_V(x) //((x=='A'? 1 : (x=='C'? 2 : (x=='G'? 3 : (x=='U'?4: 0)))))
#else     // lhuang: CONTRA: 0A 1C 2G 3U 4N
#define GET_ACGU_NUM(x)   GET_ACGU_NUM_C(x) //((x=='A'? 0 : (x=='C'? 1 : (x=='G'? 2 : (x=='U'?3: 4)))))
#endif



void initialize()
{
    _allowed_pairs[GET_ACGU_NUM('A')][GET_ACGU_NUM('U')] = true;
    _allowed_pairs[GET_ACGU_NUM('U')][GET_ACGU_NUM('A')] = true;
    _allowed_pairs[GET_ACGU_NUM('C')][GET_ACGU_NUM('G')] = true;
    _allowed_pairs[GET_ACGU_NUM('G')][GET_ACGU_NUM('C')] = true;
    _allowed_pairs[GET_ACGU_NUM('G')][GET_ACGU_NUM('U')] = true;
    _allowed_pairs[GET_ACGU_NUM('U')][GET_ACGU_NUM('G')] = true;
}

#endif
