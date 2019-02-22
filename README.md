This codebase replaces the now deprecated version: https://github.com/abentu0101/LinearFold.
This version fixes many bugs and design problems in the old version.

# LinearFold: Linear-Time Prediction for RNA Secondary Structures

This repository contains the C++ source code for the LinearFold project, the first linear-time prediction algorithm/software for RNA secondary structures.

Preprint: LinearFold: Linear-Time Prediction of RNA Secondary Structures

Liang Huang*, He Zhang**, Dezhong Deng**, Kai Zhao, Kaibo Liu, David Hendrix, David Mathews

\* corresponding author

** contributed equally

Web server: http://linearfold.eecs.oregonstate.edu


## Dependencies
gcc 4.8.5 or above; 
python2.7

## To Compile
```
make
```

## To Run
The LinearFold parser can be run with:
```
echo SEQUENCE | ./linearfold [OPTIONS]

OR

cat SEQ_OR_FASTA_FILE | ./linearfold [OPTIONS]
```
Both FASTA format and pure-sequence format are supported for input.

OPTIONS:
```
-b BEAM_SIZE
```
The beam size (default 100). Use 0 for infinite beam.
```
-V
```
Switches LinearFold-C (by default) to LinearFold-V.
```
--verbose
```
Prints out energy of each loop in the structure. (default False)
```
--sharpturn
```
Enable sharpturn in prediction. (default False)
```
--eval
```
Enable eval mode, which can calculate free energy for a given structure of a sequence. (default False)

## Example Run Predict
```
cat testseq | ./linearfold
UGAGUUCUCGAUCUCUAAAAUCG
....................... (-0.22)
AAAACGGUCCUUAUCAGGACCAAACA
.....((((((....))))))..... (4.91)
AUUCUUGCUUCAACAGUGUUUGAACGGAAU
.............................. (-0.29)
UCGGCCACAAACACACAAUCUACUGUUGGUCGA
(((((((...................))))))) (0.99)
GUUUUUAUCUUACACACGCUUGUGUAAGAUAGUUA
.....(((((((((((....))))))))))).... (6.66)

echo GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA | ./linearfold -V -b 20
GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA
(((((((..((((.......))))((((((((...)))))))).(((((.......)))))))))))).... (-31.50)
```

## Example Run Eval
```
cat testeval | ./linearfold --eval
UGAGUUCUCGAUCUCUAAAAUCG
.(((........)))........ (-1.80)
AAAACGGUCCUUAUCAGGACCAAACA
.....((((((....))))))..... (-9.30)
AUUCUUGCUUCAACAGUGUUUGAACGGAAU
(((((...(((((......))))).))))) (-6.80)
UCGGCCACAAACACACAAUCUACUGUUGGUCGA
(((((((((..............))).)))))) (-7.80)
GUUUUUAUCUUACACACGCUUGUGUAAGAUAGUUA
....((((((((((((....))))))))))))... (-13.00)

echo -e "GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA\n(((((((..((((.......))))((((((((...)))))))).(((((.......))))))))))))....\n" | ./linearfold --eval --verbose
Hairpin loop ( 13, 21) CG : 4.50
Interior loop ( 12, 22) UA; ( 13, 21) CG : -2.40
Interior loop ( 11, 23) AU; ( 12, 22) UA : -1.10
Interior loop ( 10, 24) GC; ( 11, 23) AU : -2.40
Hairpin loop ( 32, 36) UA : 5.90
Interior loop ( 31, 37) UA; ( 32, 36) UA : -0.90
Interior loop ( 30, 38) CG; ( 31, 37) UA : -2.10
Interior loop ( 29, 39) CG; ( 30, 38) CG : -3.30
Interior loop ( 28, 40) UA; ( 29, 39) CG : -2.40
Interior loop ( 27, 41) UA; ( 28, 40) UA : -0.90
Interior loop ( 26, 42) CG; ( 27, 41) UA : -2.10
Interior loop ( 25, 43) GC; ( 26, 42) CG : -3.40
Hairpin loop ( 49, 57) GC : 4.40
Interior loop ( 48, 58) GC; ( 49, 57) GC : -3.30
Interior loop ( 47, 59) GC; ( 48, 58) GC : -3.30
Interior loop ( 46, 60) UA; ( 47, 59) GC : -2.10
Interior loop ( 45, 61) CG; ( 46, 60) UA : -2.10
Multi loop ( 7, 62) GC : 1.40
Interior loop ( 6, 63) CG; ( 7, 62) GC : -2.40
Interior loop ( 5, 64) UA; ( 6, 63) CG : -2.40
Interior loop ( 4, 65) CG; ( 5, 64) UA : -2.10
Interior loop ( 3, 66) GU; ( 4, 65) CG : -2.50
Interior loop ( 2, 67) GC; ( 3, 66) GU : -1.50
Interior loop ( 1, 68) GC; ( 2, 67) GC : -3.30
External loop : -1.70
GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA
(((((((..((((.......))))((((((((...)))))))).(((((.......)))))))))))).... (-31.50)
```
