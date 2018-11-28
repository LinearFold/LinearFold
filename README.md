# LinearFold: Linear-Time Prediction for RNA Secondary Structures
This repository contains the C++ source code for the LinearFold project, the first linear-time prediction algorithm/software for RNA secondary structures.

Preprint: LinearFold: Linear-Time Prediction of RNA Secondary Structures

Liang Huang*, He Zhang, Dezhong Deng, Kai Zhao, Kaibo Liu, David Hendrix, David Mathews

*corresponding author

Web server: http://linearfold.eecs.oregonstate.edu

## Dependencies
g++4.7 or above
python2.7


## To Compile and Make script Runnable
```
make
chmod +x linearfold
```

## To Run
The LinearFold parser can be run with:
```
echo "SEQUENCE" | linearfold [OPTIONS]

OR

cat SEQ_OR_FASTA_FILE | linearfold [OPTIONS]
```
Both FASTA format and pure-sequence format are supported for input.

OPTIONS:

-b

The beam size (default 100). Use 0 for infinite beam.

-V

Switches LinearFold-C (by default) to LinearFold-V.

--verbose

Prints out energy of each loop in the structure. (default False)

--sharpturn

Enable sharpturn in prediction。 (default False)

## Example Run

cat testseq | ./linearfold 

echo "GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA" | ./linearfold -b 20 -V --verbose --sharpturn
