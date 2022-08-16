################################
# Makefile
#
# author: He Zhang
# edited by: 11/2018
################################

CC=g++
DEPS=src/LinearFoldEval.h src/LinearFold.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h 
CFLAGS=-std=c++11 -O3
.PHONY : clean linearfold
objects=bin/linearfold_v bin/linearfold_c

linearfold: src/LinearFold.cpp $(DEPS) 
		chmod +x linearfold draw_circular_plot
		mkdir -p bin
		$(CC) src/LinearFold.cpp $(CFLAGS) -Dlv -Dis_cube_pruning -Dis_candidate_list -o bin/linearfold_v 
		$(CC) src/LinearFold.cpp $(CFLAGS) -Dis_cube_pruning -Dis_candidate_list -o bin/linearfold_c

clean:
	-rm $(objects)
