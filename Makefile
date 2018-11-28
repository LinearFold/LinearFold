################################
# Makefile
#
# author: He Zhang
# edited by: 11/2018
################################

CC=g++
DEPS=LinearFold.h energy_parameter.h feature_weight.h intl11.h intl21.h intl22.h utility_v.h utility.h 
CFLAGS=-std=c++11 -O3
.PHONY : clean linearfold
objects=linearfold_v linearfold_c

linearfold: LinearFold.cpp $(DEPS) 
		chmod +x linearfold
		$(CC) LinearFold.cpp $(CFLAGS) -Dlv -Dis_cube_pruning -Dis_candidate_list -o linearfold_v 
		$(CC) LinearFold.cpp $(CFLAGS) -Dis_cube_pruning -Dis_candidate_list -o linearfold_c

clean:
	-rm $(objects)