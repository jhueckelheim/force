CC=g++
CF=-I../../include/ -I$(MPFR_HOME)/include -L$(MPFR_HOME)/lib -O3 -std=c++11 -fopenmp -lquadmath -lmpfr -DMPFR

build: ./exec_test_mxm

run: build
	./exec_test_mxm

exec_test_mxm: mxm.cpp
	$(CC) $(CF) $^ -o $@
