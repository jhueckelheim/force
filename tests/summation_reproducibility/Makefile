CC=g++
CF=-I../../include/ -I$(MPFR_HOME)/include -L$(MPFR_HOME)/lib -O3 -std=c++11 -fopenmp -lquadmath -lmpfr -DMPFR

run: ./exec_test_sum_nobranch
	./exec_test_sum_nobranch 1000000 1000

run_classic: ./exec_test_sum_classic
	./exec_test_sum_classic  1000000 1000

exec_test_sum_nobranch: test_sum.cpp
	$(CC) $(CF) $^ -o $@

exec_test_sum_classic: test_sum.cpp
	$(CC) $(CF) -D_CLASSIC_CENA_ADDITION $^ -o $@
