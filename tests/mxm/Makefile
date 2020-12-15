CC=g++
CF=-I../include/ -O3 -std=c++11 -fopenmp

run: ./exec_test_mxm_single ./exec_test_mxm_double ./exec_test_mxm_freal
	./exec_test_mxm_single
	./exec_test_mxm_double
	./exec_test_mxm_freal

exec_test_operators: test_operators.cpp
	$(CC) $(CF) $^ -o $@

exec_test_mxm_single: mxm.cpp
	$(CC) $(CF) $^ -DLEAF_SIZE=4 -Dreal1=double -Dreal2=float -DN=128 -o $@

exec_test_mxm_double: mxm.cpp
	$(CC) $(CF) $^ -DLEAF_SIZE=4 -Dreal1=double -Dreal2=double -DN=128 -o $@

exec_test_mxm_freal: mxm.cpp
	$(CC) $(CF) $^ -DLEAF_SIZE=4 -Dreal1=double -Dreal2="freal<float>" -DN=128 -o $@

