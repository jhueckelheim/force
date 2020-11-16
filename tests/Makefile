CC=g++
CF=-I../include/ -O3 -std=c++11

run: ./exec_test_sum_classic ./exec_test_sum_nobranch ./exec_test_operators ./exec_test_mxm_single ./exec_test_mxm_double ./exec_test_mxm_freal
	./exec_test_sum_classic 10000000 1000
	./exec_test_sum_nobranch 10000000 1000
	./exec_test_operators
	./exec_test_mxm_single
	./exec_test_mxm_double
	./exec_test_mxm_freal

exec_test_sum_nobranch: test_sum.cpp
	$(CC) $(CF) $^ -o $@

exec_test_sum_classic: test_sum.cpp
	$(CC) $(CF) -D_CLASSIC_CENA_ADDITION $^ -o $@

exec_test_operators: test_operators.cpp
	$(CC) $(CF) $^ -o $@

exec_test_mxm_single: mxm.cpp
	$(CC) $(CF) $^ -DLEAF_SIZE=4 -Dreal1=float -Dreal2=double -DN=128 -o $@

exec_test_mxm_double: mxm.cpp
	$(CC) $(CF) $^ -DLEAF_SIZE=4 -Dreal1=double -Dreal2=double -DN=128 -o $@

exec_test_mxm_freal: mxm.cpp
	$(CC) $(CF) $^ -DLEAF_SIZE=4 -Dreal1=float -Dreal2="freal<float>" -DN=128 -o $@

