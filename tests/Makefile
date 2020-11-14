CC=g++
CF=-I../include/ -O3 -std=c++11

run: ./exec_test_sum_classic ./exec_test_sum_nobranch
	./exec_test_sum_classic 10000000 1000
	./exec_test_sum_nobranch 10000000 1000

exec_test_sum_nobranch: test_sum.cpp
	$(CC) $(CF) $^ -o $@

exec_test_sum_classic: test_sum.cpp
	$(CC) $(CF) -D_CLASSIC_CENA_ADDITION $^ -o $@
