CC=g++
CF=-I../include/ -O3 -std=c++11

run: ./exec_test_sum
	./exec_test_sum 10000000 1000

exec_test_sum: test_sum.cpp
	$(CC) $(CF) $^ -o $@
