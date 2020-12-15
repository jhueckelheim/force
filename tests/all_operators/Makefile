CC=g++
CF=-I../../include/ -O3 -std=c++11 -fopenmp

run: ./exec_test_operators
	./exec_test_operators

exec_test_operators: test_operators.cpp
	$(CC) $(CF) $^ -o $@
