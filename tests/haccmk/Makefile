CX=g++
CXF= -O3 -Wall -std=c++11 -fopenmp -I../../include/ -lmpfr -lquadmath

all: HACCmk_sgl HACCmk_dble HACCmk_quad HACCmk_cdbl

HACCmk_sgl: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -o $@

HACCmk_dble: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -DDBLE -o $@

HACCmk_quad: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -DQUAD -o $@

HACCmk_cena: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -DCENA -o $@

HACCmk_cdbl: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -DCDBL -o $@
