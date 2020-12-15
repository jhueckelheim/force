CX=g++
CXF= -O3 -Wall -std=c++11 -fopenmp -I../../include/

all: HACCmk_sgl HACCmk_dble HACCmk_ldbl HACCmk_quad HACCmk_cdbl HACCmk_cldb HACCmk_cena HACCmk_mpfr200

perf: HACCmk_sgl_p HACCmk_dble_p HACCmk_ldbl_p HACCmk_quad_p HACCmk_cdbl_p HACCmk_cldb_p HACCmk_cena_p

test: all
	./HACCmk_sgl
	./HACCmk_dble
	./HACCmk_quad
	./HACCmk_cdbl
	./HACCmk_cena

HACCmk_sgl: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -o $@

HACCmk_dble: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DDBLE -o $@

HACCmk_quad: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -lquadmath -DQUAD -o $@

HACCmk_cena: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DCENA -o $@

HACCmk_cdbl: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DCDBL -o $@

HACCmk_cldb: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DCLDB -o $@

HACCmk_ldbl: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DLDBL -o $@

HACCmk_mpfr200: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=200 -o $@

HACCmk_sgl_p: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DPERF -o $@

HACCmk_dble_p: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DPERF -DDBLE -o $@

HACCmk_quad_p: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DPERF -lquadmath -DQUAD -o $@

HACCmk_cena_p: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DPERF -DCENA -o $@

HACCmk_cdbl_p: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DPERF -DCDBL -o $@

HACCmk_cldb_p: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DPERF -DCLDB -o $@

HACCmk_ldbl_p: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DPERF -DLDBL -o $@

NUMBERS = $(shell seq 1 1 113) 200

allmpfr: main_hacc.c mysecond.c
	$(foreach var,$(NUMBERS),$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=$(var) -o HACCmk_mpfr_$(var);)

allcmpfr: main_hacc.c mysecond.c
	$(foreach var,$(NUMBERS),$(CX) $(CXF) $^ -lmpfr -DCMPFR -DMPFRPR=$(var) -o HACCmk_cmpfr_$(var);)
