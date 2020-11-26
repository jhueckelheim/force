CX=g++
CXF= -O3 -Wall -std=c++11 -fopenmp -I../../include/

all: HACCmk_sgl HACCmk_dble HACCmk_ldbl HACCmk_quad HACCmk_cdbl HACCmk_cldb HACCmk_cena HACCmk_mpfr200

test: all
	./HACCmk_sgl
	./HACCmk_dble
	./HACCmk_quad
	./HACCmk_cdbl
	./HACCmk_cena

HACCmk_sgl: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -o $@

HACCmk_dble: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -DDBLE -o $@

HACCmk_quad: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lquadmath -DQUAD -o $@

HACCmk_cena: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -DCENA -o $@

HACCmk_cdbl: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -DCDBL -o $@

HACCmk_cldb: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -DCLDB -o $@

HACCmk_ldbl: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -DLDBL -o $@

HACCmk_mpfr200: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=200 -o $@

allmpfr: main_hacc.c mysecond.c
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=1 -o HACCmk_mpfr_1
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=2 -o HACCmk_mpfr_2
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=3 -o HACCmk_mpfr_3
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=4 -o HACCmk_mpfr_4
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=5 -o HACCmk_mpfr_5
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=6 -o HACCmk_mpfr_6
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=7 -o HACCmk_mpfr_7
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=8 -o HACCmk_mpfr_8
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=9 -o HACCmk_mpfr_9
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=10 -o HACCmk_mpfr_10
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=11 -o HACCmk_mpfr_11
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=12 -o HACCmk_mpfr_12
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=13 -o HACCmk_mpfr_13
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=14 -o HACCmk_mpfr_14
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=15 -o HACCmk_mpfr_15
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=16 -o HACCmk_mpfr_16
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=17 -o HACCmk_mpfr_17
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=18 -o HACCmk_mpfr_18
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=19 -o HACCmk_mpfr_19
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=20 -o HACCmk_mpfr_20
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=21 -o HACCmk_mpfr_21
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=22 -o HACCmk_mpfr_22
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=23 -o HACCmk_mpfr_23
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=24 -o HACCmk_mpfr_24
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=25 -o HACCmk_mpfr_25
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=26 -o HACCmk_mpfr_26
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=27 -o HACCmk_mpfr_27
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=28 -o HACCmk_mpfr_28
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=29 -o HACCmk_mpfr_29
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=30 -o HACCmk_mpfr_30
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=31 -o HACCmk_mpfr_31
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=32 -o HACCmk_mpfr_32
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=33 -o HACCmk_mpfr_33
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=34 -o HACCmk_mpfr_34
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=35 -o HACCmk_mpfr_35
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=36 -o HACCmk_mpfr_36
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=37 -o HACCmk_mpfr_37
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=38 -o HACCmk_mpfr_38
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=39 -o HACCmk_mpfr_39
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=40 -o HACCmk_mpfr_40
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=41 -o HACCmk_mpfr_41
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=42 -o HACCmk_mpfr_42
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=43 -o HACCmk_mpfr_43
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=44 -o HACCmk_mpfr_44
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=45 -o HACCmk_mpfr_45
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=46 -o HACCmk_mpfr_46
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=47 -o HACCmk_mpfr_47
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=48 -o HACCmk_mpfr_48
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=49 -o HACCmk_mpfr_49
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=50 -o HACCmk_mpfr_50
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=51 -o HACCmk_mpfr_51
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=52 -o HACCmk_mpfr_52
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=53 -o HACCmk_mpfr_53
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=54 -o HACCmk_mpfr_54
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=55 -o HACCmk_mpfr_55
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=56 -o HACCmk_mpfr_56
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=57 -o HACCmk_mpfr_57
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=58 -o HACCmk_mpfr_58
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=59 -o HACCmk_mpfr_59
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=60 -o HACCmk_mpfr_60
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=61 -o HACCmk_mpfr_61
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=62 -o HACCmk_mpfr_62
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=63 -o HACCmk_mpfr_63
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=64 -o HACCmk_mpfr_64
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=65 -o HACCmk_mpfr_65
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=66 -o HACCmk_mpfr_66
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=67 -o HACCmk_mpfr_67
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=68 -o HACCmk_mpfr_68
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=69 -o HACCmk_mpfr_69
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=70 -o HACCmk_mpfr_70
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=71 -o HACCmk_mpfr_71
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=72 -o HACCmk_mpfr_72
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=73 -o HACCmk_mpfr_73
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=74 -o HACCmk_mpfr_74
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=75 -o HACCmk_mpfr_75
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=76 -o HACCmk_mpfr_76
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=77 -o HACCmk_mpfr_77
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=78 -o HACCmk_mpfr_78
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=79 -o HACCmk_mpfr_79
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=80 -o HACCmk_mpfr_80
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=81 -o HACCmk_mpfr_81
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=82 -o HACCmk_mpfr_82
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=83 -o HACCmk_mpfr_83
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=84 -o HACCmk_mpfr_84
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=85 -o HACCmk_mpfr_85
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=86 -o HACCmk_mpfr_86
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=87 -o HACCmk_mpfr_87
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=88 -o HACCmk_mpfr_88
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=89 -o HACCmk_mpfr_89
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=90 -o HACCmk_mpfr_90
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=91 -o HACCmk_mpfr_91
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=92 -o HACCmk_mpfr_92
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=93 -o HACCmk_mpfr_93
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=94 -o HACCmk_mpfr_94
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=95 -o HACCmk_mpfr_95
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=96 -o HACCmk_mpfr_96
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=97 -o HACCmk_mpfr_97
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=98 -o HACCmk_mpfr_98
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=99 -o HACCmk_mpfr_99
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=100 -o HACCmk_mpfr_100
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=101 -o HACCmk_mpfr_101
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=102 -o HACCmk_mpfr_102
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=103 -o HACCmk_mpfr_103
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=104 -o HACCmk_mpfr_104
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=105 -o HACCmk_mpfr_105
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=106 -o HACCmk_mpfr_106
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=107 -o HACCmk_mpfr_107
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=108 -o HACCmk_mpfr_108
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=109 -o HACCmk_mpfr_109
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=110 -o HACCmk_mpfr_110
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=111 -o HACCmk_mpfr_111
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=112 -o HACCmk_mpfr_112
	$(CX) $(CXF) $^ -lmpfr -DMPFR -DMPFRPR=113 -o HACCmk_mpfr_113
