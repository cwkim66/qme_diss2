#DEFS = -DMPI

ifeq ($(DEFS), -DMPI)
FC = mpifort
CC = mpicc
CPP = mpicxx
else
FC = gfortran
CC = gcc
CPP = g++
endif

FOPT = -g -w -O2 $(DEFS)
COPT = -g -w -O2 $(DEFS)

MATH      = $(PWD)/math
SPD       = $(PWD)/spd
PTB2      = $(PWD)/ptb2
TRACE     = $(PWD)/trace
DYN       = $(PWD)/dynamics

INTEL     = /appl/intel/oneapi

INTELI    = -I $(INTEL)/include
INTELL    = -L $(INTEL)/lib/intel64

MYI       = $(PWD)/include
MYI_MATH  = $(MATH)/include
MYI_SPD   = $(SPD)/include
MYI_PTB2  = $(PTB2)/include
MYI_TRACE = $(TRACE)/include
MYI_DYN   = $(DYN)/include
MKLI      = $(INTEL)/mkl/latest/include
MYII      = -I $(MYI) -I $(MYI_MATH) -I $(MYI_SPD) -I $(MYI_PTB2) -I $(MYI_TRACE) -I $(MYI_DYN) -I $(MKLI)

MKLL      = $(INTEL)/mkl/latest/lib/intel64
OMPL      = $(INTEL)/lib/intel64
MYLL      = -L $(MKLL) -L /usr/lib/gcc/x86_64-linux-gnu/9

MKLFLAGS  = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
FLAGS     = -lgfortran $(MKLFLAGS)

all: qme.x

clean:
	rm -f *.o $(MATH)/*.o $(SPD)/*.o $(PTB2)/*.o $(TRACE)/*.o $(DYN)/*.o *.x

qme.x: qme.o set_param.o util.o\
	 $(MATH)/mathlib_all.o $(SPD)/spd_all.o $(PTB2)/ptb2_all.o $(TRACE)/trace_all.o $(DYN)/dyn_all.o
	$(CPP) -O2 -o qme.x\
	 qme.o set_param.o util.o\
	 $(MATH)/mathlib_all.o $(SPD)/spd_all.o $(PTB2)/ptb2_all.o $(TRACE)/trace_all.o $(DYN)/dyn_all.o $(MYLL) $(FLAGS)\
	 -o qme.x

qme.o: qme.c $(MYI)/qme.h
	$(CC) $(COPT) $(MYII) -c qme.c

set_param.o: set_param.c $(MYI)/set_param.h
	$(CC) $(COPT) $(MYII) -c set_param.c

util.o: util.c $(MYI)/util.h
	$(CC) $(COPT) $(MYII) -c util.c


### mathlib ###
$(MATH)/mathlib_all.o: $(MATH)/mathlib.o $(MATH)/rand.o
	ld -r -o $(MATH)/mathlib_all.o $(MATH)/mathlib.o $(MATH)/rand.o

$(MATH)/mathlib.o: $(MATH)/mathlib.c $(MYI_MATH)/mathlib.h
	$(CC) $(COPT) $(MYII) $(MKLFLAGS) -o $(MATH)/mathlib.o -c $(MATH)/mathlib.c

$(MATH)/rand.o: $(MATH)/rand.c $(MYI_MATH)/rand.h
	$(CC) $(COPT) $(MYII) $(MKLFLAGS) -o $(MATH)/rand.o -c $(MATH)/rand.c


### spd ###
$(SPD)/spd_all.o: $(SPD)/spd.o
	ld -r -o $(SPD)/spd_all.o $(SPD)/spd.o

$(SPD)/spd.o: $(SPD)/spd.c $(MYI_SPD)/spd.h
	$(CC) $(COPT) $(MYII) $(MKLFLAGS) -o $(SPD)/spd.o -c $(SPD)/spd.c


### ptb2 ###
$(PTB2)/ptb2_all.o: $(PTB2)/init_ptb2.o $(PTB2)/fret.o $(PTB2)/mrt.o $(PTB2)/opt_split.o
	ld -r -o $(PTB2)/ptb2_all.o $(PTB2)/init_ptb2.o $(PTB2)/fret.o $(PTB2)/mrt.o $(PTB2)/opt_split.o

$(PTB2)/init_ptb2.o: $(PTB2)/init_ptb2.c $(MYI_PTB2)/init_ptb2.h
	$(CC) $(COPT) $(MYII) $(MKLFLAGS) -o $(PTB2)/init_ptb2.o -c $(PTB2)/init_ptb2.c

$(PTB2)/fret.o: $(PTB2)/fret.c $(MYI_PTB2)/fret.h
	$(CC) $(COPT) $(MYII) $(MKLFLAGS) -o $(PTB2)/fret.o -c $(PTB2)/fret.c

$(PTB2)/mrt.o: $(PTB2)/mrt.c $(MYI_PTB2)/mrt.h
	$(CC) $(COPT) $(MYII) $(MKLFLAGS) -o $(PTB2)/mrt.o -c $(PTB2)/mrt.c

$(PTB2)/opt_split.o: $(PTB2)/opt_split.c $(MYI_PTB2)/opt_split.h
	$(CC) $(COPT) $(MYII) $(MKLFLAGS) -o $(PTB2)/opt_split.o -c $(PTB2)/opt_split.c


### trace ###
$(TRACE)/trace_all.o: $(TRACE)/init_tr.o $(TRACE)/calc_tr.o $(TRACE)/calc_tr_mrt_new.o
	ld -r -o $(TRACE)/trace_all.o $(TRACE)/init_tr.o $(TRACE)/calc_tr.o $(TRACE)/calc_tr_mrt_new.o

$(TRACE)/init_tr.o: $(TRACE)/init_tr.c $(MYI_TRACE)/init_tr.h
	$(CC) $(COPT) $(MYII) $(MKLFLAGS) -o $(TRACE)/init_tr.o -c $(TRACE)/init_tr.c

$(TRACE)/calc_tr.o: $(TRACE)/calc_tr.c $(MYI_TRACE)/calc_tr.h
	$(CC) $(COPT) $(MYII) $(MKLFLAGS) -o $(TRACE)/calc_tr.o -c $(TRACE)/calc_tr.c

$(TRACE)/calc_tr_mrt_new.o: $(TRACE)/calc_tr_mrt_new.c $(MYI_TRACE)/calc_tr_mrt_new.h
	$(CC) $(COPT) $(MYII) $(MKLFLAGS) -o $(TRACE)/calc_tr_mrt_new.o -c $(TRACE)/calc_tr_mrt_new.c


### dynamics ###
$(DYN)/dyn_all.o: $(DYN)/init_dyn.o $(DYN)/dyn_inco.o
	ld -r -o $(DYN)/dyn_all.o $(DYN)/init_dyn.o $(DYN)/dyn_inco.o

$(DYN)/init_dyn.o: $(DYN)/init_dyn.c $(MYI_DYN)/init_dyn.h
	$(CC) $(COPT) $(MYII) $(MKLFLAGS) -o $(DYN)/init_dyn.o -c $(DYN)/init_dyn.c

$(DYN)/dyn_inco.o: $(DYN)/dyn_inco.c $(MYI_DYN)/dyn_inco.h
	$(CC) $(COPT) $(MYII) $(MKLFLAGS) -o $(DYN)/dyn_inco.o -c $(DYN)/dyn_inco.c

