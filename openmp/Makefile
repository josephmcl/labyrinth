CC=cc
ICC=icc
CXX=icpc
CFLAGS=
CDEBUG=-W -Wall
CXXFLAGS=-g -O3

all: serial omp

serial: serial.c
	$(CC) -O3$(CFLAGS) $(CDEBUG) -o serial serial.c -lm -fopenmp 

omp: omp.c
	$(CC) -O0 -fno-unroll-loops  $(CDEBUG) -o omp omp.c -lm -fopenmp 

clean: 
	rm serial omp
