CC=cc
ICC=icc
CXX=icpc
CFLAGS=-g -O3 
CXXFLAGS=-g -O3

all: serial omp

serial: l.c
	$(CC) $(CFLAGS) -o serial l.c -lm

omp:
	$(CC) $(CFLAGS) -o omp omp.c -lm -fopenmp 

clean: 
	rm serial omp