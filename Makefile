CC=cc
ICC=icc
CXX=icpc
CFLAGS=-g -O3 
CXXFLAGS=-g -O3

all: serial omp

serial: serial.c
	$(CC) $(CFLAGS) -o serial serial.c -lm

omp:
	$(CC) $(CFLAGS) -o omp omp.c -lm -fopenmp 

clean: 
	rm serial omp
