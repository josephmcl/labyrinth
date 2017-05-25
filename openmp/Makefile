CC=cc
ICC=icc
CXX=icpc
CFLAGS=-g -O3
CDEBUG=-W -Wall -Wpedantic
CXXFLAGS=-g -O3

all: serial omp

serial: serial.c
	$(CC) $(CFLAGS) $(CDEBUG) -o serial serial.c -lm

omp: omp.c
	$(CC) $(CFLAGS) $(CDEBUG) -o omp omp.c -lm -fopenmp 

clean: 
	rm serial omp
