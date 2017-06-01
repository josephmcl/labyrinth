#!/bin/bash
	echo Parallel\n;
        for i in `seq 1 20`;
        do
            ./omp 1000
        done  
	echo Serial\n;
        for i in `seq 1 20`;
        do
            ./omp 1000
        done 
