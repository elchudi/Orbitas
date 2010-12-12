#!/bin/bash
rm ener
rm outFey
gcc -Wall -c -ansi -pedantic orbitas.c 
gcc -Wall -c -ansi -pedantic vectorAlgebra.c
gcc -Wall -c -ansi -pedantic octantFunc.c
gcc  -o exe orbitas.o vectorAlgebra.o octantFunc.o -lm
./exe
gnuplot -persist orbitas.p
