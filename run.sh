#!/bin/bash
rm ener
rm outFey
gcc -Wall -c -ansi -pedantic orbitas.c 
gcc orbitas.o -o exe -lm
./exe 
