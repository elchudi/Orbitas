#!/bin/bash
gcc -Wall -c -ansi -pedantic orbitas.c 
gcc orbitas.o -o exe -lm
./exe 
