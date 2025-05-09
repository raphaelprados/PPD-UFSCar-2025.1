#!/bin/sh

size=400

gcc seq.c -o seq
gcc seq_lin.c -o seq_lin
gcc seq_ptr.c -o seq_ptr
gcc seq_linptr.c -o seq_linptr

./seq "$size"
./seq_lin "$size"
./seq_ptr "$size"
./seq_linptr "$size"

rm seq seq_lin seq_ptr seq_linptr grid_laplace.txt

