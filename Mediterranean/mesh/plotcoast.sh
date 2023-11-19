#!/bin/bash

infile=coast.dat
outps=coast.ps

gnuplot <<EOF
set terminal postscript
set output "$outps"
plot "$infile" u 1:2 w p

EOF
gv $outps &
open $outps
