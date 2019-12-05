#!/bin/bash

#./build/MOPSI 100 0.003
./build/MOPSI $1 $2
gnuplot -e "N=$1" plot.gnu

rm *.csv
