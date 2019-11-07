#!/bin/bash

rm *.png

#./build/MOPSI 100 0.003
./build/MOPSI $1 $2
gnuplot plot.gnu

rm *.csv
