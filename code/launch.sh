#!/bin/bash

rm *.png

./build/MOPSI 100 0.003

gnuplot plot.gnu

rm *.csv
