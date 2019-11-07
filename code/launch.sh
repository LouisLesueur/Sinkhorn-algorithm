#!/bin/bash

rm *.png

./build/MOPSI

gnuplot plot.gnu

rm *.csv
