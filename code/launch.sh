#!/bin/bash

rm *.csv
rm *.png

./build/MOPSI

gnuplot plot.gnu
