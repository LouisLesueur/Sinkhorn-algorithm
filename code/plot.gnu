set term png
set key box lt -1 lw 2
set output 'simplex1.png'
plot 'simplex1.csv' with linespoints
set output 'simplex2.png'
plot 'simplex2.csv' with linespoints
set xrange [*:*]
set yrange [*:*] reverse
set output 'pi.png'
set palette defined (0 "#ffffff", 1 "#000000" )
set cbrange [0:1]
plot 'pi.csv' matrix with image
