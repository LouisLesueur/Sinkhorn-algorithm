set term png
set output 'simplex1.png'
plot 'simplex1.csv'
set output 'simplex2.png'
plot 'simplex2.csv'
set output 'pi.png'
set palette defined (0 "#ffffff", 1 "#000000" )
plot 'pi.csv' matrix with image
