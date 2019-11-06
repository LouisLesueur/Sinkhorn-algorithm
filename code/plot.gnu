set term png
set output 'simplex1.png'
plot 'simplex1.csv'
set output 'simplex2.png'
plot 'simplex2.csv'
set output 'pi.png'
plot 'pi.csv' matrix with image
