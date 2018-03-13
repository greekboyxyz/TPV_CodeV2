#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25
set output 'TiN_Fit.eps'
set key top left
set ylabel 'Refractive Index'
set xlabel 'Wavelength (nm)'
plot 'TiN_Nari_Feb_9.text' u 1:2 w l lw 4 title '{/Symbol h} - Ellipsometry', \
'TiN_Spline.txt' u ($1*1e9):2 w l lw 4 dt 2 title '{/Symbol h} - FIT', \
'TiN_Nari_Feb_9.text' u 1:3 w l lw 4 title '{/Symbol k} - Ellipsometry', \
'TiN_Spline.txt' u ($1*1e9):3 w l lw 4 dt 2 title '{/Symbol k} - FIT'
