#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25
set output 'MultiLayer_Spectra.eps'
set xlabel 'Wavelength (nm)'
set ylabel 'Emissivity'
set key bottom right
set pointsize 3
plot 'Re_Alumina_BR1_spectra.txt' u ($1*1e9):3 w l lw 5 title 'Re', \
'Rh_Alumina_BR1_spectra.txt' u ($1*1e9):3 w l lw 5 title 'Rh', \
'Ru_Alumina_BR1_spectra.txt' u ($1*1e9):3 w l lw 5 title 'Ru', \
'HfN_Alumina_BR1_spectra.txt' u ($1*1e9):3 w l lw 5 title 'HfN'
