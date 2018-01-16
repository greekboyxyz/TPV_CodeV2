#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25
set output 'HfN_Alumina1_Spectra.eps'
set xlabel 'Wavelength (nm)'
set ylabel 'Sectral Density (W/m^2)'
set pointsize 3
plot 'HfN_Alumina1_spectra.txt' u ($1*1e9):4 w l lw 5 title 'Thermal Emission', \
'HfN_Alumina1_spectra.txt' u ($1*1e9):5 w l lw 5 title 'Blackbody - 1500 K'
