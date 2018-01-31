#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25
set output 'TiN_Alumina_BR1_Spectra.eps'
set xlabel 'Wavelength (nm)'
set ylabel 'Sectral Density (W/m^2)'
set pointsize 3
plot 'TiN_Alumina_BR1_spectra.txt' u ($1*1e9):4 w l lw 5 title 'Thermal Emission', \
'TiN_Alumina_BR1_spectra.txt' u ($1*1e9):5 w l lw 5 title 'Blackbody - 1500 K'
