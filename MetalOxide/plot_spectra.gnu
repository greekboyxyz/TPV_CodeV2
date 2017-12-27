#!/usr/bin/gnuplot
  
set terminal postscript enhanced color 'Helvetica' 25
set output 'TiN_Spectrum_1.eps'
set xlabel 'Wavelength (nm)'
set ylabel 'Sectral Density (W/m^2 / nm)'
set pointsize 3
plot 'TiN_Alumina1_spectra.txt' u 1:4 w l lw 4 title 'Thermal Emission of TiN Emitter', \
'TiN_Alumina1_spectra.txt' u 1:5 w l lw 4 title 'Blackbody Spectrum at T=1500 K' 
