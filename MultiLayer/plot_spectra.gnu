#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25
set output 'MultiLayer_Spectra.eps'
set xlabel 'Wavelength (nm)'
set ylabel 'Spectral Irradiance (W / m / nm / sr)'
set key top right
set yrange [0:5e10]
set pointsize 3
plot 'TiN_Alumina_BR634_spectra.txt' u ($1*1e9):4 w l lw 5 title 'TiN, {/Symbol h} = 54%, P = 2.9e+04 W/m^2', \
'HfN_Alumina_BR347_spectra.txt' u ($1*1e9):4 w l lw 5 title 'HfN, {/Symbol h} = 51%, P= 2.6e+04 W/m^2', \
'TiN_Alumina_BR634_spectra.txt' u ($1*1e9):5 w l lw 5 lt -1 title 'B({/Symbol l},T=1473 K)'


