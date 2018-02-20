#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25
set xlabel 'Wavelength (nm)'
set ylabel 'Spectral Irradiance (W / m / nm / sr)'
set key top right
set yrange [0:5e10]
set pointsize 3
set output 'MultiLayer_Only_Spectra.eps'
plot 'TiN_Alumina_BR1_MLspectra.txt' u ($1*1e9):4 w l lw 5 title 'TiN, {/Symbol h} = 31%, P = 3.3e+04 W/m^2', \
'HfN_Alumina_BR9_MLspectra.txt' u ($1*1e9):4 w l lw 5 title 'HfN, {/Symbol h} = 32%, P = 2.85e+04 W/m^2', \
'TiN_Alumina_BR1_MLspectra.txt' u ($1*1e9):5 w l lw 5 lt -1 title 'B({/Symbol l},T=1473 K)'

