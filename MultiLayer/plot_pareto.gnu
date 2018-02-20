#!/usr/bin/gnuplot

set key top lef
set yrange [0:8.2e4]
set terminal postscript enhanced color 'Helvetica' 25
set output 'MultiLayer_Pareto.eps'
set xlabel 'Spectral Efficiency'
set ylabel 'Sectral Density (W/m^2)'
set pointsize 3 
plot 'HfN_Alumina_BR_Pareto.txt' u 6:7 w p title 'HfN - Full Structure', \
'TiN_Alumina_BR_Pareto.txt' u 6:7 w p title 'TiN - Full Structure', \
'Scan_HfN_MultiLayer.txt' u 7:8 w p title 'HfN - Multi-layer Only', \
'Scan_TiN_MultiLayer.txt' u 7:8 w p title 'TiN - Multi-layer Only'
