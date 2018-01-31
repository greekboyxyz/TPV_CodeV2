#!/usr/bin/gnuplot

set key bottom left
set terminal postscript enhanced color 'Helvetica' 25
set output 'HfN_Alumina_BR_Pareto.eps'
set xlabel 'Spectral Efficiency'
set ylabel 'Sectral Density (W/m^2)'
set pointsize 3 
plot 'HfN_Alumina_BR_Pareto.txt' u 6:7 w p title 'HfN'
