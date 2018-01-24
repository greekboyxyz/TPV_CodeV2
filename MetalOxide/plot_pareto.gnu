#!/usr/bin/gnuplot

set key bottom left
set terminal postscript enhanced color 'Helvetica' 25
set output 'Pareto.eps'
set xlabel 'Spectral Efficiency'
set ylabel 'Sectral Density (W/m^2)'
set pointsize 3 
plot 'Pareto.txt' u 6:7 w p title 'Pareto', /
'Ru_Alumina_Pareto.txt' u 6:7 w p title 'Ru', /
'Rh_Alumina_Pareto.txt' u 6:7 w p title 'Rh', /
'Re_Alumina_Pareto.txt' u 6:7 w p title 'Re', /
'HfN_Alumina_Pareto.txt' u 6:7 w p title 'HfN'
