#!/usr/bin/gnuplot

set terminal png
set xlabel 'TPV Efficiency'
set ylabel 'Spectral Efficiency'
set output 'MultiLayer_PO_SE_vs_ETA.png'
plot 'TiN_Alumina_BR_1ML_Pareto.txt' u 8:6 w p title '1 ML', \
'TiN_Alumina_BR_2ML_Pareto.txt' u 8:6 w p title '2 ML', \
'TiN_Alumina_BR_3ML_Pareto.txt' u 8:6 w p title '3 ML', \
'TiN_Alumina_BR_4ML_Pareto.txt' u 8:6 w p title '4 ML', \
'TiN_Alumina_BR_5ML_Pareto.txt' u 8:6 w p title '5 ML', \
'TiN_Alumina_BR_6ML_Pareto.txt' u 8:6 w p title '6 ML'

set yrange [0:8.2e4]
#set terminal postscript enhanced color 'Helvetica' 25
set terminal png
set output 'MultiLayer_Pareto.png'
set xlabel 'Spectral Efficiency'
set ylabel 'Sectral Density (W/m^2)'
set pointsize 3 
plot 'TiN_Alumina_BR_1ML_Pareto.txt' u 6:7 w p title '1 ML', \
'TiN_Alumina_BR_2ML_Pareto.txt' u 6:7 w p title '2 ML', \
'TiN_Alumina_BR_3ML_Pareto.txt' u 6:7 w p title '3 ML', \
'TiN_Alumina_BR_4ML_Pareto.txt' u 6:7 w p title '4 ML', \
'TiN_Alumina_BR_5ML_Pareto.txt' u 6:7 w p title '5 ML', \
'TiN_Alumina_BR_6ML_Pareto.txt' u 6:7 w p title '6 ML'

set xlabel 'TPV Efficiency'
set output 'MultiLayer_PO_SD_vs_ETA.png'
plot 'TiN_Alumina_BR_1ML_Pareto.txt' u 8:7 w p title '1 ML', \
'TiN_Alumina_BR_2ML_Pareto.txt' u 8:7 w p title '2 ML', \
'TiN_Alumina_BR_3ML_Pareto.txt' u 8:7 w p title '3 ML', \
'TiN_Alumina_BR_4ML_Pareto.txt' u 8:7 w p title '4 ML', \
'TiN_Alumina_BR_5ML_Pareto.txt' u 8:7 w p title '5 ML', \
'TiN_Alumina_BR_6ML_Pareto.txt' u 8:7 w p title '6 ML'

