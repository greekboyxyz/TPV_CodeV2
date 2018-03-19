#!/usr/bin/gnuplot
reset
set yrange [0.3:0.65]
set xrange [0.09:0.2]
set terminal png
set xlabel 'TPV Efficiency'
set ylabel 'Spectral Efficiency'
set output 'MultiLayer_PO_SE_vs_ETA.png'
plot 'TiN_Alumina_BR_0ML_Pareto.txt' u 8:6 w p title '0 ML', \
'TiN_Alumina_BR_1ML_Pareto.txt' u 8:6 w p title '1 ML', \
#'TiN_Alumina_BR_2ML_Pareto.txt' u 8:6 w p title '2 ML', \
#'TiN_Alumina_BR_3ML_Pareto.txt' u 8:6 w p title '3 ML', \
#'TiN_Alumina_BR_4ML_Pareto.txt' u 8:6 w p title '4 ML', \
#'TiN_Alumina_BR_5ML_Pareto.txt' u 8:6 w p title '5 ML', \
#'W_Only.txt' u 4:2 w p title 'W Only'

<<<<<<< HEAD
set key top left
=======

reset 
>>>>>>> 623a9780826e2745891d211d1f8babdca850e69d
set yrange [0:8.2e4]
#set terminal postscript enhanced color 'Helvetica' 25
set terminal png
set output 'MultiLayer_Pareto.png'
set xlabel 'Spectral Efficiency'
set ylabel 'Sectral Density (W/m^2)'
set pointsize 3 
<<<<<<< HEAD
plot 'TiN_Alumina_BR_1ML_Pareto.txt' u 6:7 w p title 'TiN - 1 ML', \
'TiN_Alumina_BR_2ML_Pareto.txt'u 6:7 w p title 'TiN - 2 ML', \
'TiN_Alumina_BR_3ML_Pareto.txt' u 6:7 w p title 'TiN - 3 ML', \
'TiN_Alumina_BR_4ML_Pareto.txt' u 6:7 w p title 'TiN - 4 ML', \
'TiN_Alumina_BR_5ML_Pareto.txt' u 6:7 w p title 'TiN - 5 ML', \
'TiN_Alumina_BR_6ML_Pareto.txt' u 6:7 w p title 'TiN - 6 ML', \
'TiN_Alumina_BR_7ML_Pareto.txt' u 6:7 w p title 'TiN - 7 ML', \
'TiN_Alumina_BR_8ML_Pareto.txt' u 6:7 w p title 'TiN - 8 ML', \
'TiN_Alumina_BR_9L_Pareto.txt' u 6:7 w p title 'TiN - 9 L', \
'TiN_Alumina_BR_10ML_Pareto.txt' u 6:7 w p title 'TiN - 10 ML', \
=======
plot 'TiN_Alumina_BR_0ML_Pareto.txt' u 6:7 w p title '0 ML', \
'TiN_Alumina_BR_1ML_Pareto.txt' u 6:7 w p title '1 ML', \
#'TiN_Alumina_BR_2ML_Pareto.txt' u 6:7 w p title '2 ML', \
#'TiN_Alumina_BR_3ML_Pareto.txt' u 6:7 w p title '3 ML', \
#'TiN_Alumina_BR_4ML_Pareto.txt' u 6:7 w p title '4 ML', \
#'TiN_Alumina_BR_5ML_Pareto.txt' u 6:7 w p title '5 ML', \
#'W_Only.txt' u 2:3 w p title 'W Only' 

set xlabel 'TPV Efficiency'
set output 'MultiLayer_PO_SD_vs_ETA.png'
plot 'TiN_Alumina_BR_0ML_Pareto.txt' u 8:7 w p title '0 ML', \
'TiN_Alumina_BR_1ML_Pareto.txt' u 8:7 w p title '1 ML', \
<<<<<<< HEAD
'TiN_Alumina_BR_2ML_Pareto.txt' u 8:7 w p title '2 ML', \
'TiN_Alumina_BR_3ML_Pareto.txt' u 8:7 w p title '3 ML', \
'TiN_Alumina_BR_4ML_Pareto.txt' u 8:7 w p title '4 ML', \
'TiN_Alumina_BR_5ML_Pareto.txt' u 8:7 w p title '5 ML', \
'W_Only.txt' u 4:3 w p title 'W Only'
>>>>>>> 623a9780826e2745891d211d1f8babdca850e69d
=======
#'TiN_Alumina_BR_2ML_Pareto.txt' u 8:7 w p title '2 ML', \
#'TiN_Alumina_BR_3ML_Pareto.txt' u 8:7 w p title '3 ML', \
#'TiN_Alumina_BR_4ML_Pareto.txt' u 8:7 w p title '4 ML', \
#'TiN_Alumina_BR_5ML_Pareto.txt' u 8:7 w p title '5 ML', \
#'W_Only.txt' u 4:3 w p title 'W Only'
>>>>>>> de96695c3370aef97f8f422a3a5e8d0f1f6b6423
