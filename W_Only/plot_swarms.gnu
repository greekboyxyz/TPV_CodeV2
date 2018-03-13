#!/usr/bin/gnuplot

set terminal png
set output 'PowerDensity_vs_EtaTPV_1ML.png'
set xlabel 'TPV Efficiency'
set ylabel 'Sectral Density (W/m^2)'
set pointsize 3 
plot 'Scan_1ML.txt' u 9:8 w p title '1 ML', \
'TiN_Alumina_BR_1ML_Pareto.txt' u 8:7 w p title 'Pareto Front'

set output 'PowerDensity_vs_EtaTPV_2ML.png'
plot 'Scan_2ML.txt' u 9:8 w p title '2 ML', \
'TiN_Alumina_BR_2ML_Pareto.txt' u 8:7 w p title 'Pareto Front'

set output 'PowerDensity_vs_EtaTPV_3ML.png'
plot 'Scan_3ML.txt' u 9:8 w p title '3 ML', \
'TiN_Alumina_BR_3ML_Pareto.txt' u 8:7 w p title 'Pareto Front'

set output 'PowerDensity_vs_EtaTPV_4ML.png'
plot 'Scan_4ML.txt' u 9:8 w p title '4 ML', \
'TiN_Alumina_BR_4ML_Pareto.txt' u 8:7 w p title 'Pareto Front'

set output 'PowerDensity_vs_EtaTPV_5ML.png'
plot 'Scan_5ML.txt' u 9:8 w p title '5 ML', \
'TiN_Alumina_BR_5ML_Pareto.txt' u 8:7 w p title 'Pareto Front'

set output 'SpectralEfficiency_vs_EtaTPV_1ML.png'
set xlabel 'TPV Efficiency'
set ylabel 'Sectral Efficiency'
set pointsize 3
plot 'Scan_1ML.txt' u 9:7 w p title '1 ML', \
'TiN_Alumina_BR_1ML_Pareto.txt' u 8:6 w p title 'Pareto Front'

set output 'SpectralEfficiency_vs_EtaTPV_2ML.png'
plot 'Scan_2ML.txt' u 9:7 w p title '2 ML', \
'TiN_Alumina_BR_2ML_Pareto.txt' u 8:6 w p title 'Pareto Front'

set output 'SpectralEfficiency_vs_EtaTPV_3ML.png'
plot 'Scan_3ML.txt' u 9:7 w p title '3 ML', \
'TiN_Alumina_BR_3ML_Pareto.txt' u 8:6 w p title 'Pareto Front'

set output 'SpectralEfficiency_vs_EtaTPV_4ML.png'
plot 'Scan_4ML.txt' u 9:7 w p title '4 ML', \
'TiN_Alumina_BR_4ML_Pareto.txt' u 8:6 w p title 'Pareto Front'

set output 'SpectralEfficiency_vs_EtaTPV_5ML.png'
plot 'Scan_5ML.txt' u 9:7 w p title '5 ML', \
'TiN_Alumina_BR_5ML_Pareto.txt' u 8:6 w p title 'Pareto Front'
