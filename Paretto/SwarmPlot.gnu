#!/usr/bin/gnuplot
#set terminal png
set terminal postscript enhanced color 'Helvetica' 23 

set pointsize 2
#set yrange [0:1.4e6]
set key top left
f = 100*100
set yrange [0:800000/f]
set output 'SwarmNew.eps'
set xlabel 'Spectral Efficiency'
set ylabel 'Useful Power Density (W/cm^2)'
plot 'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/SWARM_18_MG_Alloy_SiO2_TiO2_Static.txt' u 6:($7/f) w p title 'BR + W + Alloy', \
'Substrate_and_BraggReflector/Temperature_Scan/SWARM_18_Vendor_BR.txt' u 7:($8/f) w p title 'BR + W', \
'Substrate/Temperature_Scan/W_FOM_vs_T.txt' u 2:($3/f) w p title 'W', \
'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/Pareto_18_MG_SiO2_TiO2_Static_Temp.txt' u 6:($7/f) w p title 'Pareto Front', \
'Step_FOM.txt' u 2:($3/f) w l lw 4 lt 16 dt 2 title 'Step Function', \
'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/Select.txt' u 6:($7/f) w p lt 16 pt 5 title 'Target BR + W + Alloy', \
'Substrate_and_BraggReflector/Temperature_Scan/Target_W_BR.txt' u 7:($8/f) w p lt 16 pt 6 title 'Target BR + W + Alloy' 

set output 'SwarmNew2.eps'
set xlabel 'Spectral Efficiency'
set ylabel 'Useful Power Density (W/cm^2)'
plot 'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/SWARM_18_MG_Alloy_SiO2_TiO2_Static.txt' u 6:($7/f) w p title 'BR + W + Alloy', \
'Substrate_and_BraggReflector/Temperature_Scan/SWARM_18_Vendor_BR.txt' u 7:($8/f) w p title 'BR + W', \
'Substrate/Temperature_Scan/W_FOM_vs_T.txt' u 2:($3/f) w p title 'W', \
'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/Pareto_18_MG_SiO2_TiO2_Static_Temp.txt' u 6:($7/f) w p title 'Pareto Front', \
'Step_FOM.txt' u 2:($3/f) w l lw 4 lt 16 dt 2 title 'Step Function', \
'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/Select_NonOpt.txt' u 1:($2/f) w p lt 16 pt 5 title 'Target BR + W + Alloy', \
'/home/jay/SOFTWARE/TPV_Code/Spectral_Plot/Spectral_Plot_Vendor_Actual_BR/W_BR_Only/Select_NonOpt_W_Br_Only.txt' u 1:($2/f) w p lt 16 pt 6 title 'Target BR + W + Alloy'


#set yrange [0:800000]
#set output 'SwarmCompare_MG_to_BG.eps'
#set xlabel 'Spectral Efficiency'
#set ylabel 'Spectral Density (W/m^2)'
#plot 'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/Swarm_18_MG_Alloy_Vendor_BR.txt' u 9:10 w p title 'Bulk Alloy', \
#'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/Swarm_Bruggenman_Alloy_Vendor_BR_Temp.txt' u 9:10 w p lt 5 pt 5 title 'ALD Alloy', \
#'Substrate_and_BraggReflector/Temperature_Scan/SWARM_18_Vendor_BR.txt' u 7:8 w p lt 2 pt 2 title 'BR + W', \
#'Substrate/Temperature_Scan/W_FOM_vs_T.txt' u 2:3 w p lt 3 pt 3 title 'W', \
#'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/FabricationTarget.txt' u 1:2 w p lt 16 pt 5 title 'Fab Alloy + BR + W', \
#'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/FabricationTarget_BR_Only.txt' u 1:2 w p lt 17 pt 6 title 'Fab BR + W'




#set yrange [0:500000]
#set output 'SwarmCompare_HfO2.eps'
#set xlabel 'Spectral Efficiency'
#set ylabel 'Spectral Density (W/m^2)'
#plot 'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/SWARM_18.txt' u 9:10 w p title 'BR + W + Alloy', \
#'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/SWARM_18_Alumina_HfO2.txt' u 9:10 w p title 'Al_2O_3/HfO_2 + W + Alloy', \
#'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/Pareto_18.txt' u 6:7 w p title 'Pareto Front', \
#'Substrate_and_BraggReflector_and_Alloy/W_alloy_SiO2_TiO2_W_Sub/Temperature_Scan/Pareto_18_Alumina_HfO2.txt' u 6:7 w p title 'Al_2O_3/HfO_2 Pareto Front'
##
