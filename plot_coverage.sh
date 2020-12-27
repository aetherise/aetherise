#!/bin/bash

if [[ $# -eq 0 ]]; then
	echo  "$0 <filename> [title] [output-filename.svg]"
	exit 1
fi

title=$2

if [ "$3" = "" ]; then
	term="qt 0 font 'Sans,9' size 800,600"
else	
	if [ "$4" = "" ]; then		
		# fontsize 8pt and width 15.9cm on a 96dpi display
		# 10.67*72/96=8pt, 640/96*2.54=16.9cm
		term="svg size 768,320 font 'Sans,10.67' enhanced linewidth 0.5 background '#ffffff'; set pointsize 0.5; set output '$3'"		
	else
		term="$4""; set output '$3'"
	fi

fi



if [[ $LANG == de* ]]
then
	xaxis="Signalamplitude / λ"
	yaxis="Sternzeit / h"
	title_amp_w="KHS-Dipol, Winter"
	title_amp_s="KHS-Dipol, Sommer"
	title_cov_w="Überdeckung, Winter"
	title_cov_s="Überdeckung, Sommer"	
	decimalsign=","
else	
	xaxis="Signal Amplitude / λ"
	yaxis="Sidereal Time / h"
	title_amp_w="CMB-Dipole, Winter"
	title_amp_s="CMB-Dipole, Summer"
	title_cov_w="Coverage, Winter"
	title_cov_s="Coverage, Summer"
	decimalsign="."
fi

gnuplot -persist -e "aether_data_file='$1'; aether_x_axis='$xaxis';aether_y_axis='$yaxis';aether_title='$title'; title_amp_w='$title_amp_w';title_amp_s='$title_amp_s';title_cov_s='$title_cov_s';title_cov_w='$title_cov_w'; set terminal $term;set decimalsign '$decimalsign'" coverage.gnuplot

