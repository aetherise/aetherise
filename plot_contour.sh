#!/bin/bash


if [[ $# -eq 0 ]]; then
	echo  "$0 <filename> [title] [output-filename] [gnuplot-terminal]"
	exit 1
fi



title=$2


asize=0.78
origin_x=0.119
origin_y=0.122

if [ "$3" = "" ]; then
	term="qt 0 font 'Sans,9' size 1000,600"
else	
	if [ "$4" = "" ]; then
		# gnuplot 5.2
		asize=0.77
		origin_x=0.129
		origin_y=0.135
		# fontsize 8pt and width 15.9cm on a 96dpi display
		# 10.67*72/96=8pt, 640/96*2.54=16.9cm
		term="svg size 640,420 font 'Sans,10.67' enhanced linewidth 0.5 background '#ffffff'; set pointsize 0.5; set output '$3'"
		#term="pngcairo font 'Sans,9' size 1000,600; set output '$3'"
	else
		term="$4""; set output '$3'"
	fi

fi




if [[ $LANG == de* ]]
then
	echo "WARNUNG: Multiplots müssen von Hand ausgerichtet werden."
	xaxis="v / (km/s)"
	yaxis="Rektaszension α / h"
	y2axis="Deklination δ"
	decimalsign=","
else	
	echo "WARNING: Multiplots have to be aligned manually."
	xaxis="v / (km/s)"
	yaxis="Right Ascension α / h"
	y2axis="Declination δ"
	decimalsign="."
fi

gnuplot -persist -e "aether_size=$asize;aether_origin_x=$origin_x;aether_origin_y=$origin_y;aether_data_file='$1';aether_x_axis='$xaxis';aether_y_axis='$yaxis';aether_y2_axis='$y2axis';aether_title='$title';set terminal $term;set decimalsign '$decimalsign'" contour.gnuplot

