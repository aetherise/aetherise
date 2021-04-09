#!/bin/bash



if [[ $# -eq 0 ]]; then
	echo  "$0 <filename> [title] [output-filename] [gnuplot-terminal]"
	exit 1
fi

title=$2

asize=0.78
origin_x=0.111
origin_y=0.125

if [ "$3" = "" ]; then
	term="qt 0 font 'Sans,9' size 1200,600"
else
	if [ "$4" = "" ]; then		
		asize=0.746
		origin_x=0.139
		origin_y=0.154
		# fontsize 8pt and width 15.9cm on a 96dpi display
		# 10.67*72/96=8pt, 768/96*2.54=20.3cm
		term="svg font 'Sans,10.67' size 768,320 linewidth 0.5 background '#ffffff'; set pointsize 0.5; set output '$3'"
		#term="pngcairo font 'Sans,9' size 1200,600; set output '$3'"
	else
		term="$4""; set output '$3'"
	fi
fi




if [[ $LANG == de* ]]
then
	echo "WARNUNG: Multiplots müssen von Hand ausgerichtet werden."
	xaxis="Rektaszension 𝛼 / h"
	yaxis="Deklination 𝛿"
	decimalsign=","
else	
	echo "WARNING: Multiplots have to be aligned manually."
	xaxis="Right Ascension 𝛼 / h"	
	yaxis="Declination 𝛿"
	decimalsign="."
fi

gnuplot -persist -e "aether_size_x=$asize;aether_size_y=$asize;aether_origin_x=$origin_x;aether_origin_y=$origin_y;aether_data_file='$1';aether_x_axis='$xaxis';aether_y_axis='$yaxis';aether_title='$title';set terminal $term;set decimalsign '$decimalsign'" contour_apex.gnuplot

