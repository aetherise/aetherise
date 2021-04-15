#!/bin/bash


if [[ $# -eq 0 ]]; then
	echo  "$0 <filename> [title] [output-filename.svg]"
	exit 1
fi



if [ "$3" = "" ]; then
	term="qt 0 size 640,300"
else	
	if [ "$4" = "" ]; then
		# fontsize 6pt and width 8.5cm on a 96dpi display
		# 8*72/96=6pt, 320/96*2.54=8.5cm
		term="svg font 'Sans,8' size 320,150 linewidth 0.5 background '#ffffff'; set pointsize 0.5; set output '$3'"
	else
		term="$4""; set output '$3'"
	fi

fi


title=$2

xaxis="Azimuth Index"
yaxis="Fringes"
data_title="Data"
theorie_title="Theory"
decimalsign="."


gnuplot -persist -e "aether_data_file='$1'; aether_x_axis='$xaxis';aether_y_axis='$yaxis';aether_data_title='$data_title';aether_theorie_title='$theorie_title';aether_title='$title';set decimalsign '$decimalsign';set terminal $term" fringes.gnuplot

