#!/bin/bash


if [[ $# -eq 0 ]]; then
	echo  "$0 <filename> [title] [output-filename] [gnuplot-terminal]"
	exit 1
fi


title=$2


if [ "$3" = "" ]; then
	term="qt 0 font 'Sans,9'"
else
	if [ "$4" = "" ]; then
		# fontsize 8pt and width 15.9cm on a 96dpi display
		# 10.67*72/96=8pt, 600/96*2.54=15.9cm
		term="svg font 'Sans,8' size 320,200 linewidth 0.5 background '#ffffff'; set output '$3'"
	else
		term="$4""; set output '$3'"
	fi
fi


if [[ $LANG == de* ]]
then
	xaxis="Sternzeit / h"
	yaxis="Signalamplitude / λ"
	data_title="Daten"
	theorie_title="Theorie"
	decimalsign=","
else
	xaxis="Sidereal Time / h"
	yaxis="Signal Amplitude / λ"
	data_title="Data"
	theorie_title="Theory"
	decimalsign="."
fi

gnuplot -persist -e "aether_data_file='$1'; aether_x_axis='$xaxis';aether_y_axis='$yaxis';aether_data_title='$data_title';aether_theorie_title='$theorie_title';aether_title='$title';set terminal $term; set decimalsign '$decimalsign'" theory_amp.gnuplot

