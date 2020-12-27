#!/bin/bash


if [[ $# -eq 0 ]]; then
	echo  "$0 <filename> [title]"
	exit 1
fi


title=$2

if [[ $LANG == de* ]]
then	
	xaxis="Sternzeit / h"
	yaxis="Signalamplitude"
	data_title="Daten"
	theorie_title="Theorie"
	model_title="Modell"
else	
	xaxis="Sidereal Time / h"
	yaxis="Signal Amplitude"
	data_title="Data"
	theorie_title="Theory"
	model_title="Model"
fi

gnuplot -persist -e "aether_data_file='$1'; aether_x_axis='$xaxis';aether_y_axis='$yaxis';aether_data_title='$data_title';aether_theorie_title='$theorie_title';aether_model_title='$model_title';aether_title='$title'" sidereal.gnuplot

