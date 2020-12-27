#!/bin/bash


if [[ $# -eq 0 ]]; then
	echo  "$0 <filename> [title]"
	exit 1
fi

title=$2

if [[ $LANG == de* ]]
then
	xaxis="Azimut"
	yaxis="Abstand in 1/10 Streifen"
else
	xaxis="Azimuth"
	yaxis="Distance in 1/10 fringe"
fi

gnuplot -persist -e "aether_data_file='$1'; aether_x_axis='$xaxis';aether_y_axis='$yaxis';aether_title='$title'" raw.gnuplot

