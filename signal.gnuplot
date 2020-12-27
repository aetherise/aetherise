set encoding utf8

set title aether_title
set xlabel aether_x_axis
set ylabel aether_y_axis
set key box

set datafile missing "?"

set xtics 1,4,17
set xrange [0:18]
set yrange [*<-0.04:0.04<*]

plot	aether_data_file using 1:2:3 title aether_data_title w yerrorbars lc 1 pt 7,\
	'' w lines lc 1 dt 2 notitle,\
	aether_data_file using 1:4 title aether_theorie_title w lines lc 2,\
	aether_data_file using 1:5 notitle w lines lc 3


