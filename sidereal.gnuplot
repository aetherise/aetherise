set encoding utf8

set title aether_title
set xlabel aether_x_axis
set ylabel aether_y_axis
set key box

set y2range [-10:30]

set xtics 0,4,24
set xrange [0:24]

plot	aether_data_file using 1:2 title aether_data_title lc 1 pt 7,\
	'' w lines lc 1 dt 2 notitle,\
	aether_data_file using 1:3 title aether_theorie_title lc 2 pt 1,\
	aether_data_file using 1:3 w lines lc 2 notitle,\
	aether_data_file using 1:4 title aether_model_title w lines lc 3,\
	aether_data_file using 1:5 title "TD" axes x1y2 w lines lc 4
	


