set encoding utf8

set title aether_title
set xlabel aether_x_axis
set ylabel aether_y_axis
set zlabel aether_z_axis rotate by 90
set key off # does not work with version 5.3 on windows

set xrange [0:24]
set xtics 0,4,24
set yrange [0:18]
set ytics 1,4,17
set zrange [*<-0.02:0.02<*]
set ztics 0.01
set view 33,136

splot 	aether_data_file using 6:1:2 title aether_data_title w lines lc 1,\
	aether_data_file using 6:1:4 title aether_theorie_title w lines lc 2


