
set title aether_title
set xlabel aether_x_axis
set ylabel aether_y_axis
set key box

set datafile missing "?"

set xtics 0,4,24
set xrange [0:24]
set yrange [*<-0.03:0.03<*]
set ytics 0.01

plot	aether_data_file using 6:4 title aether_theorie_title w lines lc 2


