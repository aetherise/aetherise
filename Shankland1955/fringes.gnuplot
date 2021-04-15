set encoding utf8

set title aether_title
set xlabel aether_x_axis
set ylabel aether_y_axis
unset key

set datafile missing "?"

set xtics 1,1,9
set ytics -0.05,0.05,0.05

set xrange [0:10]
set yrange [-0.08:0.08]

plot 	aether_data_file using 1:2:3 title aether_data_title w yerrorbars lc 1 pt 7,\
	'' w lines lc 1 dt 2 notitle,\
	aether_data_file using 1:($1<10 ? $4*2 : 1/0) title aether_theorie_title w linespoints lc 2 pt 7,\
	aether_data_file using 1:5 notitle w lines lc 3


