set encoding utf8

set title aether_title
set xlabel aether_x_axis
set ylabel aether_y_axis
set key off

set datafile missing "?"

set xtics out nomirror
#set boxwidth 0.7 absolute # does not work
#set style fill solid 1.0
set yrange [0:0.04<*]
set xrange [0:*]

plot	aether_data_file using 1:2 title aether_data_title w impulses lw 1


