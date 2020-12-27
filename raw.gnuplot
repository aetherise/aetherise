
set title aether_title
set xlabel aether_x_axis
set ylabel aether_y_axis
set key box outside

set datafile missing "?"

set xtics 1,4,17
set xrange [0:18]

plot for [i=0:*] aether_data_file index i using 1:2 w linespoints pt i+1 lc i+1 title columnheader(2)

