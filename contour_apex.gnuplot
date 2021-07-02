
set size ratio 1 
set origin 0,0
set encoding utf8

#set style line 1 lc rgb '#000000' 
#set style line 2 lc rgb '#104983'
#set style line 3 lc rgb '#2c81d6'
#set style line 4 lc rgb '#7aafe3'
#set style line 5 lc rgb '#cccccc'
set style line 6 dt 3 lc rgb '#bbbbbb' # gridstyle
#set style increment user

set multiplot title aether_title


# -------------- plot grid and contour --------------------

set size ratio 0.5 

set key at screen 0.94,0.55
set xlabel aether_x_axis
set ylabel aether_y_axis
set xtics 24,-2,0 scale 0.75
set ytics -90,30,90 scale 0.5
set format y "%.0f°"
set grid ls 6



set view map
unset surface
set contour
set cntrparam levels discrete 1,2.3,6.2,11.8

set xrange [24:0]
set yrange [-90:90]
splot aether_data_file index 4 u (12/3.1416*$1):(180/3.1416*$2):3 w lines title '∆χ²'


set surface
unset contour


# -------------- plot symbols --------------------

# alway commented, to draw border over grid
#unset border # comment this line to check if both plots are aligned
unset xtics  
unset ytics  
unset xlabel
unset ylabel
unset grid

set size ratio 0.5 aether_size_x,aether_size_y
set origin aether_origin_x,aether_origin_y
set key at screen 0.94,0.85

set view equal xy
#set object ellipse center 12,0 size 24,180 fs empty border rgb "black"

plot  	for [i=0:1] aether_data_file index i using 1:2 w points pt (7-i) lc 1 title columnhead(4)



# -------------- plot star symbols --------------------


set key at screen 0.94,0.725

plot 	'-' u 1:2 w points pt 5 lc -1 title "Regulus",\
	'-' u 1:2 w points pt 13 lc -1 title "Denebola",\
	'-' u 1:2 w points pt 1 lc -1 title "Wolf 359"
10.13953	11.967	"Regulus"
e
11.81766	14.572	"Denebola"
e
10.9413	7.014	"Wolf 359"
e

unset multiplot		

