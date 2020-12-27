
set size ratio 1
set origin 0,0
set termoption enhanced
set encoding utf8

set multiplot title aether_title


# -------------- plot contour --------------------


set key at screen 0.95,0.65

unset border # comment this line to check if both plots are aligned
unset xtics  # comment this line to check if both plots are aligned
unset ytics  # comment this line to check if both plots are aligned
set view map
unset surface
set contour
set cntrparam levels discrete 1,2.3,6.2,11.8

set xrange [-550:550]
set yrange [-550:550]
splot aether_data_file index 2 u (cos($2)*($1/1000)):(sin($2)*($1/1000)):3 w lines title '∆χ²'
unset key
splot aether_data_file index 3 u (cos($2)*($1/1000)):(sin($2)*($1/1000)):3 w lines notitle 


set surface
unset contour


# -------------- plot grid and symbols --------------------


#set xlabel aether_x_axis
set label aether_x_axis at graph 1,0.54 right
set ylabel aether_y_axis
set y2label aether_y2_axis

set size ratio 1 aether_size,aether_size
set origin aether_origin_x,aether_origin_y
set key at screen 0.95,0.85

set style line 6 dt 3 lc rgb '#cccccc' # gridstyle

set polar 
set angles degrees
set grid ls 6 polar 15
unset border
unset xtics
unset ytics
unset raxis
set yzeroaxis lt -1
set rrange [0:550]
set view equal xy
set rtics axis scale 0.5 ("100" 100, "200" 200, "" 300, "" 400, "500" 500)
set label "75°" at graph 0.61,0.95
set label "45°" at graph 0.825,0.825 
set label "15°" at graph 0.94,0.625
set label "-15°" at graph 0.94,0.375
set label "-75°" at graph 0.61,0.05
set label "-45°" at graph 0.825,0.175 
set label "8" at graph 0.25,0.9 
set label "10" at graph 0.06,0.725 
set label "12" at graph 0.0,0.5 
set label "14" at graph 0.06,0.275 
set label "16" at graph 0.25,0.09

plot  	for [i=0:1] aether_data_file index i using ($1*15):($3/1000) w points pt (i+7) lc 1 title columnhead(1),\
	for [i=0:1] aether_data_file index i using 2:($3/1000) w points pt (i+7) lc 2 title columnhead(2)

unset multiplot		
