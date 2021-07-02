set encoding utf8

set title aether_title
#set xlabel aether_x_axis
set label aether_x_axis at graph 0.84,0.54 left
set ylabel aether_y_axis offset -3,0


set key outside Left
set style line 6 dt 3 lc rgb '#bbbbbb' # gridstyle
set polar 
set angles degrees
set grid ls 6 polar 15
unset border
unset xtics
unset ytics
set rrange [0:0.03]
set view equal xy
set size ratio 1
set rtics axis 0.01,0.01 scale 0.75

set margins 2,2,2,2
set label "2" at graph 0.95,0.75 
set label "4" at graph 0.75,0.95 
set label "6" at graph 0.5,1.02 
set label "8" at graph 0.23,0.95 right 
set label "10" at graph 0.04,0.75 right
set label "12" at graph -0.02,0.5 right
set label "14" at graph 0.04,0.25 right
set label "16" at graph 0.23,0.05 right
set label "18" at graph 0.49,-0.02
set label "20" at graph 0.75,0.04
set label "22" at graph 0.95,0.25    

		
plot aether_data_file index 0 using ($1/24*360):($2*0.03) w filledcurve fc "#aa9999ee" title title_cov_w,\
aether_data_file index 1 using ($1/24*360):($2*0.03) w filledcurve fc "#aaee9999" title title_cov_s,\
aether_data_file2 index 0 using ($1/24*360):($3) w lines lc 3 title title_amp_w,\
aether_data_file2 index 1 using ($1/24*360):($3) w lines lc 7 title title_amp_s
