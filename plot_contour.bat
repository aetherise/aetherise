@echo off
FOR /F "tokens=3" %%a IN ('reg query "HKCU\Control Panel\International" /v LocaleName ^| find "LocaleName"') DO set locale=%%a


if "%1"=="" (	
	echo %0 ^<filename^> [title] [output-filename.svg]
	exit /B 1
)

echo Changing code page to utf-8
chcp 65001

set "title=%~2"


set asize=0.78
set origin_x=0.115
set origin_y=0.121

if "%3"=="" (
	set "term=wxt 0 font 'Arial,9' size 1200,600"	
) else (	
	rem gnuplot 5.3	
	set asize=0.765
	set origin_x=0.127
	set origin_y=0.138
	rem fontsize 8pt and width 15.9cm on a 96dpi display
	rem 10.67*72/96=8pt, 640/96*2.54=16.9cm
	set "term=svg font 'Arial,10.67' size 640,420 linewidth 0.5 background '#ffffff'; set pointsize 0.5; set output '%~3'"			
)



if "%locale:~0,2%" EQU "de" (	
	echo WARNUNG: Multiplots müssen von Hand ausgerichtet werden.
	rem one must quote the assignment!
	set "xaxis=v / (km/s)"
	set "yaxis=Rektaszension \U+03b1 / h"
	set "y2axis=Deklination \U+03b4"
	set "decimalsign=,"
) else (
	echo WARNING: Multiplots have to be aligned manually.	
	set "xaxis=v / (km/s)"
	rem set "yaxis=Right Ascension α / h" test in utf8 is not correctly displayed
	rem set "y2axis=Declination δ"
	set "yaxis=Right Ascension \U+03b1 / h"
	set "y2axis=Declination \U+03b4"
	set "decimalsign=."
)

start /B gnuplot -persist -e "aether_size=%asize%;aether_origin_x=%origin_x%;aether_origin_y=%origin_y%;aether_data_file='%1';aether_x_axis='%xaxis%';aether_y_axis='%yaxis%';aether_y2_axis='%y2axis%';aether_title='%title%';set terminal %term%;set decimalsign '%decimalsign%'" contour.gnuplot
