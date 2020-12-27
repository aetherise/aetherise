@echo off
FOR /F "tokens=3" %%a IN ('reg query "HKCU\Control Panel\International" /v LocaleName ^| find "LocaleName"') DO set locale=%%a


echo Changing code page to utf-8
chcp 65001


if "%1"=="" (	
	echo %0 ^<filename^> [title] [output-filename.svg]
	exit /B 1
)

set "title=%~2"


set asize=0.781
set origin_x=0.116
set origin_y=0.122

if "%3"=="" (
	set "term=wxt 0 font 'Arial,9' size 1200,600"	
) else (
	set asize=0.746	
	set	origin_x=0.139
	set	origin_y=0.155
	rem fontsize 8pt and width 15.9cm on a 96dpi display
	rem 10.67*72/96=8pt, 768/96*2.54=20.3cm
	set "term=svg font 'Arial,10.67' size 768,320 linewidth 0.5 background '#ffffff'; set pointsize 0.5; set output '%~3'"		
)



if "%locale:~0,2%" EQU "de" (	
	echo WARNUNG: Multiplots müssen von Hand ausgerichtet werden.
	rem one must quote the assignment!	
	set "xaxis=Rektaszension \U+03b1 / h"
	set "yaxis=Deklination \U+03b4"	
	set "decimalsign=,"
) else (	
	echo WARNING: Multiplots have to be aligned manually.
	set "xaxis=Right Ascension \U+03b1 / h"
	set "yaxis=Declination \U+03b4"	
	set "decimalsign=."
)

start /B gnuplot -persist -e "aether_size_x=%asize%;aether_size_y=%asize%;aether_origin_x=%origin_x%;aether_origin_y=%origin_y%;aether_data_file='%1';aether_x_axis='%xaxis%';aether_y_axis='%yaxis%';aether_title='%title%';set terminal %term%;set decimalsign '%decimalsign%'" contour_apex.gnuplot

