@echo off
FOR /F "tokens=3" %%a IN ('reg query "HKCU\Control Panel\International" /v LocaleName ^| find "LocaleName"') DO set locale=%%a

echo Changing code page to utf-8
chcp 65001

set "title=%~2"


if "%3"=="" (
	set "term=wxt 0 font 'Arial,9' size 640,480"	
) else (		
	rem 8*72/96=6pt, 320/96*2.54=8.5cm
	set "term=svg font 'Arial,10.67' size 640,480 linewidth 0.5 background '#ffffff'; set pointsize 0.7; set output '%~3'"			
)

if "%locale:~0,2%" EQU "de" (
	set "xaxis=Sternzeit / h"
	set "yaxis=Azimut"
	set "zaxis=Verschiebung / \U+03bb"
	set "data_title=Daten"
	set "theorie_title=Theorie"
	set "decimalsign=,"
) else (
	set "xaxis=Sidereal Time / h"
	set "yaxis=Azimuth"
	set "zaxis=Displacement / \U+03bb"
	set "data_title=Data"
	set "theorie_title=Theory"
	set "decimalsign=."
)

start /B gnuplot -e "aether_data_file='%1'; aether_x_axis='%xaxis%';aether_y_axis='%yaxis%';aether_z_axis='%zaxis%'; aether_data_title='%data_title%';aether_theorie_title='%theorie_title%';aether_title='%title%'; set decimalsign '%decimalsign%';set terminal %term%" signal3d.gnuplot -
