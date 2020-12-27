@echo off
FOR /F "tokens=3" %%a IN ('reg query "HKCU\Control Panel\International" /v LocaleName ^| find "LocaleName"') DO set locale=%%a

set "title=%~2"

if "%locale:~0,2%" EQU "de" (
	set "xaxis=Azimut"
	set "yaxis=Abstand in 1/10 Streifen"	
) else (
	set "xaxis=Azimuth"
	set "yaxis=Distance in 1/10 fringe"	
)

start /B gnuplot -persist -e "aether_data_file='%1'; aether_x_axis='%xaxis%';aether_y_axis='%yaxis%';aether_title='%title%'" raw.gnuplot


