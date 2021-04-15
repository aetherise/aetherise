@echo off
FOR /F "tokens=3" %%a IN ('reg query "HKCU\Control Panel\International" /v LocaleName ^| find "LocaleName"') DO set locale=%%a


if "%1"=="" (	
	echo %0 ^<filename^> [title] [output-filename.svg]
	exit /B 1
)

echo Changing code page to utf-8
chcp 65001

set "title=%~2"


if "%3"=="" (
	set "term=wxt 0 font 'Arial,9' size 640,480"	
) else (		
	rem 8*72/96=6pt, 320/96*2.54=8.5cm
	set "term=svg font 'Arial,8' size 400,300 linewidth 0.5 background '#ffffff'; set pointsize 0.5; set output '%~3'"			
)

if "%locale:~0,2%" EQU "de" (
	set "xaxis=Azimut"
	set "yaxis=Abstand in 1/10 Streifen"	
	set "decimalsign=,"
) else (
	set "xaxis=Azimuth"
	set "yaxis=Distance in 1/10 fringe"
	set "decimalsign=."	
)

start /B gnuplot -persist -e "aether_data_file='%1'; aether_x_axis='%xaxis%';aether_y_axis='%yaxis%';aether_title='%title%';set decimalsign '%decimalsign%';set terminal %term%" raw.gnuplot


