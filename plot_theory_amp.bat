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
	set "term=wxt 0"		
) else (
	rem  fontsize 8pt and width 15.9cm on a 96dpi display
	rem 10.67*72/96=8pt, 600/96*2.54=15.9cm	
	set "term=svg font 'Arial,8' size 320,200 linewidth 0.5 background '#ffffff'; set output '%~3'"
)

if "%locale:~0,2%" EQU "de" (
	set "xaxis=Sternzeit / h"
	set "yaxis=Signalamplitude / \U+03bb"	
	set "theorie_title=Theorie"
	set "decimalsign=,"
) else (
	set "xaxis=Sidereal Time / h"
	set "yaxis=Signal Amplitude / \U+03bb"	
	set "theorie_title=Theory"
	set "decimalsign=."
)

start /B gnuplot -persist -e "aether_data_file='%1'; aether_x_axis='%xaxis%';aether_y_axis='%yaxis%'; aether_theorie_title='%theorie_title%';aether_title='%title%'; set terminal %term%;set decimalsign '%decimalsign%'" theory_amp.gnuplot

