@echo off
FOR /F "tokens=3" %%a IN ('reg query "HKCU\Control Panel\International" /v LocaleName ^| find "LocaleName"') DO set locale=%%a



if "%1"=="" (	
	echo %0 ^<coverage-data-filename^> ^<amplitude-data-filename^> [title] [output-filename.svg]
	exit /B 1
)

echo Changing code page to utf-8
chcp 65001

set "title=%~3"



if "%4"=="" (
	set "term=wxt 0 font 'Arial,9' size 800,600"	
) else (	
	
	rem fontsize 8pt and width 15.9cm on a 96dpi display
	rem 10.67*72/96=8pt, 640/96*2.54=16.9cm
	set "term=svg font 'Arial,10.67' size 768,320 linewidth 0.5 background '#ffffff'; set pointsize 0.7; set output '%~4'"			
)



if "%locale:~0,2%" EQU "de" (		
	rem one must quote the assignment!
	set "xaxis=Signalamplitude / \U+03bb"
	set "yaxis=Sternzeit / h"
	set "title_amp_w=KHS-Dipol, Winter"
	set "title_amp_s=KHS-Dipol, Sommer"
	set "title_cov_w=\U+00dc\berdeckung, Winter"
	set "title_cov_s=\U+00dc\berdeckung, Sommer"
	set "decimalsign=,"
) else (		
	set "xaxis=Signal Amplitude / \U+03bb"
	set "yaxis=Sidereal Time / h"
	set "title_amp_w=CMB-Dipole, Winter"
	set "title_amp_s=CMB-Dipole, Summer"
	set "title_cov_w=Coverage, Winter"
	set "title_cov_s=Coverage, Summer"
	set "decimalsign=."
)

start /B gnuplot -persist -e "aether_data_file='%1'; aether_data_file2='%2'; aether_x_axis='%xaxis%';aether_y_axis='%yaxis%';aether_title='%title%'; title_amp_w='%title_amp_w%';title_amp_s='%title_amp_s%';title_cov_s='%title_cov_s%';title_cov_w='%title_cov_w%'; set decimalsign '%decimalsign%'; set terminal %term%" coverage.gnuplot
