#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

source("aetherise.R")

if (length(args) == 0) {
	message("histogram.R [options...] <file>")
	message("  Show a histogram for all azimuths.")
	aether.showHelp()
	quit()
}

if (length(aether.files) != 1) {
	stop("expected one file")
}

datasheet <- aether.readDatasheet(aether.files[[1]])
reduced <- aether.reduce(datasheet)

X11()

for (i in 1:aether.azimuths()) {
	hist(reduced[[i]],main=paste("Azimuth",i,"/",aether.azimuths()),xlab="Fringe Displacement in 1/10 Fringe")
	
	result = shapiro.test(reduced[[i]])
	cat(i,"/",aether.azimuths()," Test for normality, p=",format(result[[2]]),"\n")

	message("Press return to continue")
	invisible(readLines("stdin", n=1))
}



