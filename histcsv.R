#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
	message("histcsv.R <file> [column] [xdesc] [ydesc] [output-file.pdf]")
	message("  Show a histogram for the values in the given column")
	quit()
}


# parse args
col <- 1;
if (length(args) > 1) {
	col <- as.integer(args[2])
}

xdesc <- "Sample"
if (length(args) > 2) {
	xdesc <- args[3]
}

ydesc <- "Frequency"
if (length(args) > 3) {
	ydesc <- args[4]
}

if (length(args) > 4) {
	outfile <- args[5]
	pdf(file=outfile, width=11/2.54, height=11/2.54, pointsize=12);
}



# read csv and display/render histogramm

rs <- read.table(file=args[1]);

if (length(args) <= 4) {
	X11();
}
hist(rs[[col]],main="",xlab=xdesc,ylab=ydesc,col="#eeeeee");
shapiro.test(rs[[col]]);

if (length(args) <= 4) {
	message("Press return to continue")
	invisible(readLines("stdin", n=1))
}

if (length(args) > 4) {
	dev.off();
}

