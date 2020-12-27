#
# Provides initialisation and functions
#
# include with: source("aetherise.R")
#


aether.args = commandArgs(trailingOnly=TRUE)

aether.files <- list()

# options
aether.options.single = FALSE
aether.options.ignoreall = FALSE


#
# Parse command line arguments
#
if (length(aether.args) > 0) {
	for (i in 1:length(aether.args)) {
		if (startsWith(aether.args[[i]],"-")) {
			name = aether.args[i]
			if (name == "-single") {
				aether.options.single = TRUE
			}
			else if (name == "-ignore_all") {
				aether.options.ignoreall = TRUE
			}
			else {
				stop("Unknown option: ",name,call.=FALSE)
			}
		
		}
		else {
			aether.files[length(aether.files)+1] <- aether.args[[i]]
		}
	}
}



#
#
#
aether.showHelp <- function()
{
	message("Options:")
	message("  -single        reduce double period to single period")
	message("  -ignore_all    ignore all new row attributes")
}



#
# Read a data sheet (CSV file) into a table.
# Also applies row attributes.
#
aether.readDatasheet <- function(filename) 
{
	ds <- read.table(file=filename, skip=4, sep=";")
	
	# collect line numbers to remove
	ignore_lines <- vector()
	for (i in 1:length(ds[[1]])) {
		if ( grepl("b",ds[i,18]) && !aether.options.ignoreall
		|| grepl("C",ds[i,18]) && aether.options.ignoreall
		|| 	grepl("c",ds[i,18])	) 
		{
			ignore_lines[length(ignore_lines)+1] <- i
		}
	}
	if (length(ignore_lines)>0) 
		ds <- ds[-ignore_lines,] # remove rows
	
	# reverse sign	
	for (i in 1:length(ds[[1]])) {
		if ( grepl("i",ds[i,18]) && !aether.options.ignoreall
		|| grepl("R",ds[i,18]) && aether.options.ignoreall
		|| grepl("r",ds[i,18]) ) 
		{
			for (j in 1:17) 
				ds[i,j] <- -ds[i,j]
		}
	}
	return (ds)
}



#
# Number of Azimuths
#
aether.azimuths <- function() {
	if (aether.options.single) {
		n <- 8
	}
	else {
		n <- 16
	}
	return (n)
}



#
# Reduce Data
#
aether.reduce <- function(ds) {

	# drift
	for (i in 1:length(ds[[1]]) ) {
		delta <- ds[i,17]-ds[i,1]
		for (j in 1:17) {	
			ds[i,j] <- ds[i,j] - delta/16*(j-1)
		}
	}
	
	

	# offset
	for (i in 1:length(ds[[1]]) ) {
		mean <- 0
		for (j in 1:16) {	
			mean <- mean + ds[i,j]
		}
		mean <- mean/16

		for (j in 1:17) {	
			ds[i,j] <- ds[i,j] - mean
		}				
	}
	

	# to single period
	if (aether.options.single) {
		for (i in 1:length(ds[[1]]) ) {
			for (j in 1:8) {	
				ds[i,j] <- (ds[i,j]+ds[i,j+8])/2
			}
			ds[i,9] = (ds[i,9]+ds[i,17])/2		
			# copy
			for (j in 2:9) { 
				ds[i,j+8] <- ds[i,j]
			}
		}
	}
	
			
	return (ds)
}



#
# Agresti-Coull interval for binomial distributions
#
aether.confidenceInterval <- function(k,n)
{
	c <- 1.96 # 95% confidence
	
	k_ <- k+c*c/2
	n_ <- n+c*c
	p_ <- k_/n_
	u_ <- c*sqrt(p_*(1-p_)/n_)
	return (c(p_,u_))
}



# Wait for window close
#while (!is.null(dev.list())) Sys.sleep(1)

