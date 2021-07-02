#
# Provides initialization and functions
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
	for (j in 1:length(ds[[1]])) {
		if ( grepl("b",ds[j,18]) && !aether.options.ignoreall
		|| grepl("C",ds[j,18]) && aether.options.ignoreall
		|| 	grepl("c",ds[j,18])	) 
		{
			ignore_lines[length(ignore_lines)+1] <- j
		}
	}
	if (length(ignore_lines)>0) 
		ds <- ds[-ignore_lines,] # remove rows
	
	# reverse sign	
	for (j in 1:length(ds[[1]])) {
		if ( grepl("i",ds[j,18]) && !aether.options.ignoreall
		|| grepl("R",ds[j,18]) && aether.options.ignoreall
		|| grepl("r",ds[j,18]) ) 
		{
			for (i in 1:17) 
				ds[j,i] <- -ds[j,i]
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
	for (j in 1:length(ds[[1]]) ) {
		delta <- ds[j,17]-ds[j,1]
		for (i in 1:17) {	
			ds[j,i] <- ds[j,i] - delta/16*(i-1)
		}
	}
	
	

	# offset
	for (j in 1:length(ds[[1]]) ) {
		mean <- 0
		for (i in 1:16) {	
			mean <- mean + ds[j,i]
		}
		mean <- mean/16

		for (i in 1:17) {	
			ds[j,i] <- ds[j,i] - mean
		}				
	}
	

	# to single period
	if (aether.options.single) {
		for (j in 1:length(ds[[1]]) ) {
			for (i in 1:8) {	
				ds[j,i] <- (ds[j,i]+ds[j,i+8])/2
			}
			ds[j,9] = (ds[j,9]+ds[j,17])/2		
			# copy
			for (i in 2:9) { 
				ds[j,i+8] <- ds[j,i]
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

