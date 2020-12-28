#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# lib
source("aetherise.R")

if (length(args) == 0) {
	message("test_for_normality.R [options...] <files...>")
	message("  Test for normality for the given data sheets")
	message("  using the Shapiro-Wilk test.")
	aether.showHelp()
	quit()
}



#
# Percentage of rejected, with uncertainty (95% level).
# Returns a formatted string.
#
quota <- function(r,n)
{
	ci <- aether.confidenceInterval(r,n)
	return (paste(format((ci[1])*100,digits=1,nsmall=1),"+/-",format(ci[2]*100,digits=1,nsmall=1),"%"))
}




# New tables of coefficients and percentage points for the w test for normality.
# Parrish 1992, Table 3
table50 <- c(
	    0,     0, 0.933, 0.914, 0.915, 0.919, 0.923, 0.927, 0.932, 0.935,
	0.939, 0.942, 0.944, 0.947, 0.949, 0.951, 0.953, 0.955, 0.956, 0.958,
	0.959, 0.961, 0.962, 0.963, 0.964, 0.965, 0.966, 0.967, 0.968, 0.969,
	0.969, 0.970, 0.971, 0.971, 0.972, 0.973, 0.973, 0.974, 0.974, 0.975,
	0.975, 0.976, 0.976, 0.976, 0.977, 0.977, 0.978, 0.978, 0.798, 0.979
	)

table25 <- c(
	    0,     0, 0.854, 0.860, 0.867, 0.877, 0.885, 0.892, 0.899, 0.906,
	0.911, 0.916, 0.920, 0.924, 0.928, 0.930, 0.934, 0.936, 0.939, 0.941,
	0.943, 0.945, 0.947, 0.948, 0.950, 0.951, 0.953, 0.954, 0.955, 0.956,
	0.957, 0.958, 0.959, 0.960, 0.961, 0.962, 0.963, 0.964, 0.964, 0.965,
	0.966, 0.967, 0.967, 0.968, 0.968, 0.969, 0.969, 0.970, 0.970, 0.971
)



#
# Percentage Point for the given number of samples n.
# 50% level
# 
#
percentage_point_level_50 <- function(n)
{		
	if (n>50) w50 <- 0.979 + 0.001*(n-50)*0.3
	else w50 <- table50[n]
	return (w50)
}



#
# Percentage Point for the given number of samples n.
# 25% level
# 
#
percentage_point_level_25 <- function(n)
{		
	if (n>50) w25 <- 0.971 + 0.001*(n-50)*0.5
	else w25 <- table25[n]
	return (w25)
}


n_1 <- 0
n_5 <- 0
n_10 <- 0
n_25 <- 0
n_50 <- 0
n_x <- 0

total_samples <- 0

for (i in 1:length(aether.files)) {
	message(aether.files[[i]])
	
	datasheet <- aether.readDatasheet(aether.files[[i]])
	reduced <- aether.reduce(datasheet)
	samples <- length(reduced[[1]])
	total_samples <- total_samples + samples
	
	w50 <- percentage_point_level_50(samples)
	w25 <- percentage_point_level_25(samples)
	
	for (j in 1:aether.azimuths()) {		
		result <- shapiro.test(reduced[[j]]);
		p <- result["p.value"] # not accurate for p>0.1
		W <- result["statistic"]
		
		if (p < 0.01) {
			n_1 <- n_1+1
			n_5 <- n_5+1	
			n_10 <- n_10+1
			n_25 <- n_25+1
			n_50 <- n_50+1
		}	
		else if (p < 0.05) {
			n_5 <- n_5+1	
			n_10 <- n_10+1
			n_25 <- n_25+1
			n_50 <- n_50+1
		}
		else if (p < 0.10) {
			n_10 <- n_10+1	
			n_25 <- n_25+1
			n_50 <- n_50+1
		}
		#else if (p < 0.25) { # gives almost same result
		else if (W < w25) {
			n_25 <- n_25+1
			n_50 <- n_50+1
		}
		#else if (p < 0.50) {
		else if (W < w50) { 
			n_50 <- n_50+1
		}
		
		# n_x
		if (p < 0.20) {
			n_x <- n_x+1
		}
	}	
}

total <- length(aether.files)*aether.azimuths()


cat("Mean number of samples per azimuth:",format(total_samples/length(aether.files),digits=1,nsmall=1),"\n")
cat("Rejected\n")
#cat(" at xx% level: ",format(n_x, width=4),"/",total,"  ",format(quota(n_x, total),width=16,justify="right"),"\n")
cat(" at 50% level: ",format(n_50,width=4),"/",total,"  ",format(quota(n_50,total),width=16,justify="right"),"\n")
cat(" at 25% level: ",format(n_25,width=4),"/",total,"  ",format(quota(n_25,total),width=16,justify="right"),"\n")
cat(" at 10% level: ",format(n_10,width=4),"/",total,"  ",format(quota(n_10,total),width=16,justify="right"),"\n")
cat(" at  5% level: ",format(n_5, width=4),"/",total,"  ",format(quota(n_5, total),width=16,justify="right"),"\n")
cat(" at  1% level: ",format(n_1, width=4),"/",total,"  ",format(quota(n_1, total),width=16,justify="right"),"\n")

#warnings()

