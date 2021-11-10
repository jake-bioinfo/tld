#! /usr/bin/env Rscript
# Get min med mean and max statistics 
c <- scan("stdin", quiet=TRUE)

# Perform calculation and print
cat("\n","min:",min(c),"\n",
	"max:",max(c),"\n",
	"median:",median(c),"\n",
	"mean:",mean(c),"\n\n")
