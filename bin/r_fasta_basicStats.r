#! /usr/bin/env Rscript
c<-scan("stdin", quiet=TRUE)
cat("\n","min:",min(c),"\n",
	"max:",max(c),"\n",
	"median:",median(c),"\n",
	"mean:",mean(c),"\n\n")
