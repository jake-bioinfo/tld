#!/usr/bin/env Rscript
# Import libraries
suppressPackageStartupMessages(library("optparse"))                                
suppressPackageStartupMessages(library("stats"))

suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(doParallel))

# Import custom functions
source('~/pfalci/201908_telomere_lengths/github/tld/fxns/truncation.R')
source('~/pfalci/201908_telomere_lengths/github/tld/fxns/stopQuietly.R')

# Import options
option_list <- list(make_option(c("-v", "--verbose"), action = "store_true", default = TRUE, 
                                help = "Print extra output [default]"),                                      
                    make_option(c("-q", "--quietly"), action = "store_false", 
                                dest = "verbose", help = "Print little output"),
                    make_option(c("-i", "--in_fasta1"), type = "character",  
                                help = "full path to telomere fasta file", 
                                metavar = "file"),
                    make_option(c("-s", "--in_sub_fasta"), type = "character",  
                                help = "full path to telomere sub sampled fasta file", 
                                metavar = "file"),
                    make_option(c("-o", "--out_path"), type = "character", 
                                help = "full path to folder where results should be stored"),
                    make_option(c("-p", "--prefix"), type = "character",
                                help = "prefix for naming output files, must match from module 1"),
                    make_option(c("-t", "--threads"), type = "integer", 
                                help = "number of processors to allocate, DEFAULT=max-2")
)         

# get command line options, if help option encountered print help and exit,          
# otherwise if options not found on command line then set defaults,                  
opt <- parse_args(OptionParser(option_list=option_list))                             

# print some progress messages to stderr if \"quietly\" wasn't requested             
if ( opt$verbose ) {                                                                 
  write(paste("\nStarting to process telomere reads at:", 
              Sys.time(), collapse = ""), stderr())         
}                                                                                    

# Determining number of threads
if ( is.null(opt$threads) ) {
  cat("\nNumber of parallel processors automatically set to max-2.\n")
  noCores <- detectCores() - 2
} else {
  noCores <- opt$threads
}

# Confirm all inputs
cat("\n These are the options you submitted: \n",
    paste("\tInput fasta1:", opt$in_fasta1, collapse = ""), "\n",
    paste("\tInput fasta2 sub sample:", opt$in_sub_fasta, collapse = ""), "\n",
    paste("\tOut path:", opt$out_path, collapse = ""), "\n",
    paste("\tPrefix:", opt$prefix, collapse = ""), "\n",
    paste("\tProcessors:", noCores, collapse = ""), "\n",
    paste("Are these correct? (y/n)", collapse = ""), "\n"
)

# Exit if user does not confirm variables
opt_check <- readLines(con = "stdin", n = 1)
if ( opt_check == "y") {
  cat("\nOptions have been confirmed, continuing... \n\n")
} else {
  cat("\nOptions have not been confirmed, exiting. \n\n")
  q(save = "no", status = 1, runLast = FALSE)
}

# Setup parallel environment
cl <- makeCluster(noCores)
registerDoParallel()

# Paths
result.path <- opt$out_path

# Import data, truncate sequences and blast truncated sequences
# Import fasta and basic name processing
# Telomere length results, normalized data only
result.df.f <- paste(c(result.path, "/", opt$prefix, ".result.df.Rda"), collapse = "")
result.df <- readRDS(result.df.f)
result.df <- result.df[result.df$norm=="normalized", ]

# wt nl
pf3d7.f <- opt$in_fasta
pf3d7.seq <- readDNAStringSet(pf3d7.f, format = "fasta")

# irradiated nl
pfIrr.f <- opt$in_sub_fasta
pfIrr.seq <- readDNAStringSet(pfIrr.f, format = "fasta")

# Truncate sequences
res.seq.ls <- list()
seq.ls <- DNAStringSetList(pf3d7.seq, pfIrr.seq)
seq.ls <- unlist(seq.ls)
names(seq.ls) <- gsub(pattern = " RQ.*", replacement = "", names(seq.ls))
names(seq.ls) <- gsub(pattern = "/", replacement = "", names(seq.ls))
names(seq.ls) <- gsub(pattern = " .*", replacement = "", names(seq.ls))
res.seq.df <- data.frame()

res.seq.df <- ddply(result.df,
                    .(r.name, s.name, trunc.start, 
                      trunc.end, telomere.end, tel.length, 
                      read.length),
                    truncation, seq.ls, .parallel = TRUE)

dna.trunc <- DNAStringSet(res.seq.df$trunc.seq)
names(dna.trunc) <- res.seq.df$r.name

# Write all truncated reads to fasta
writeXStringSet(dna.trunc, paste(c(result.path, "/", opt$prefix, ".dna.trunc.fa.gz"), collapse = ""), 
                compress = TRUE, format = "fasta")

# Stop Cluster
stopCluster(cl)
