#!/usr/bin/env Rscript
# Import libraries
suppressPackageStartupMessages(require("optparse"))                                
suppressPackageStartupMessages(require("stats"))

suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(require(doParallel))

# Import custom functions
source('tld/fxns/truncation.R')
source('tld/fxns/stopQuietly.R')

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
  write(paste("\nStarting to truncate telomere reads for alignment at:", 
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
    paste("\tProcessors:", noCores, collapse = ""), "\n"
)

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
s1.f <- opt$in_fasta
s1.seq <- readDNAStringSet(s1.f, format = "fasta")

# irradiated nl
s2.f <- opt$in_sub_fasta
s2.seq <- readDNAStringSet(s2.f, format = "fasta")

# Truncate sequences
# First remove some characters from sequence header
res.seq.ls <- list()
seq.ls <- DNAStringSetList(s1.seq, s2.seq)
seq.ls <- unlist(seq.ls)
names(seq.ls) <- gsub(pattern = " RQ.*", replacement = "", names(seq.ls))
names(seq.ls) <- gsub(pattern = "/", replacement = "", names(seq.ls))
names(seq.ls) <- gsub(pattern = " .*", replacement = "", names(seq.ls))
res.seq.df <- data.frame()

# Truncate sequences
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