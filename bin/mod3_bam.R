#!/usr/bin/env Rscript
# Import libraries
suppressPackageStartupMessages(require("optparse"))                                
suppressPackageStartupMessages(require("stats"))

suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(GenomicAlignments))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(doParallel))
suppressPackageStartupMessages(require(seqinr))

# Sourcing useful functions
source('/tld/fxns/end_bam_chr.R')
source('/tld/fxns/stopQuietly.R')

# Import options
option_list <- list(make_option(c("-v", "--verbose"), action = "store_true", default = TRUE, 
                                help = "Print extra output [default]"),                                      
                    make_option(c("-q", "--quietly"), action = "store_false", 
                                dest = "verbose", help = "Print little output"),
                    make_option(c("-i", "--in_bam"), type = "character",  
                                help = "full path to telomere bam file", 
                                metavar = "file"),
                    make_option(c("-r", "--reference"), type = "character",
                                help = "full path to reference used for alignment",
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
  write(paste("\nStarting to assign telomere ends based on bam at:", 
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
    paste("\tInput bam:", opt$in_bam, collapse = ""), "\n",
    paste("\tReference:", opt$reference, collapse = ""), "\n",
    paste("\tOut path:", opt$out_path, collapse = ""), "\n",
    paste("\tPrefix:", opt$prefix, collapse = ""), "\n",
    paste("\tProcessors:", noCores, collapse = ""), "\n",
    paste("Are these correct? (y/n)", collapse = ""), "\n"
)

# Setup parallel environment
cl <- makeCluster(noCores)
registerDoParallel()

# Import reference
ref.path <- opt$reference
ref <- read.fasta(ref.path, seqtype = "DNA")

# Import csv files
bam.path <- opt$in_bam
param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE,
                                         isSecondaryAlignment = FALSE))
bam <- readGAlignments(bam.path, index = bam.path, use.names = TRUE, param = param)
result.path <- opt$out_path

# Import result.df data frame
result.df.f <- paste(c(result.path, "/", opt$prefix, ".result.df.Rda"), collapse = "")
result.df <- readRDS(result.df.f)
result.df <- result.df[result.df$norm=="normalized", ]

# Convert bam to bam dataframe
bam.df <- as.data.frame(bam)

bam.df$rnames <- names(bam)

bam.df <- data.frame(bam.df$rnames, bam.df$seqnames,
                     bam.df$strand, bam.df$qwidth,
                     bam.df$start, bam.df$end,
                     bam.df$width)
colnames(bam.df) <- c("r.name", "seqnames", "strand", "qwidth",
                      "start", "end", "width")

# Add endedness based on bam alignment
# Combine data frame with https://rdrr.io/cran/plyr/man/join.html
result.bam.df <- merge(result.df[result.df$norm=="normalized", ], bam.df, by = "r.name")

# Only take longest read alignment from bam
r.name.ls <- list(unique(result.bam.df$r.name))
result.comb.bam.df <- data.frame()

for (r in r.name.ls[[1]]){
  tmp.df <- result.bam.df[result.bam.df$r.name==r, ]
  tmp.df <- tmp.df[order(-tmp.df$width), ][1,]
  result.comb.bam.df <- rbind(result.comb.bam.df, tmp.df)
  tmp.df <- data.frame()
}

# Update basic stats
read_count.Abam <- ddply(result.comb.bam.df, .(s.name, r.type, norm),
                         summarize, 
                         reads.after.bam = length(unique(r.name)), 
                         mean.tel.length.Abam = mean(tel.length), 
                         .parallel = TRUE)

# Create chr list
chr_nm_ls <- names(ref)
chr_ln_ls <- list()
i <- 1

while(i < length(chr_nm_ls)) {
  chr_ln_ls[paste(chr_nm_ls[i])] <- length(ref[[i]])
  i <- i + 1
}

# Add endedness
tmp.bam.df <- data.frame()

tmp.bam.df <- ddply(result.comb.bam.df,.(s.name, r.name, r.type,
                                         norm, threshold, seqnames, 
                                         strand, start, end, width, 
                                         tel.length), 
                    end_bam_chr, .parallel = TRUE)

end.bam.df <- na.omit(tmp.bam.df)

# Reformat chromosome names for further processing
end.bam.df$seqnames <- gsub("chr_", "", end.bam.df$seqnames)

# Add to basic stats after end assignment
read_count.Abam <- ddply(end.bam.df,.(s.name, r.type, norm), summarize, 
                         reads.after.en = length(unique(r.name)),
                         mean.tel.length.Aen = mean(tel.length), .parallel = TRUE)

# Save Read count, bam, telomere csv
read_count.Abam.f <- paste(c(result.path, "/", opt$prefix, ".read_count.Abam.Rda"), collapse = "")
end.bam.f <- paste(c(result.path, "/", opt$prefix, ".end.bam.csv"), collapse = "")

saveRDS(read_count.Abam, file = read_count.Abam.f)
write.csv(end.bam.df, file = end.bam.f)

# Close parallel environment
stopCluster(cl)