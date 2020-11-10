#!/usr/bin/env Rscript
# Import libraries
suppressPackageStartupMessages(library("optparse"))                                
suppressPackageStartupMessages(library("stats"))

suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(ggplot2))

# Sourcing useful functions
source('/home/jaker/rscripts_tmp/R/telo_modules/fxns/modes.R')
source('/home/jaker/rscripts_tmp/R/telo_modules/fxns/ddply_thresh.R')
source('/home/jaker/rscripts_tmp/R/telo_modules/fxns/stopQuietly.R')

# Import options
# option_list <- list(make_option(c("-v", "--verbose"), action = "store_true", default = TRUE, 
#                                 help = "Print extra output [default]"),                                      
#                     make_option(c("-q", "--quietly"), action = "store_false", 
#                                 dest = "verbose", help = "Print little output"),
#                     make_option(c("-a", "--auto_par"), action = "store_true", default = TRUE, 
#                                 help = "automatically determine number of processors"),
#                     make_option(c("-i", "--in_csv"), type = "character",  
#                                 help = "full path to telomere comma separated sheet", 
#                                 metavar = "file"), 
#                     make_option(c("-o", "--out_path"), type = "character", 
#                                 help = "full path to folder where results should be stored"),
#                     make_option(c("-p", "--prefix"), type = "character", 
#                                 help = "prefix for output filenames"),
#                     make_option(c("-t", "--threads"), type = "integer", 
#                                 help = "number of processors to allocate, DEFAULT=max-2"),
#                     make_option(c("-m", "--median_values"), type = "character", 
#                                 help = "Provide comma separated list of median values, ex: \"2355,3577,8977\""),
#                     make_option(c("-r", "--rename_samples"), type = "character",
#                                 help = "do samples need to be renamed? If yes, provide list, ex: \"wt,irr\"")
# )         
# 
# # get command line options, if help option encountered print help and exit,          
# # otherwise if options not found on command line then set defaults,                  
# opt <- parse_args(OptionParser(option_list=option_list))                             
# 
# # print some progress messages to stderr if \"quietly\" wasn't requested             
# if ( opt$verbose ) {                                                                 
#   write(paste("\nStarting to process telomere reads at:", 
#               Sys.time(), collapse = ""), stderr())         
# }                                                                                    
# 
# # Determining number of threads
# if ( is.null(opt$threads) ) {
#   cat("\nNumber of parallel processors automatically set to max-2.\n")
#   noCores <- detectCores() - 2
# } else {
#   noCores <- opt$threads
# }
# 
# # Confirm all inputs
# cat("\n These are the options you submitted: \n",
#     paste("\tInput CSV:", opt$in_csv, collapse = ""), "\n", 
#     paste("\tOut path:", opt$out_path, collapse = ""), "\n",
#     paste("\tFile prefix:", opt$prefix, collapse = ""), "\n",
#     paste("\tMedian list:", opt$median_values, collapse = ""), "\n",
#     paste("\tSample list:", opt$rename_samples, collapse = ""), "\n",
#     paste("\tProcessors:", noCores, collapse = ""), "\n",
#     paste("Are these correct? (y/n)", collapse = ""), "\n"
# )
# 
# # Exit if user does not confirm variables
# opt_check <- readLines(con = "stdin", n = 1)
# if ( opt_check == "y") {
#   cat("\nOptions have been confirmed, continuing... \n")
# } else {
#   cat("\nOptions have not been confirmed, exiting. \n")
#   q(save = "no", status = 1, runLast = FALSE)
# }

# Import csv files 
# p_perTable <- opt$in_csv
# result.path <- opt$out_path
p_perTable <- "/home/jaker/pfalci/201908_telomere_lengths/data/202009_telo_homer/200sw_telomere_ranges.perc.sorted.csv"
result.path <- "/home/jaker/pfalci/201908_telomere_lengths/results/202009_homer"
noCores <- 6

df_in <- read.csv(p_perTable, comment.char = '#', stringsAsFactors = FALSE)

# Create new column with nl and assign sliding window
df_in$r.type <- "nl"
df_in$s.win <- "200"

# Convert list to one large df and rename columns and sample names, add median seq pool read lengths
df <- df_in

names <- c("s.name", "r.name", "st.range", "en.range", 
           "win.length", "t.count", "per.tel.repeat",
           "r.type", "s.win")
colnames(df) <- names

# Rename samples if needed
# if (is.character(opt$rename_samples) == TRUE) {
# sample_names <- as.list(strsplit(opt$rename_samples, ","))
#   cat(paste("\nSamples need to be renamed, looping over unique samples.", 
#             "\n*Letters and #s only.", 
#             "\n*Samples need to be in same order as median values submitted",
#             collapse = ""))
#   
  # 
  # # For loop over number of samples
  # for (sn in 1:length(unique(df$s.name))) {
  #   cat(paste("\nCurrent sample:", unique(df$s.name)[sn],
  #             "\nNew name:", sample_names[[1]][sn], collapse = ""))
  #   df$s.name[df$s.name==unique(df$s.name)[sn]] <- sample_names[[1]][sn]
  # }
  # 
  # cat("\n\nAre these correct (y/n)")
  # sm_nm_ch <- readLines(con = "stdin", n = 1)
  # if(sm_nm_ch == "y") {
  #   cat("\n Samples names are correct, continuing... ")
  # } else {
  #   cat("\nSample names are incorrect, exiting.\n")
  #   q(save = "no", status = 2, runLast = FALSE)
  # }
  # 
  # # Order data frame based on supplied sample names
  # df <- df %>% arrange(factor(s.name, levels = sample_names))

#} else {
#  cat("\nSamples are named correctly.")
#}

# Match samples and medians
#med_vec <- as.vector(strsplit(opt$median_values, ","))

# # Ensure sample # matches number of median values
# if (length(med_vec[[1]]) == length(unique(df$s.name))) {
#   cat("\n Number of median values matches number of samples, continuing... \n")
# } else {
#   cat("\n Number of median values does not match number of samples, exiting.\n")
#   q(save = "no", status = 3, runLast = FALSE)
# }
# 
# # For loop to assign median values
# for (i in 1:length(unique(df$s.name))) {
#   cat(paste("\nAssigning sample:", unique(df$s.name)[i], 
#             "\n\t Median value:", med_vec[[1]][i], 
#             "\nIs this correct (y/n)?", collapse = ""))
#   
#   # Check that median values are being assigned to correct sample
#   nm_lp_ch <- readLines(con = "stdin", n = 1)
#   if (nm_lp_ch == "y") {
#     cat("\n\t Correct median match, continuing... \n\n")
#     df$read.median[df$s.name==unique(df$s.name)[i]] <- med_vec[[1]][i]
#   } else {
#     cat("\n\t Incorrect median match, exiting.\n\n")
#     q(save = "no", status = 4, runLast = FALSE)
#   }
#   
# }

# Convert median to numeric
sample_names <- c("wt", "irr")
med_vec <- c(8306,6825)
df$s.name[df$s.name==">3d7_2s"] <- sample_names[1]
df$s.name[df$s.name==">irr_2s"] <- sample_names[2]
df$read.median[df$s.name=="wt"] <- med_vec[1]
df$read.median[df$s.name=="irr"] <- med_vec[2]
df$read.median <- as.numeric(df$read.median)

# Remove ranges from readname
df$r.name <- gsub(pattern = "seqRng.*", replacement = "", x = df$r.name)
read_count <- ddply(df,.(s.name, r.type, s.win), summarize, 
                    reads.before.Tlen=length(unique(r.name)))

# Setup parallel environment
cl <- makeCluster(noCores)
registerDoParallel()

# Downsample number of reads
wt_r_ls <- sample(as.list(unique(df$r.name[df$s.name=="wt"])), size=200)
irr_r_ls <- sample(as.list(unique(df$r.name[df$s.name=="irr"])), size=200)
df_wt <- data.frame()
df_irr <- data.frame()

df_wt <- foreach (i=1:200, .combine = rbind) %dopar% {
  df_wt_t <- df[df$r.name==wt_r_ls[[i]], ]
  df_wt_t
}

df_irr <- foreach (i=1:200, .combine = rbind) %dopar% {
  df_irr_t <- df[df$r.name==irr_r_ls[[i]], ]
  df_irr_t
}

df_s <- rbind(df_wt, df_irr)

# Run telomere length analysis at all thresholds, iot determine appropriate 
# threshold
# and return read range of non telomeric repeats
result.df <- data.frame()

for (et in seq (10, 90, 2)) {
  cat("\nProcessing end threshold:", as.character(et))
  tmp.df <- foreach (st=seq (10, 90, 2), .combine = rbind) %dopar% {
    cat("\n\tProcessing start threshold:", as.character(st))
    tmp.result <- ddply_thresh(df_s, et, st)
    tmp.result
  }
  result.df <- rbind(result.df, tmp.df)
}

# Determine appropriate threshold
st_thresh_pl <- ddply(result.df[result.df$norm=="normalized",], 
                      .(s.name, s.thresh),
                      summarize, tel.read=length(unique(r.name)), 
                      .parallel = TRUE)

st_thresh <- ggplot(st_thresh_pl, aes(x=wt, y=mpg)) + 
  geom_point()+
  geom_smooth(method=lm)

e_thresh_pl <- ddply(result.df[result.df$norm=="normalized",], 
                     .(s.name, e.thresh),
                     summarize, tel.length=mean(tel.length), 
                     .parallel = TRUE)

# Updating basic stats 
res_num_reads <- ddply(result.df[result.df$norm=="normalized",],
                       .(s.name, r.type, s.win), summarize, 
                       reads=length(unique(r.name)), .parallel=TRUE)
means.df <- ddply(result.df[result.df$norm=="normalized",],
                  .(s.name, r.type, s.win, threshold), summarize,
                  means.tel.length = mean(tel.length), .parallel=TRUE)
read_count$reads.after.Tlen <- res_num_reads$reads
read_count$thresh <- means.df$threshold
read_count$mean.tel.length <- means.df$means.tel.length

# Ampls for result.df
result.tel.stats <- ddply(result.df,.(s.name, r.type, s.win, threshold, norm),
                          summarize,
                          mode.1 = amps(tel.length)$Peaks[,1][which(amps(tel.length)$Peaks[,2] >=
                                                                      sort(amps(tel.length)$Peaks[,2], TRUE)[2])][1],
                          mode.2 = amps(tel.length)$Peaks[,1][which(amps(tel.length)$Peaks[,2] >=
                                                                      sort(amps(tel.length)$Peaks[,2], TRUE)[2])][2],
                          antimode = amps(tel.length)$Antimode[,1][which(amps(tel.length)$Antimode[,2] ==
                                                                           min(amps(tel.length)$Antimode[,2]))],
                          bimod_coeff = bimodality_coefficient(tel.length, TRUE), .parallel=TRUE)


# Create percent telomere repeats means per sample type
per.tel.stats <- ddply(df,.(s.name, r.type, s.win),
                       summarize,
                       mode.1 = amps(per.tel.repeat)$Peaks[,1][which(amps(per.tel.repeat)$Peaks[,2] >=
                                                                       sort(amps(per.tel.repeat)$Peaks[,2], TRUE)[2])][1],
                       mode.2 = amps(per.tel.repeat)$Peaks[,1][which(amps(per.tel.repeat)$Peaks[,2] >=
                                                                       sort(amps(per.tel.repeat)$Peaks[,2], TRUE)[2])][2],
                       antimode = amps(per.tel.repeat)$Antimode[,1][which(amps(per.tel.repeat)$Antimode[,2] ==
                                                                            min(amps(per.tel.repeat)$Antimode[,2]))],
                       bimod_coeff = bimodality_coefficient(per.tel.repeat, TRUE), .parallel=TRUE)

# Save dataframes
df.f <- paste(c(result.path, "/", opt$prefix, ".df.Rda"), collapse = "")
result.df.f <- paste(c(result.path, "/", opt$prefix, ".result.df.Rda"), collapse = "")
per.tel.stats.f <- paste(c(result.path, "/", opt$prefix, ".per.tel.stats.Rda"), collapse = "")
result.tel.stats.f <- paste(c(result.path, "/", opt$prefix, ".result.tel.stats.Rda"), collapse = "")
read_count.f <- paste(c(result.path, "/", opt$prefix, ".read_count.Rda"), collapse = "")

saveRDS(df, file = df.f)
saveRDS(result.df, file = result.df.f)
saveRDS(per.tel.stats, file = per.tel.stats.f)
saveRDS(result.tel.stats, file = result.tel.stats.f)
saveRDS(read_count, file = read_count.f)

#Stop Cluster
stopCluster(cl)
