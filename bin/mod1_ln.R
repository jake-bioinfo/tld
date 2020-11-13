#!/usr/bin/env Rscript
## Make changes from mod1_ln_tmp as option for ont data!!
# Import libraries
suppressPackageStartupMessages(require("optparse"))                                
suppressPackageStartupMessages(require("stats"))

suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(foreach))
suppressPackageStartupMessages(require(doParallel))

# Sourcing useful functions
source('/tld/fxns/modes.R')
source('/tld/fxns/ddply_thresh.R')
source('/tld/fxns/stopQuietly.R')

# source('~/tmp/docker/test_data/tld/fxns/modes.R')
# source('~/tmp/docker/test_data/tld/fxns/ddply_thresh.R')
# source('~/tmp/docker/test_data/tld/fxns/stopQuietly.R')


# Import options
option_list <- list(make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
                                help = "Print extra output [default]"),
                    make_option(c("-q", "--quietly"), action = "store_false",
                                dest = "verbose", help = "Print little output"),
                    make_option(c("-a", "--auto_par"), action = "store_true", default = TRUE,
                                help = "automatically determine number of processors"),
                    make_option(c("-i", "--in_csv"), type = "character",
                                help = "full path to telomere comma separated sheet",
                                metavar = "file"),
                    make_option(c("-o", "--out_path"), type = "character",
                                help = "full path to folder where results should be stored"),
                    make_option(c("-p", "--prefix"), type = "character",
                                help = "prefix for output filenames"),
                    make_option(c("-f", "--platform"), type = "character",
                                help = "sequencing platform"),
                    make_option(c("-t", "--threads"), type = "integer",
                                help = "number of processors to allocate, DEFAULT=max-2"),
                    make_option(c("-m", "--median_values"), type = "character",
                                help = "Provide comma separated list of median values, ex: \"2355,3577,8977\""),
                    make_option(c("-r", "--rename_samples"), type = "character",
                                help = "do samples need to be renamed? If yes, provide list, ex: \"wt,irr\"")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

# print some progress messages to stderr if \"quietly\" wasn't requested
if ( opt$verbose ) {
  write(paste("\nStarting to process telomere reads at:",
              Sys.time(), collapse = ""), stderr())
}

# opt <- list()
# opt$in_csv <- '~/tmp/docker/test_data/tld/data/w_dir/200sw_telomere_ranges.perc.sorted.csv'
# opt$out_path <- '~/tmp/docker/test_data/tld/data/o_dir'
# opt$prefix <- 'dockTest'
# opt$median_values <- '1,1'
# opt$rename_samples <- 'Spara,Scere'
# opt$platform <- 'pb'
# opt$threads <- 7

# Determining number of threads
if ( is.null(opt$threads) ) {
  cat("\nNumber of parallel processors automatically set to max-2.\n")
  noCores <- detectCores() - 2
} else {
  noCores <- opt$threads
}

# Confirm all inputs
cat("\n These are the options you submitted: \n",
    paste("\tInput CSV:", opt$in_csv, collapse = ""), "\n",
    paste("\tOut path:", opt$out_path, collapse = ""), "\n",
    paste("\tFile prefix:", opt$prefix, collapse = ""), "\n",
    paste("\tMedian list:", opt$median_values, collapse = ""), "\n",
    paste("\tSample list:", opt$rename_samples, collapse = ""), "\n",
    paste("\tSequencing platform:", opt$platform, collapse = ""), "\n",
    paste("\tProcessors:", noCores, collapse = ""), "\n",
    paste("Are these correct? (y/n)", collapse = ""), "\n"
)

# Exit if user does not confirm variables
opt_check <- readLines(con = "stdin", n = 1)
if ( opt_check == "y") {
  cat("\nOptions have been confirmed, continuing... \n")
} else {
  cat("\nOptions have not been confirmed, exiting. \n")
  q(save = "no", status = 1, runLast = FALSE)
}

# Import csv files
p_perTable <- opt$in_csv
result.path <- opt$out_path

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
if (is.character(opt$rename_samples) == TRUE) {
  sample_names <- as.list(strsplit(opt$rename_samples, ","))
  cat(paste("\nSamples need to be renamed, looping over unique samples.",
            "\n*Letters and #s only.",
            "\n*Samples need to be in same order as median values submitted",
            collapse = ""))
  
  for (sn in 1:length(unique(df$s.name))) {
    cat(paste("\nCurrent sample:", unique(df$s.name)[sn],
              "\nNew name:", sample_names[[1]][sn], collapse = ""))
    df$s.name[df$s.name==unique(df$s.name)[sn]] <- sample_names[[1]][sn]
  }
  
  cat("\n\nAre these correct (y/n)")
  sm_nm_ch <- readLines(con = "stdin", n = 1)
  if(sm_nm_ch == "y") {
    cat("\n Samples names are correct, continuing... ")
  } else {
    cat("\nSample names are incorrect, exiting.\n")
    q(save = "no", status = 2, runLast = FALSE)
  }
  
  # Order data frame based on supplied sample names
  df <- df %>% arrange(factor(s.name, levels = sample_names))
  
} else {
  cat("\nSamples are named correctly.")
}

# Match samples and medians
med_vec <- as.vector(strsplit(opt$median_values, ","))

# Ensure sample # matches number of median values
if (length(med_vec[[1]]) == length(unique(df$s.name))) {
  cat("\n Number of median values matches number of samples, continuing... \n")
} else {
  cat("\n Number of median values does not match number of samples, exiting.\n")
  q(save = "no", status = 3, runLast = FALSE)
}

# For loop to assign median values
for (i in 1:length(unique(df$s.name))) {
  cat(paste("\nAssigning sample:", unique(df$s.name)[i],
            "\n\t Median value:", med_vec[[1]][i],
            "\nIs this correct (y/n)?", collapse = ""))
 #df$read.median[df$s.name==unique(df$s.name)[i]] <- med_vec[[1]][i]
  
  # Check that median values are being assigned to correct sample
  nm_lp_ch <- readLines(con = "stdin", n = 1)
  if (nm_lp_ch == "y") {
    cat("\n\t Correct median match, continuing... \n\n")
    df$read.median[df$s.name==unique(df$s.name)[i]] <- med_vec[[1]][i]
  } else {
    cat("\n\t Incorrect median match, exiting.\n\n")
    q(save = "no", status = 4, runLast = FALSE)
  }
}

# Convert median to numeric
df$read.median <- as.numeric(df$read.median)

# Remove ranges from readname
df$r.name <- gsub(pattern = "seqRng.*", replacement = "", x = df$r.name)
read_count <- ddply(df,.(s.name, r.type, s.win), summarize, 
                    reads.before.Tlen=length(unique(r.name)))

# Wrap theshold determiniation in if statement if platform is pacbio
if (opt$platform == "pb") { 
  
  # Setup parallel environment
  cl <- makeCluster(noCores)
  registerDoParallel()
  
  # Downsample number of reads
  sample_size <- round((0.15 * (length(unique(df$r.name)))), digits = -1)
  df_sr <- sample(unique(df$r.name), sample_size)
  df_s <- data.frame()
  
  function_g <- function(df, int, df_sr) {
    df_s_tmp <- df[grep(df_sr[int], df$r.name), ]
    return(df_s_tmp)
  }
  
  cat("\nStarting subsampling for all threshold processing.\n")
  df_s <- foreach(s=1:length(df_sr), .combine = rbind) %dopar% {
    df_s_tmp <- function_g(df, s, df_sr)
    df_s_tmp
  }
  cat("\nSubsampling completed.\n")
  
  # Create data frame from 1000 downsampled reads to determine threshold to use 
  # going forward
  # Create different test data frames for ont and pb bc ont data has messier 
  # telomeres and often does not meet the criteria at lower and higher end 
  # thresholds
  et.s <- 5
  st.s <- 5
  et.e <- 85
  st.e <- 85
  
  # Add total for progress bar 
  total <- et.e
  
  # create progress bar
  pb <- txtProgressBar(min = et.s, max = total, styl = 3)
  
  cat("\nStarting to process all thresholds for threshold determination.\n")
  s.result.df <- data.frame()
  for (et in et.s:et.e) {
    setTxtProgressBar(pb, et)  
    s.result.df.tmp <- foreach(st=st.s:st.e, .combine=rbind) %dopar% {
      tmp.result <- ddply_thresh(df_s, et, st)
      tmp.result
    }
    s.result.df <- rbind(s.result.df.tmp, s.result.df)
  }
  close(pb)
  cat("\nAll threshold processing completed.\n")
  
  # Producing end threshold data frame to determine end threshold
  et.df <- ddply(s.result.df[s.result.df$norm=="normalized", ],
                 .(s.name, e.thresh, s.thresh), summarize,
                 tel.ln.mean=mean(tel.length), .parallel = TRUE)
  
  # Subset in order to determine asymptote, convert data to numeric
  ets.df <- et.df[et.df$s.thresh==50&et.df$s.name==unique(s.result.df$s.name)[1], ]
  ets.df$e.thresh <- as.numeric(ets.df$e.thresh)
  
  # Producing start threshold data frame to determine start threshold
  st.df <- ddply(s.result.df[s.result.df$norm=="normalized", ],
                 .(s.name, s.thresh, e.thresh), summarize,
                 r.per.st=length(unique(r.name)), .parallel = TRUE)
  
  # Subset start threshold data frame in order to determine max
  sts.df <- st.df[st.df$e.thresh==50&st.df$s.name==unique(s.result.df$s.name)[1], ]
  sts.df$r.per.st <- as.numeric(sts.df$r.per.st)
  
  # Determine asymptote of end threshold and max of start threshold
  # subsets
  # If reads are ccs reads then find vertical asymptote values
  x <- grepl("ccs", head(df$r.name, 1))
  if ( x ) {
    et.res.vec <- vector()
    for (et in et.s:et.e) {
      #cat("\nThis is step:", et)
      
      if(et==et.s){
        #cat("\n\tThis is first step continuing without calculation")
      } else {
        
        p.et <- et -1
        #cat("\n\tThis is previous end threshold:", p.et)
        p.tl <- ets.df$tel.ln.mean[ets.df$e.thresh==p.et]
        #cat("\n\tThis is previous telomere length:", p.tl)
        
        t.tl <- ets.df$tel.ln.mean[ets.df$e.thresh==et]
        #cat("\n\tThis is tmp.tl:", t.tl)
        
        diff <- p.tl - t.tl
        #cat("\n\tThis is difference:", diff)
        
        per.diff <- diff/median(ets.df$tel.ln.mean)
        #cat("\n\tThis is percent difference of current length:", per.diff)
        
        if( per.diff > 0.006 ) {
          #cat("\n\tPercent difference is greather than 0.6%, adding to result vector.")
          et.res.vec <- c(et.res.vec, et)
        } 
      }
    }
    # For loop for sts.df for ccs reads to determine vertical asymptote
    st.res.vec <- vector()
    for (st in st.s:st.e) {
      #cat("\nThis is step:", st)
      
      if(st==st.s){
        #cat("\n\tThis is first step continuing without calculation")
      } else {
        
        p.st <- st - 1
        #cat("\n\tThis is previous start threshold:", p.st)
        pr.n <- sts.df$r.per.st[sts.df$s.thresh==p.st]
        #cat("\n\tThis is previous read number:", pr.n)
        
        t.rn <- sts.df$r.per.st[sts.df$s.thresh==st]
        #cat("\n\tThis is temporary read number:", t.rn)
        
        diff <- pr.n - t.rn
        #cat("\n\tThis is difference:", diff)
        
        per.diff <- diff/median(sts.df$r.per.st)
        #cat("\n\tThis is percent difference of current length:", per.diff)
        
        if( per.diff > 0.006 ) {
          #cat("\n\tPercent difference is greater than 0.6%, adding to result vector.")
          st.res.vec <- c(st.res.vec, st)
        }
        
      }
    }
    
  } else {
  
  # For loop for ets.df to determine horizontal asymptote
  et.res.vec <- vector()
  for (et in et.s:et.e) {
    #cat("\nThis is step:", et)
    
    if(et==et.s){
      #cat("\n\tThis is first step continuing without calculation")
    } else {
      
      p.et <- et -1
      #cat("\n\tThis is previous end threshold:", p.et)
      p.tl <- ets.df$tel.ln.mean[ets.df$e.thresh==p.et]
      #cat("\n\tThis is previous telomere length:", p.tl)
      
      t.tl <- ets.df$tel.ln.mean[ets.df$e.thresh==et]
      #cat("\n\tThis is tmp.tl:", t.tl)
      
      diff <- p.tl - t.tl
      #cat("\n\tThis is difference:", diff)
      
      per.diff <- diff/median(ets.df$tel.ln.mean)
      #cat("\n\tThis is percent difference of current length:", per.diff)
        
      if( per.diff > 0.002 | per.diff <= 0 ) {
        #cat("\n\tPercent difference is greather than 0.2% or negative, continuing.")
      } else {
        #cat("\n\tPercent difference is less than 0.2%, adding to vec")
        et.res.vec <- c(et.res.vec, et)
      }
      
    }
  } 
  
  # For loop for sts.df to determine horizontal asymptote
  st.res.vec <- vector()
  for (st in st.s:st.e) {
    #cat("\nThis is step:", st)
    
    if(st==st.s){
      #cat("\n\tThis is first step continuing without calculation")
    } else {
      
      p.st <- st - 1
      #cat("\n\tThis is previous start threshold:", p.st)
      pr.n <- sts.df$r.per.st[sts.df$s.thresh==p.st]
      #cat("\n\tThis is previous read number:", pr.n)
      
      t.rn <- sts.df$r.per.st[sts.df$s.thresh==st]
      #cat("\n\tThis is temporary read number:", t.rn)
      
      diff <- pr.n - t.rn
      #cat("\n\tThis is difference:", -diff)
      
      per.diff <- -diff/median(sts.df$r.per.st)
      #cat("\n\tThis is percent difference of current length:", per.diff)
      
      if(per.diff > 0.006 | per.diff < 0 ) {
        #cat("\n\tPercent difference is greather than 0.6% or negative, continuing.")
      } else {
        #cat("\n\tPercent difference is less than 0.6%, adding to vec")
        st.res.vec <- c(st.res.vec, st)
      }
      
    }
  }
  
  }
  
  # Determine threshold values based on median values of result vector
  end.threshold <- round(quantile(et.res.vec)[2])
  if (is.na(end.threshold)) {
    end.threshold <- 25
  }
  
  # Determine threshold values based on median values of result vector
  start.threshold <- round(quantile(st.res.vec)[2])
  if (is.na(start.threshold)) {
    start.threshold <- 60
  }
  
  } else {
  # These are ont reads
    end.threshold <- 35
    start.threshold <- 40
}

# Determine telomere lengths based on thresholds
result.df <- ddply_thresh(df, start.threshold, end.threshold)

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
                                                                      sort(amps(tel.length)$Peaks[,2], TRUE)[1])][1], 
                          .parallel = TRUE)

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
s.result.df.f <- paste(c(result.path, "/", opt$prefix, ".s.result.df.Rda"), collapse = "")

saveRDS(df, file = df.f)
saveRDS(result.df, file = result.df.f)
saveRDS(per.tel.stats, file = per.tel.stats.f)
saveRDS(result.tel.stats, file = result.tel.stats.f)
saveRDS(read_count, file = read_count.f)
saveRDS(s.result.df, file = s.result.df.f)


#Stop Cluster
stopCluster(cl)