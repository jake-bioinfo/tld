ddply_thresh <- function(df, s.thresh, e.thresh ) {#, cores) {
  
  require(plyr)
  #require(parallel)
  #require(doParallel)
  
  # Setup parallel environment
  #noCores <- cores
  #cl <- makeCluster(noCores)
  #registerDoParallel()
  
  source('~/rscripts_tmp/R/telo_modules/fxns/findTelLength.R')
  iter.label <- paste(c(s.thresh, "-", e.thresh), collapse = "")
  iter.df <- ddply(df,.(s.name, r.name, s.win, r.type, read.median),
                  findTelLength, s.thresh, e.thresh)#, .parallel=TRUE)
  
  iter.df$norm <- "raw"
  iter.df$tel.length <- as.numeric(iter.df$tel.length)/1000
  iter.df$read.length <- as.numeric(iter.df$read.length)/1000
  iter.df <- iter.df[iter.df$tel.length>0,]
  med_1 <- unique(iter.df$read.median)[1]
  med_2 <- unique(iter.df$read.median)[2]

  if (nrow(iter.df)>0){
    iter.df$s.thresh <- s.thresh
    iter.df$e.thresh <- e.thresh
    iter.df$threshold <- iter.label
  #  print(iter.label)

    # Create normalized data frame
    iter.norm.df <- ddply(iter.df,.(s.name, r.name, s.win,
                                    r.type, read.median, read.length,
                                    trunc.start, trunc.end, telomere.end,
                                    s.thresh, e.thresh, threshold),
                          summarize, tel.length=round((tel.length/(read.median/1000))*
                                                        (med_1+med_2)/2/1000, digits = 3),
                          norm="normalized")#, .parallel=TRUE)
    
    # Combine iter.df and normalized df
    iter.df <- rbind(iter.df, iter.norm.df)
  }
  
  return(iter.df)
  #stopCluster(cl)
}
