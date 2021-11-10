ddply_thresh <- function(df, s.thresh, e.thresh ) {
  
  # Imports
  require(plyr)
  
  # Source findtelomere function
  source('/tld/fxns/findTelLength.R')
  
  # Find telomere length by read
  iter.label <- paste(c(s.thresh, "-", e.thresh), collapse = "")
  iter.df <- ddply(df,.(s.name, r.name, r.type, read.median),
                  findTelLength, s.thresh, e.thresh)
  
  iter.df$norm <- "raw"
  iter.df$tel.length <- as.numeric(iter.df$tel.length)/1000
  iter.df$read.length <- as.numeric(iter.df$read.length)/1000
  iter.df <- iter.df[iter.df$tel.length>0,]
  med_1 <- unique(iter.df$read.median)[1]
  med_2 <- unique(iter.df$read.median)[2]

  # If median is the same, assign both median one
  if (is.na(med_2)) {
    med_2 <- med_1
  }
  
  # Check that first calculation returned a result 
  if (nrow(iter.df) > 0){
    iter.df$s.thresh <- s.thresh
    iter.df$e.thresh <- e.thresh
    iter.df$threshold <- iter.label

    # Create normalized data frame for telomere length by read
    iter.norm.df <- ddply(iter.df,.(s.name, r.name, r.type, read.median, 
                                    read.length, trunc.start, trunc.end, 
                                    telomere.end, s.thresh, e.thresh, 
                                    threshold),
                          summarize, tel.length=round((tel.length/
                                                         (read.median/1000))*
                                                        (med_1+med_2)/2/1000, 
                                                      digits = 5),
                          norm="normalized")
    
    # Combine iter.df and normalized df
    iter.df <- rbind(iter.df, iter.norm.df)
  }
  
  return(iter.df)
}
