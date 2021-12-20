findTelLength <- function(df, st.thresh, en.thresh){
  
  require(plyr)
  
  # Create initial values 
  start.tel <- df$per.tel.repeat[df$st.range==0] > st.thresh
  end.tel <- df$per.tel.repeat[df$en.range==max(df$en.range)] > st.thresh
  read.length <- max(df$en.range)
  
  # Initialized tel length to 0
  tel.length <- 0
  
  # First condidtion, telomere starts at the beginning of the read
  if(start.tel & !end.tel){
    tel.end="L"
    
    # Assign first row of subset of read df that is below the end threshhold (telomere end window)
    window <- df[which(df$per.tel.repeat < en.thresh)[1],]
    
    # Calculate amount of final window which is composed of a telomere
    tel.length <- window$st.range + window$t.count
    
    # Assign value of 0, if entire read is telomere repeats
    if(is.na(window[1]) || length(window$t.count)==0) {
      tel.length <- 0
      tel.end <- NA
    }
    
    if (length(tel.length) == 0 || tel.length == 0) {
	    tel.length <- 0
	    tel.end <- NA
    }

    # Return telomere value
    trunc.s <- tel.length + 1
    trunc.e <- read.length
    
    return(c(tel.length=round(tel.length), read.length=read.length,
               trunc.start=trunc.s, trunc.end=trunc.e,
               telomere.end=tel.end))

   }

  # Second condition, telomere starts at the end of the read
  if(end.tel & !start.tel){
    
    tel.end="R"
    window <- df[which(df$per.tel.repeat < en.thresh)[length(which(df$per.tel.repeat < en.thresh))],]
    tel.length <- read.length - window$en.range + window$t.count
    
    # Assign value of 0, if entire read is telomere repeat
    if(is.na(window[1]) || length(window$t.count)==0) {
      tel.length <- 0
    }
    
    if (length(tel.length) == 0 || tel.length == 0) {
	    tel.length <- 0
	    tel.end <- NA
    }

    # Return telomere value
    trunc.s <- 1
    trunc.e <- read.length - tel.length
    
    return(c(tel.length=round(tel.length), read.length=read.length, 
               trunc.start=trunc.s, trunc.end=trunc.e, 
               telomere.end=tel.end))
    
  }
  
  # Return 0 no matter what
  trunc.s <- tel.length + 1
  trunc.e <- read.length
  tel.end <- NA

  return(c(tel.length=round(tel.length), read.length=read.length,      
               trunc.start=trunc.s, trunc.end=trunc.e,
               telomere.end=tel.end))

}
