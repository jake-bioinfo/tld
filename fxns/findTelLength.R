findTelLength <- function(df, st.thresh, en.thresh){
  
  require(plyr)
  # Create initial values 
  start.tel <- df$per.tel.repeat[df$st.range==0] > st.thresh
  end.tel <- df$per.tel.repeat[df$en.range==max(df$en.range)] > st.thresh
  read.length <- max(df$en.range)
  
  # print(head(df))
  # cat("\nHello",
  #     "\nThis is read:", unique(df$r.name)) #,
  #     "\n\tThis is start and end thresh:", st.thresh, "-", en.thresh,
  #     "\n\tThis is start.tel:", start.tel,
  #     "\n\tThis is end.tel:", end.tel,
  #     "\n\tThis is read.length:", read.length,
  #     "\n\tThis is start and end threshold:", st.thresh, "-", en.thresh,
  #     "\n\tThis is s.name, s.win, and r.type:", unique(df$s.name), "-", unique(df$s.win), "-", unique(df$r.type))

  # Initialized tel length to null
  tel.length <- NULL
  
  # # Null condition
  # if(!start.tel & !end.tel) {
  #   # cat("\n\nNeither start or end of the read meet the telomere start criteria, assigning zero and exiting")
  # }
  # 
  # # First condition, both ends meet condition for telomere start, return and exit
  # if(start.tel & end.tel){
  #   # cat("\n\nStart and end meet telomere start threshold, setting telomere length to 0 and exiting.")
  #   exists()
  # }
  
  # Second condition, telomere starts at the beginning of the read
  if(start.tel & !end.tel){
    tel.end="L"
    # cat("\n\nStart or read meets criteria for telomere start.")
    
    # Assign first row of subset of read df that is below the end threshhold (telomere end window)
    window <- df[which(df$per.tel.repeat < en.thresh)[1],]
    # cat("\n\tThis is window:\n")
    # print(window)
    
    # Calculate amount of final window which is composed of a telomere
    tel.length <- window$st.range + window$t.count
    
    # Assign 0 if value is na?? WHY WOULD IT BE NA? or LENGTH 0 or window is empty
    if(is.na(window[1]) || length(window$t.count)==0) {
      if(length(window$t.count) > 0) {
        if(is.na(window[1])) {
          cat("\n\nTHIS IS NA:", paste(df[1,]))
          tel.length <- NULL
        }
      }
    }
    
    # Return telomere value
    # cat("\n\tThis is telomere and read length:", round(tel.length), "-", read.length)
    trunc.s <- tel.length + 1
    trunc.e <- read.length
    # cat("\n\tThis is the range of the read which does not contain telomeric repteats:", trunc.s, "-", trunc.e, "\n\n")
    if (length(tel.length) == 0 || tel.length == 0) {tel.length = NULL}
    if (!is.null(tel.length)) {
      return(c(tel.length=round(tel.length), read.length=read.length, 
             trunc.start=trunc.s, trunc.end=trunc.e, 
             telomere.end=tel.end))
    }
  }
  
  # Third condition, telomere starts at the end of the read
  if(end.tel & !start.tel){
    
    tel.end="R"
    # cat("\n\nEnd of read meets criteria for telomere start.")
    window <- df[which(df$per.tel.repeat < en.thresh)[length(which(df$per.tel.repeat < en.thresh))],]
    # cat("\nThis is window:\n")
    # print(window)
    tel.length <- read.length - window$en.range + window$t.count
    
    # Assign 0 if value is na
    if(is.na(window[1]) || length(window$t.count)==0) {
      if(length(window$t.count) > 0) {
        if(is.na(window[1])) {
          cat("\n\nTHIS IS NA:", paste(df[1,]))
          tel.length <- NULL
        }
      }
    }
    
    # Return telomere value
    # cat("\n\tThis is telomere and read length:", round(tel.length), "-", read.length)
    trunc.s <- 1
    trunc.e <- read.length - tel.length
    # cat("\n\tThis is the range of the read which does not contain telomeric repteats:", trunc.s, "-", trunc.e, "\n\n")
    if (length(tel.length) == 0 || tel.length == 0) {tel.length = NULL}
    if (!is.null(tel.length)) {
      return(c(tel.length=round(tel.length), read.length=read.length, 
               trunc.start=trunc.s, trunc.end=trunc.e, 
               telomere.end=tel.end))
    }
  }
  
  # Check if no condition was met
  if(is.null(tel.length)) {
    # cat("\n\n\nNO CONDITION WAS MET!! FOR READ:", paste(df[1,]))
  }
}