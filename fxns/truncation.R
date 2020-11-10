truncation <- function(read.df, vl) {
  require(Biostrings)
  
  # Get initial values
  read.name <- read.df$r.name
  st <- read.df$trunc.start
  en <- read.df$trunc.end
  name.index <- which(names(vl) == read.name)
  #cat("\nThis is read.name:", read.name,
  #   "\n\tThis is trunc start:", st,
  #   "\n\tThis is trunc end:", en,
  #   "\n\tThis is read.name index:", name.index, "\n\n")

  #print(head(vl))
  #print(vl[name.index])
  #cat("\n\nThis is the read line:\n")
  #print(read.df)
  
  read.seq <- DNAStringSet(vl[name.index], start = as.numeric(st), 
                           end = as.numeric(en), use.names = TRUE)
  read.seq <- data.frame(read.seq)
  
  #cat("\n\nThis is the subsetted read sequence:\n")
  colnames(read.seq) <- "trunc.seq"
  read.seq$trunc.id <- names(vl[name.index])
  #print(read.seq$trunc.seq)
  #cat("\n\nThis is dim:\n")
  #print(dim(read.seq))
  
  return(read.seq)
}