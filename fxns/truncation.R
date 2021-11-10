truncation <- function(read.df, vl) {
  require(Biostrings)
  
  # Get initial values
  read.name <- read.df$r.name
  st <- read.df$trunc.start
  en <- read.df$trunc.end
  name.index <- which(names(vl) == read.name)
  
  read.seq <- DNAStringSet(vl[name.index], start = as.numeric(st), 
                           end = as.numeric(en), use.names = TRUE)
  read.seq <- data.frame(read.seq)
  
  colnames(read.seq) <- "trunc.seq"
  read.seq$trunc.id <- names(vl[name.index])
 
  return(read.seq)
}