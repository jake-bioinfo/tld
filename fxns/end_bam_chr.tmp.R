end_bam_chr <- function(df) {
  
  end <- df$telomere.end
  
  # Determine is the end falls on the left or right side of the chromosome
  left <- (end == "L")
  right <- (end == "R")
  
  
  # Assign endedness and return value
  if (left) {
    endedness <- "5'"
  }
  
  if (right) {
    endedness <- "3'"
  }
  
  if (!left & !right) {
    endedness <- NA
  }
  
  return(c(chr.end = endedness))
}