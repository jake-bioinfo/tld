end_bam_chr <- function(df, mp) {
  
  # Set initial values
  chr <- as.numeric(gsub(pattern = "chr_", replacement = "", df$seqnames))
  st <- as.numeric(df$start)
  mp <- mp
  
  # Check if irradiation status exists, if so then set mp by that calculation
  if ("i.status" %in% names(df)) {
    irr <- df$s.name=="irr"
    wt <- df$s.name=="3d7"
    
    # Case for irradiated samples blast hits
    if (irr) {
      chr.3d7.l <- list(504382, 746179, 986745, 1214710, 1348989,
                        1422185, 1448692 , 1481828 , 1497267, 1702824,
                        2040918, 2287845, 2943957, 3295951, 1447516, 
                        16897, 67395, 11705, 61840, 19426, 63493)
      c <- df$i.d
      
      if (length(c)!=0&(c=="16"|c=="17"|c=="18"|c=="19"|c=="20"|c=="21")) {
        mp <- 0.5*chr.3d7.l[[chr]]
      } 
    }
    
    # Case for non irradiated assembly
    if (wt) {
      chr.3d7.l <- list(574122, 774001, 930220, 1235345, 
                        1331859, 1268092, 1444814, 1479982, 
                        1525650, 1268186, 1894971, 2300568, 
                        2900879, 3292472)
      
    }
  
    } else {
    
    chr.3d7.l <- list(640851, 947102, 1067971, 1200490, 1343557, 
                      1418242, 1445207, 1472805, 1541735, 1687656, 
                      2038340, 2271494, 2925236, 3291936)
    
  }
  
  chr.l <- chr.3d7.l[[chr]]
  #mp <- 0.5*as.numeric(chr.l)
  left <- (st < mp)
  right <- (st > (chr.l - mp))
  # cat("\n\nchr.l:", chr.l,
  #     "\n\tleft:", left,
  #     "\n\tst:", st,
  #     "\n\tmp:", mp, "\n\n")
  
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