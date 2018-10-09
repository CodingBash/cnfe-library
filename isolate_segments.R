isolate <- function(segtable, criteria_function, probes){
  segtable_as_list <- split(segtable, seq(nrow(segtable)))
  selected_segtable <- segtable[sapply(segtable_as_list, criteria_function),]
  
  isolationResults <- NA
  if (probes == TRUE){
    # TODO: Decide between abs or chrom
    # TODO: Allow user to choose segmedian or maxzmean
    isolationResults <- selected_segtable[,c("chrom", "abs.pos.start", "abs.pos.end", "segmedian")]
  } else {
    isolationResults <- selected_segtable[,c("chrom", "start", "end", "segmedian")]
  }
  return(isolationResults)
}