isolate <- function(segtable, criteria_function){
  segtable_as_list <- split(segtable, seq(nrow(segtable)))
  selected_segtable <- segtable[sapply(segtable_as_list, criteria_function),]
  
  isolationResults <- NA
  return(isolationResults)
}