cnprep_process <- function(){
  normal_samples <- NA # List of strings containing IDs of normal diploid sequences
  reference <- NA # ID of normal reference
  cytobands <- NA # Retrieve cytobands
  chromosome_sizes <- NA # Retrieve chromosome sizes (preferably dynamically)
  
  normal_segments <- NA # Segment data of all IDs in normal_samples
  
  target_samples <- NA # List of string containing IDs of target diploid sequences
  
  for(target_samples.i in seq_along(target_samples)){
    segment_data <- NA # segment data for sample.i
    ratio_data <- NA # ratio data for sample.i
    
    # . . .
    
  }
  
  # . . .
  
  segtable <- NA # segtable of all results
  return(segtable)
}