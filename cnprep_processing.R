cnprep_process <- function(standardCopyNumberMap){
  normal_samples <- retrieveTargetSampleList(standardCopyNumberMap, class="N")
  reference <- standardCopyNumberMap@reference
  cytobands <- NA # Retrieve cytobands
  chromosome_sizes <- NA # Retrieve chromosome sizes (preferably dynamically) - user should send in genome (i.e. hg19)
  
  normal_segments <- retrieveTargetSampleListSegments(standardCopyNumberMap, class="N")
  
  target_samples <- retrieveTargetSampleList(standardCopyNumberMap, class=c("T", "F", "M"))
  
  for(target_samples.i in seq_along(target_samples)){
    sample <- target_samples[target_samples.i]
    segment_data <- standardCopyNumberMap@map[[sample]]@segments
    ratio_data <- standardCopyNumberMap@map[[sample]]@ratio
    
    # TODO: Get seginput, ratinput
    
    # . . .
    
  }
  
  # . . .
  
  segtable <- NA # segtable of all results
  return(segtable)
}