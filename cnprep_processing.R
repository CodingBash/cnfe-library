#
# This is for generating the chromosomeSizes - somehow the genome needs to be passed in by the user.
#
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)

cnprep_process <- function(standardCopyNumberMap){
  normal_samples <- retrieveTargetSampleList(standardCopyNumberMap, class="N") # May not even need this
  reference <- standardCopyNumberMap@reference
  normal_segments <- retrieveTargetSampleListSegments(standardCopyNumberMap, class="N")
  
  target_samples <- retrieveTargetSampleList(standardCopyNumberMap, class=c("T", "F", "M"))
  
  for(target_samples.i in seq_along(target_samples)){
    sample <- target_samples[target_samples.i]
    standardCopyNumberProfile <- standardCopyNumberMap@map[[sample]]
    
    seginput <- transformSegmentsForCNprep(standardCopyNumberProfile, sample)
    ratinput <- transformRatioForCNprep(standardCopyNumberProfile, sample)
    # TODO: Get seginput, ratinput
    
    # . . a
    
  }
  
  # . . .
  
  segtable <- NA # segtable of all results
  return(segtable)
}

transformSegmentsForCNprep <- function(standardCopyNumberProfile, sample){
  seginput <- data.frame(stringsAsFactors = FALSE)
  chromosomalSegments <- standardCopyNumberProfile@chromosomalSegments
  absoluteSegments <- standardCopyNumberProfile@chromosomalSegments 
    
  # Iterate through each segment
  for(segment_data.index in seq(1, nrow(standardCopyNumberProfile@chromosomalSegments))){
     
    # TODO: Support probes by including in adapter and object design (need to really consider object design (need to improve S4 class skills))   
    
    seginput.entry <- data.frame(ID = sample, seg.median = absoluteSegments[[4]], 
                                 chrom = absoluteSegments[[1]], chrom.pos.start = chromosomalSegments[[2]], 
                                 chrom.pos.end = chromosomalSegments[[3]], abs.pos.start = absoluteSegments[[2]],
                                 abs.pos.end = absoluteSegments[[3]])
    
    seginput <- rbind(seginput, seginput.entry)
  }
  return(seginput)
}

transformRatioForCNprep <- function(standardCopyNumberProfile, sample){
  ratinput <- data.frame(standardCopyNumberProfile@absoluteRatio[[4]])
  names(ratinput) <- c(sample)
  return(ratinput)
}

test_cnprep_process <- function(){
  setwd("C:/Users/bbece/Documents/Git-Projects/Git-Research-Projects/cnfe-library")
  source("./class_definitions.R")
  source("./cn_tool_adapter.R")
  facetsCopyNumberResults <- retrieve_legacy_facets_objects()
  standardCopyNumberMapList <- facets_adapter(facetsCopyNumberResults)
  standardCopyNumberMap <- standardCopyNumberMapList[['hN30']]
  cnprep_process(standardCopyNumberMap)
}
test_cnprep_process()
