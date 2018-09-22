#
# This is for generating the chromosomeSizes - somehow the genome needs to be passed in by the user.
#
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
library(CNprep)

cnprep_process <- function(standardCopyNumberMap){
  normal_samples <- retrieveTargetSampleList(standardCopyNumberMap, class="N") # May not even need this
  reference <- standardCopyNumberMap@reference
  normal_segments <- retrieveTargetSampleListSegments(standardCopyNumberMap, class="N", getAbsolute = TRUE)
  
  target_samples <- retrieveTargetSampleList(standardCopyNumberMap, class=c("T", "F", "M"))
  
  for(target_samples.i in seq_along(target_samples)){
    sample <- target_samples[target_samples.i]
    standardCopyNumberProfile <- standardCopyNumberMap@map[[sample]]
    
    # TODO: Create one huge seginput (I have code somewhere on drug-response-prediction)
    # TODO: Allow users to choose between parallel and syncrynous CNprep run -> for fault continuation if one of samples are bad.
    # TODO: Among whole package, need to focus on exception handling.
    
    seginput <- transformSegmentsForCNprep(standardCopyNumberProfile, sample)
    ratinput <- transformRatioForCNprep(standardCopyNumberProfile, sample)
    # TODO: Get seginput, ratinput
    
    try({
      # TODO: Add parameters as function inputs - allow all parameters
      segtable <- runCNpreprocessing(seginput = seginput, ratinput = ratinput, norminput = norminput, modelNames = mclust_model, minjoin = minjoin, ntrial = ntrial) #TODO: Is there a distrib="Grid"?
      print(paste("Produced segtable for sample", sample))
      
      write.table(segtable, paste("./output/", output_dir,"/", sample, "_segtable.tsv", sep = ""), row.names = F, sep = "\t", quote = FALSE)
      print(paste("Wrote output for sample", sample))
    }, silent = TRUE)
    # . . a
    
  }
  
  # . . .
  
  segtable <- NA # segtable of all results
  return(segtable)
}

transformSegmentsForCNprep <- function(standardCopyNumberProfile, sample){
  seginput <- data.frame(stringsAsFactors = FALSE)
  chromosomalSegments <- standardCopyNumberProfile@chromosomalSegments
  absoluteSegments <- standardCopyNumberProfile@absoluteSegments 
    
  # Iterate through each segment
  for(segment_data.index in seq(1, nrow(standardCopyNumberProfile@chromosomalSegments))){
     
    # TODO: Support probes by including in adapter and object design (need to really consider object design (need to improve S4 class skills))   
    
    # TODO: Support for combinatorial missings of standardCopyNumberProfile object - some users might not have all data available - only require what is necessary (through S4 object validation)
    seginput.entry <- data.frame(ID = sample, seg.median = absoluteSegments[[4]], 
                                 chrom = absoluteSegments[[1]], chrom.pos.start = chromosomalSegments[[2]], 
                                 chrom.pos.end = chromosomalSegments[[3]], abs.pos.start = absoluteSegments[[2]],
                                 abs.pos.end = absoluteSegments[[3]])
    
    seginput <- rbind(seginput, seginput.entry)
  }
  return(seginput)
}

transformRatioForCNprep <- function(standardCopyNumberProfile, sample){
  ratio <- NA
  if(nrow(standardCopyNumberProfile@absoluteRatio) == 0){
    ratio <- standardCopyNumberProfile@chromosomalRatio    
  } else {
    ratio <- standardCopyNumberProfile@absoluteRatio
  }
  ratinput <- data.frame(ratio[[4]])
  names(ratinput) <- c(sample)
  return(ratinput)
}

#
# Wrapper function for CNprep::CNpreprocessing()
#
runCNpreprocessing <- function(seginput, ratinput, norminput,
                               blsize=50, minjoin=0.25, cweight=0.4, bstimes=50,
                               chromrange=1:22, distrib="Rparallel", njobs=4, modelNames="E", ntrial= 10){
  segtable<-CNpreprocessing(segall=seginput,ratall=ratinput,"ID","start","end",
                            chromcol="chrom",bpstartcol="chrom.pos.start",bpendcol="chrom.pos.end",blsize=blsize,
                            minjoin=minjoin,cweight=cweight,bstimes=bstimes,chromrange=chromrange,distrib=distrib,njobs=njobs,
                            modelNames=modelNames,normalength=norminput[,1],normalmedian=norminput[,2], ntrial=ntrial)
  return(segtable)
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
