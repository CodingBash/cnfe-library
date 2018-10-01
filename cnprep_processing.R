#
# This is for generating the chromosomeSizes - somehow the genome needs to be passed in by the user.
#
library(CNprep)

cnprep_process <- function(standardCopyNumberMap, parallel = TRUE, mclust_model = "E", minjoin = 0, ntrial = 10){
  normal_samples <- retrieveTargetSampleList(standardCopyNumberMap, class="N") # May not even need this
  reference <- standardCopyNumberMap@reference
  normal_segments <- retrieveTargetSampleListSegments(standardCopyNumberMap, class="N", getAbsolute = TRUE)
  norminput <- retrieveNormInput(normal_segments)
  norminput <- filterNormInput(norminput=norminput, length_threshold = 10000000)
  target_samples <- retrieveTargetSampleList(standardCopyNumberMap, class=c("T", "F", "M"))
  
  cnprep_inputs <- lapply(target_samples, function(target_sample){
    standardCopyNumberProfile <- standardCopyNumberMap@map[[target_sample]]
    seginput <- transformSegmentsForCNprep(standardCopyNumberProfile, target_sample)
    ratinput <- transformRatioForCNprep(standardCopyNumberProfile, target_sample)
    return(list("seginput" = seginput, "ratinput" = ratinput))
  })
  names(cnprep_inputs) <- target_samples
  
  segtable_binded <- data.frame(stringsAsFactors = FALSE)
  if(parallel == TRUE){
    seginput_binded <- do.call(rbind, lapply(cnprep_inputs, function(cnprep_input){cnprep_input[["seginput"]]}))
    ratinput_binded <- do.call(cbind, lapply(cnprep_inputs, function(cnprep_input){cnprep_input[["ratinput"]]}))
    try({
      print("CNprep: Starting processing for all samples in parallel")
      segtable_binded <- runCNpreprocessing(seginput = seginput_binded, ratinput = ratinput_binded, norminput = norminput, modelNames = mclust_model, minjoin = minjoin, ntrial = ntrial)
      print("CNprep: Completed processing for all samples")
      }) # TODO: Verify exception handling is performed correctly
  } else {
    segtable_binded <- do.call(rbind, lapply(cnprep_inputs, function(cnprep_input){
      seginput <- cnprep_input[["seginput"]]
      ratinput <- cnprep_input[["ratinput"]]
      
      print(paste0("CNprep: Starting processing for sample: ",  seginput[1, 1])) # TODO: Get sample from list name, not table
      segtable <- data.frame(stringsAsFactors = FALSE)
      try({
        segtable <- runCNpreprocessing(seginput = seginput, ratinput = ratinput, norminput = norminput, modelNames = mclust_model, minjoin = minjoin, ntrial = ntrial)
        print(paste0("CNprep: Processed sample: ",  seginput[1, 1])) # TODO: Get sample from list name, not table
        return(segtable)
      }) # TODO: Verify exception handling is performed correctly
      print(paste0("CNprep: Sample: ",  seginput[1, 1], " failed!")) # TODO: Get sample from list name, not table
      return(segtable)
    }))
  }
  return(segtable_binded)
}

transformSegmentsForCNprep <- function(standardCopyNumberProfile, sample){
  chromosomalSegments <- standardCopyNumberProfile@chromosomalSegments
  absoluteSegments <- standardCopyNumberProfile@absoluteSegments 
  probes <- standardCopyNumberProfile@probes
  
  seginput <- data.frame(ID = sample, start = probes[[2]], end = probes[[3]], seg.median = absoluteSegments[[4]], 
                         chrom = absoluteSegments[[1]], chrom.pos.start = chromosomalSegments[[2]], 
                         chrom.pos.end = chromosomalSegments[[3]], abs.pos.start = absoluteSegments[[2]],
                         abs.pos.end = absoluteSegments[[3]], stringsAsFactors = FALSE)
  
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
# Retrieve the norminput argument for CNprep::CNpreprocessing()
#
retrieveNormInput <- function(normal_segments){
  norminput <- do.call(rbind, lapply(seq(1, nrow(normal_segments)), function(normal_segments.index){
    norminput.entry <- data.frame(length = normal_segments[normal_segments.index, 3] - normal_segments[normal_segments.index, 2], segmedian = normal_segments[normal_segments.index, 4], stringsAsFactors = FALSE)
    return(norminput.entry)
  }))
  return(norminput)
}

#
# Filter norminput from artifacts
#
filterNormInput <- function(norminput, length_threshold=10000000){
  # TODO: Algorithm to dynamically determine length_threshold
  filtered_norminput <- norminput[norminput$length > length_threshold & abs(norminput$segmedian) < 0.5,]
  return(filtered_norminput)
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