source("./slicing_resources/rfunctions/peaks.valleys.AK1.ssc")
source("./slicing_resources/rfunctions/set.tiers.new.R")
source("./slicing_resources/TCGA/newMatchedMasking")
source("./slicing_resources/slicing/sliceHT.R")
source("./slicing_resources/slicing/convertSlicingOutputToBed.R")
slice <- function(segtable, probes, amplification){
  # Compute slices
  jointslices <- sliceHT(segtable)
  
  # Format slices
  slicesFormatted <- formatSlicingOutput(jointslices)
  
  # Concatenate slices from all samples with user specifications (unit and slice type)
  slicesResults <- do.call(rbind, lapply(slicesFormatted, function(sampleSlices){
    targetSlices <- if (probes == TRUE) sampleSlices$probeSlices else sampleSlices$genomeSlices
    target.lesion.type <- if (amplification == TRUE) 1 else 0
    return(targetSlices[targetSlices$lesion.type == target.lesion.type, ])
  }))
  
  return(slicesResults)
}