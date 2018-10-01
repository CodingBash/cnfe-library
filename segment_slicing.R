slice_segments <- function(){
  sample_list <- NA # IDs of all samples included in the segtable
  segtable <- NA # segtable concatenation of all CNprep results for all samples in sample_list
  
  # . . . slicing
  
  organoidSlices <- NA # slicing results
  
  # . . . convert slicing output to bed
  
  # Slicing BED files 
  probeBed <- NA
  genomeBed <- NA
}

source("./slicing_resources/rfunctions/peaks.valleys.AK1.ssc")
source("./slicing_resources/rfunctions/set.tiers.new.R")
source("./slicing_resources/TCGA/newMatchedMasking")
source("./slicing_resources/slicing/sliceHT.R")
source("./slicing_resources/slicing/convertSlicingOutputToBed.R")
slice <- function(segtable){
  jointslices <- sliceHT(segtable)
}