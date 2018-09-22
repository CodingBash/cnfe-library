#
# This is for generating the chromosomeSizes - somehow the genome needs to be passed in by the user.
#
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)

cnprep_process <- function(standardCopyNumberMap){
  normal_samples <- retrieveTargetSampleList(standardCopyNumberMap, class="N") # May not even need this
  reference <- standardCopyNumberMap@reference
  cytobands <- NA # TODO: See if I can do CNprep without cytobands
  chromosome_sizes <- generateChromosomeSizes(genome = BSgenome.Hsapiens.UCSC.hg19) # TODO Temporary placement
  
  normal_segments <- retrieveTargetSampleListSegments(standardCopyNumberMap, class="N")
  
  target_samples <- retrieveTargetSampleList(standardCopyNumberMap, class=c("T", "F", "M"))
  
  for(target_samples.i in seq_along(target_samples)){
    sample <- target_samples[target_samples.i]
    segment_data <- standardCopyNumberMap@map[[sample]]@segments
    ratio_data <- standardCopyNumberMap@map[[sample]]@ratio
    
    seginput <- transformSegmentsForCNprep(segment_data, sample, chromosome_sizes)
    ratinput <- transformRatioForCNprep(ratio_data, sample)
    # TODO: Get seginput, ratinput
    
    # . . a
    
  }
  
  # . . .
  
  segtable <- NA # segtable of all results
  return(segtable)
}

#
# Retrieve information of cytobands (downloaded for genome hg19)
#
retrieveCytobands <- function(dir = "cytoBand.txt"){
  cytobands <- read.table(dir, header=F, sep = "\t", stringsAsFactors = F)
  names(cytobands) <- c("chrom", "start", "end", "cytoloc", "stain")
  return(cytobands)
}

#
# With the chromosome and chromosome location, retrieve the cytoband location
#
findCytolocation <- function(cytobands, chrom, chrom.position){
  chrom = if(chrom == 23) "X" else if (chrom == 24 ) "Y" else chrom
  row <- cytobands[cytobands$chrom == paste("chr", chrom, sep = "") & cytobands$start <= chrom.position & cytobands$end >= chrom.position, ]
  returnme <- data.frame(row)$cytoloc
  return(returnme)
}


#
# Retrieve the seginput argument for CNprep::CNpreprocessing()
#
# TODO: Put in its own library?
transformSegmentsForCNprep <- function(segment_data, sample, chromosomeSizes){
  seginput <- data.frame(stringsAsFactors = FALSE)
  # Iterate through each segment
  for(segment_data.index in seq(1, nrow(segment_data))){
    # Get absolute position of segment
    abs_position <- chromsomeToAbsoluteBPConversionForSingleEntry(segment_data[segment_data.index,]$X.chrom., segment_data[segment_data.index,]$X.start., segment_data[segment_data.index,]$X.end., chromosomeSizes)
    
    probes.start = 0
    probes.end = 0
    for(i in seq(1, segment_data.index)){
      if(i != segment_data.index){
        probes.start <- probes.start + segment_data[i,]$X.num.mark.
      }
      probes.end <- probes.end + segment_data[i,]$X.num.mark.
    }
    probes.start <- probes.start + 1
    
    
    seginput.entry <- data.frame(ID = sample, start = probes.start, end = probes.end, 
                                 num.probes = segment_data[segment_data.index,]$X.num.mark., seg.median = segment_data[segment_data.index,]$X.cnlr.median., 
                                 chrom = segment_data[segment_data.index,]$X.chrom., chrom.pos.start = segment_data[segment_data.index,]$X.start., 
                                 chrom.pos.end = segment_data[segment_data.index,]$X.end., cytoband.start = cytoband.my.start, 
                                 cytoband.end = cytoband.my.end, abs.pos.start = abs_position$start,
                                 abs.pos.end = abs_position$end)
    
    
    seginput <- rbind(seginput, seginput.entry)
  }
  return(seginput)
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
