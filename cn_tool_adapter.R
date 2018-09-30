setwd("C:/Users/bbece/Documents/Git-Projects/Git-Research-Projects/cnfe-library")
source("./class_definitions.R")
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)

# Input: Object containing a list of FACETS objects (facetsG5XX, facetsG5OO, facetsG5Fit) mapped to sample ID
# Ouptput: Datastructure that CNprep can resolve
facets_adapter <- function(facetsCopyNumberResults, chromosomeSizes){
  generateProbes <- function(fit){
    probes.df <- data.frame(stringsAsFactors = FALSE)
    for(fit.index in seq(1, nrow(fit))){
      probes.chrom <- fit[fit.index, ]$chrom
      probes.start <- 0
      probes.end <- 0
      for(i in seq(1, fit.index)){
        if(i != fit.index){
          probes.start <- probes.start + fit[i, ]$num.mark
        }
        probes.end <- probes.end + fit[i, ]$num.mark
      }
      probes.start <- probes.start + 1
      probes.df.row <- data.frame(chrom = probes.chrom, start = probes.start, end = probes.end)
      probes.df <- rbind(probes.df, probes.df.row)
    }
    return(probes.df)
  }
  standardCopyNumberMapList <- lapply(names(facetsCopyNumberResults), function(reference){
    standardCopyNumberMap <- new("ReferencedCopyNumberMap", reference=reference, chromosomeSizes=chromosomeSizes)  
    sapply(ls(facetsCopyNumberResults[[reference]]@map), function(target){
      print(paste0("Converting target = ", target))
      facetsProfile <- facetsCopyNumberResults[[reference]]@map[[target]]
      
      ratio <- select(facetsProfile@xx, chrom, maploc, cnlr)
      colnames(ratio) <- c("chr", "start", "cnlr")
      ratio$end <- ratio$start
      ratio <- select(ratio, chr, start, end, cnlr)
      
      segments <- select(facetsProfile@fit, chrom, start, end, cnlr.median)
      colnames(segments) <- c("chr", "start", "end", "cnlr")
      
      addSegments(map=standardCopyNumberMap, target=target, segments=segments, isAbsolute=FALSE)
      standardCopyNumberMap@map[[target]]@chromosomalRatio=ratio # TODO: Figure out solution for genomic conversion
      
      
      # solveProbes(map=standardCopyNumberMap, target=target) # TODO: This solves probes based on ratio/segments matching
      standardCopyNumberMap@map[[target]]@probes <- generateProbes(facetsProfile@fit) # TODO: This solves probes based on facets num.mark (shown to be inaccurate)
      
    })
    return(standardCopyNumberMap)
  })
  names(standardCopyNumberMapList) <- names(facetsCopyNumberResults)
  return(standardCopyNumberMapList)
}

test_facets_adapter <- function(){
  facetsCopyNumberResults <- retrieve_legacy_facets_objects()
  chromosomeSizes <- generateChromosomeSizes(genome = BSgenome.Hsapiens.UCSC.hg19)
  standardCopyNumberMapList <- facets_adapter(facetsCopyNumberResults, chromosomeSizes)
}