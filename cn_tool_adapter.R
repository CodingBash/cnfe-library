setwd("C:/Users/bbece/Documents/Git-Projects/Git-Research-Projects/cnfe-library")
source("./class_definitions.R")
library(dplyr)

retrieve_legacy_facets_objects <- function(){
  #
  # Take directory of FACETS objects and write to tsv
  #
  
  setwd(paste0("~/Git-Projects/Git-Research-Projects/FACETS_write_files/")) 
  source("helperFunctions.R")
  source("facetsAnalysisLibrary.R")
  library(stringr)
  
  fileDf<- read.table("./resources/listR1File.txt", sep="\t")
  fileList <- fileDf$V1
  names(fileList) <- rownames(fileList)
  
  #
  # Determine reference sample
  #
  reg <- "_([^/]+)/"
  normal_ref <- 40 # This is the reference sample chosen
  ref_match <- str_extract(fileList[[normal_ref]], reg)[1]
  ref_sample <- substr(ref_match, 2, nchar(ref_match) - 1)
  xxNew <- readRDS(paste0("./resources/xxJoinsegNew_r", normal_ref, ".rds"))
  facets_dir <- paste0("output/", "FACETS_Reference_", ref_sample, "_8_2_18_1/")
  
  facetsCopyNumberResults <- list()
  facetsCopyNumberResults[[ref_sample]] <- new("ReferencedFacetsOutputMap", reference=ref_sample)
  
  #dir.create(file.path(facets_dir), showWarnings = FALSE)
  for(fileList.index in seq(1, length(fileList))){
    #
    # Pre-work
    #
    if(fileList.index == normal_ref){
      next
    }
    target_match <- str_extract(fileList[[fileList.index]], reg)[1]
    target_sample <- substr(target_match, 2, nchar(target_match) - 1)
    target_class <- substr(target_sample, 2,2)
    
    res_dir <- if(target_class == "N") "facetsRNormal" else "facetsRCancer"
    res_dir <- paste0("resources/", res_dir)
    
    
    
    #
    # Save to directory
    #
    out_dir <- paste0(facets_dir, "Sample_", target_sample)
    #dir.create(file.path(out_dir), showWarnings = FALSE)
    xxOutputFilename <- paste0(out_dir, "/", target_sample, "--", ref_sample, ".procSample-jseg.cnv.facets.v0.5.2.txt")
    fitOutputFilename <- paste0(out_dir, "/", target_sample, "--", ref_sample, ".cnv.facets.v0.5.2.txt")
    try({
      #
      # Load SNPs
      #
      saveSnps <- NA
      if(target_class == "N"){
        xxFilename <- paste0(res_dir, "/facetsG5XX_", normal_ref, "_", fileList.index, ".rds")
        xx <- readRDS(xxFilename)  
        saveSnps <- xx$jointseg
      } else if (target_class %in% c("T", "F", "M")){
        saveSnps <- xxNew[[paste0("n", normal_ref, "_", fileList.index)]]    
      } else {
        print(paste0("WARNING: Target class not recognized: ", target_class))
        return()
      }
      
      addRatio(map=facetsCopyNumberResults[[ref_sample]], target=target_sample, ratio=saveSnps)
      #write.table(saveSnps, file = xxOutputFilename, quote=TRUE, sep = "\t", row.names = FALSE  )
      print(paste0("Wrote SNPS to: ", xxOutputFilename))
    }, silent = TRUE)
    
    try({
      #
      # Load segment fit
      #
      fitFilename <- paste0(res_dir, "/facetsG5Fit_", normal_ref, "_", fileList.index, ".rds")
      fit <- readRDS(fitFilename)
      saveFit <- fit$cncf
      
      addSegments(map=facetsCopyNumberResults[[ref_sample]], target=target_sample, segments=saveFit)
      #write.table(saveFit, file = fitOutputFilename, quote=TRUE, sep = "\t", row.names = FALSE  )
      print(paste0("Wrote FIT to: ", fitOutputFilename))
    }, silent = TRUE)
  }
  
  return(facetsCopyNumberResults)
  
}

# Input: Object containing a list of FACETS objects (facetsG5XX, facetsG5OO, facetsG5Fit) mapped to sample ID
# Ouptput: Datastructure that CNprep can resolve
facets_adapter <- function(facetsCopyNumberResults){
  standardCopyNumberMapList <- lapply(names(facetsCopyNumberResults), function(reference){
    standardCopyNumberMap <- new("ReferencedCopyNumberMap", reference=reference)  
    invisible(sapply(ls(facetsCopyNumberResults[[reference]]@map), function(target){
      facetsProfile <- facetsCopyNumberResults[[reference]]@map[[target]]
      
      ratio <- select(facetsProfile@xx, chrom, maploc, cnlr)
      colnames(ratio) <- c("chr", "start", "cnlr")
      ratio$end <- ratio$start
      ratio <- select(ratio, chr, start, end, cnlr)
      
      segments <- select(facetsProfile@fit, chrom, start, end, cnlr.median)
      colnames(segments) <- c("chr", "start", "end", "cnlr")
      
      profile <- new("CopyNumberProfile", ratio=ratio, segments=segments)
      addEntry(map=standardCopyNumberMap, target=target, profile=profile)  
    }))
    return(standardCopyNumberMap)
  })
  names(standardCopyNumberMapList) <- names(facetsCopyNumberResults)
  return(standardCopyNumberMapList)
}

test_facets_adapter <- function(){
  facetsCopyNumberResults <- retrieve_legacy_facets_objects()
  standardCopyNumberMapList <- facets_adapter(facetsCopyNumberResults)
}
test_facets_adapter()