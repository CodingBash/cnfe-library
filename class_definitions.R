setwd("C:/Users/bbece/Documents/Git-Projects/Git-Research-Projects/cnfe-library")
source("./utility_functions.R")

defineGenerics <- function(){
  setGeneric("addEntry", function(map, target, profile) {
    standardGeneric("addEntry")
  })
  setGeneric("addRatio", function(map, target, ratio) {
    standardGeneric("addRatio")
  })
  setGeneric("addRatio", function(map, target, ratio, isAbsolute) {
    standardGeneric("addRatio")
  })
  setGeneric("addRatio", function(map, target, ratio, chromosomeSizes, isAbsolute) {
    standardGeneric("addRatio")
  })
  setGeneric("addSegments", function(map, target, segments) {
    standardGeneric("addSegments")
  })
  setGeneric("addSegments", function(map, target, segments, isAbsolute) {
    standardGeneric("addSegments")
  })
  setGeneric("addSegments", function(map, target, segments, chromosomeSizes, isAbsolute) {
    standardGeneric("addSegments")
  })
  setGeneric("removeEntry", function(map, target) {
    standardGeneric("removeEntry")
  })
  setGeneric("solveProbes", function(map, target) {
    standardGeneric("solveProbes")
  })
  setGeneric("retrieveTargetSampleList", function(map) {
    standardGeneric("retrieveTargetSampleList")
  })
  setGeneric("retrieveTargetSampleList", function(map, class) {
    standardGeneric("retrieveTargetSampleList")
  })
  setGeneric("retrieveTargetSampleListSegments", function(map, class, getAbsolute) {
    standardGeneric("retrieveTargetSampleListSegments")
  })
}

defineFacetsClass <- function(){
  setClass("FacetsOutputContainer", representation(xx="data.frame", fit="data.frame"))
  setClass("ReferencedFacetsOutputMap", representation(map="environment", reference="character"),
           prototype(map=new.env(), reference=NA_character_))
  setMethod("addEntry", signature(map = "ReferencedFacetsOutputMap", target="character", profile="FacetsOutputContainer"), function(map, target, profile){ 
    map@map[[target]] <- profile
  })
  setMethod("addRatio", signature(map = "ReferencedFacetsOutputMap", target="character", ratio="data.frame"), function(map, target, ratio){ 
    if(is.null(map@map[[target]])) {
      outputContainer <- new("FacetsOutputContainer", xx=ratio)
      map@map[[target]] <- outputContainer
    } else {
      outputContainer <- map@map[[target]]
      outputContainer@xx <- ratio
      map@map[[target]] <- outputContainer
    }
  })
  setMethod("addSegments", signature(map = "ReferencedFacetsOutputMap", target="character", segments="data.frame"), function(map, target, segments){ 
    if(is.null(map@map[[target]])) {
      outputContainer <- new("FacetsOutputContainer", fit=segments)
      map@map[[target]] <- outputContainer
    } else {
      outputContainer <- map@map[[target]]
      outputContainer@fit <- segments
      map@map[[target]] <- outputContainer
    }
  })
  setMethod("removeEntry", signature(map = "ReferencedFacetsOutputMap", target="character"), function(map, target){ 
    remove(list=c(target), envir=map@map)
  })
} 
defineStandardClass <- function(){
  # TODO: Add validation on BED data.frame
  setClass("CopyNumberProfile", representation(probes="data.frame", chromosomalSegments="data.frame", chromosomalRatio="data.frame", absoluteSegments="data.frame", absoluteRatio="data.frame"))
  setClass("ReferencedCopyNumberMap", representation(map="environment", reference="character", chromosomeSizes="data.frame"),
           prototype(map=new.env(), reference=NA_character_))
  setMethod("addEntry", signature(map = "ReferencedCopyNumberMap", target="character", profile="CopyNumberProfile"), function(map, target, profile){ 
    map@map[[target]] <- profile
  })
  
  setMethod("removeEntry", signature(map = "ReferencedCopyNumberMap", target="character"), function(map, target){ 
    remove(list=c(target), envir=map@map)
  })
  
  setMethod("addRatio", signature(map = "ReferencedCopyNumberMap", target="character", ratio="data.frame", chromosomeSizes="data.frame", isAbsolute="logical"), function(map, target, ratio, chromosomeSizes, isAbsolute){ 
    chromosomalRatio <- NA
    absoluteRatio <- NA
    
    
    if(isAbsolute == TRUE){
      absoluteRatio <- ratio
      chromosomalRatio <- genomicConversion(absoluteRatio, chromosomeSizes, chromosomeToAbsolute = FALSE)
    } else {
      chromosomalRatio <- ratio
      absoluteRatio <- genomicConversion(chromosomalRatio, chromosomeSizes, chromosomeToAbsolute = TRUE)
    }
    
    if(is.null(map@map[[target]])) {
      profile <- new("CopyNumberProfile", chromosomalRatio=chromosomalRatio, absoluteRatio=absoluteRatio)
      map@map[[target]] <- profile
    } else {
      profile <- map@map[[target]]
      profile@chromosomalRatio <- chromosomalRatio
      profile@absoluteRatio <- absoluteRatio
      map@map[[target]] <- profile
    }
  })
  
  setMethod("addRatio", signature(map = "ReferencedCopyNumberMap", target="character", ratio="data.frame", isAbsolute="logical"), function(map, target, ratio, isAbsolute){ 
    if(!is.null(map@chromosomeSizes)){
      chromosomeSizes = map@chromosomeSizes
      addRatio(map = map, target = target, ratio = ratio, chromosomeSizes = chromosomeSizes, isAbsolute = isAbsolute)
    } else {
      stop("Add chromosomeSizes before adding ratios") # TODO: Bad design: should not require chromosomeSizes to be added beforehand
    }
  })
  
  
  setMethod("addSegments", signature(map = "ReferencedCopyNumberMap", target="character", segments="data.frame", chromosomeSizes="data.frame", isAbsolute="logical"), function(map, target, segments, chromosomeSizes, isAbsolute){ 
    chromosomalSegments <- NA
    absoluteSegments <- NA
    
    if(isAbsolute == TRUE){
      absoluteSegments <- segments
      chromosomalSegments <- genomicConversion(absoluteSegments, chromosomeSizes, chromosomeToAbsolute = FALSE)
    } else {
      chromosomalSegments <- segments
      absoluteSegments <- genomicConversion(chromosomalSegments, chromosomeSizes, chromosomeToAbsolute = TRUE)
    }
    
    if(is.null(map@map[[target]])) {
      profile <- new("CopyNumberProfile", chromosomalSegments=chromosomalSegments, absoluteSegments=absoluteSegments)
      map@map[[target]] <- profile
    } else {
      profile <- map@map[[target]]
      profile@chromosomalSegments <- chromosomalSegments
      profile@absoluteSegments <- absoluteSegments
      map@map[[target]] <- profile
    }
  })
  
  setMethod("addSegments", signature(map = "ReferencedCopyNumberMap", target="character", segments="data.frame", isAbsolute="logical"), function(map, target, segments, isAbsolute){ 
    if(!is.null(map@chromosomeSizes)){
      chromosomeSizes = map@chromosomeSizes
      addSegments(map = map, target = target, segments = segments, chromosomeSizes = chromosomeSizes, isAbsolute = isAbsolute)
    } else {
      stop("Add chromosomeSizes before adding segments") # TODO: Bad design: should not require chromosomeSizes to be added beforehand - need verification of chromosomeSizes when creating object
    }
  })
  
  setMethod("solveProbes", signature(map = "ReferencedCopyNumberMap", target="character"), function(map, target){ 
    # TODO: Either both segments and ratio are chromosomal or absolute, but not different - add support for absolute units
    segments <- map@map[[target]]@chromosomalSegments
    ratio <-  map@map[[target]]@chromosomalRatio
      
    if(is.na(segments) || is.na(ratio)){
      stop("Segments and ratio must be set before solving probes")
    } 
    
    probes.df <- do.call(rbind, sapply(seq(1, nrow(segments)), function(segments.index){
      probes.chrom <- segments[segments.index, 1]
      probes.start <- which(ratio[[2]] == segments[segments.index, 2])
      probes.end <- which(ratio[[3]] == segments[segments.index, 3])
      probes.row <- data.frame(chrom=probes.chrom, start=probes.start, end=probes.end, stringsAsFactors = FALSE)
      return(probes.row)
    }))
    
    map@map[[target]]@probes <- probes.df
  })
  
  setMethod("retrieveTargetSampleList", signature(map = "ReferencedCopyNumberMap"), function(map){ 
    return(ls(map@map))
  })
  
  setMethod("retrieveTargetSampleList", signature(map = "ReferencedCopyNumberMap", class = "character"), function(map, class){ 
    targetSampleList <- retrieveTargetSampleList(map) 
    matches <- c()
    if(length(class) > 1){
      matches <- unlist(lapply(class, function(single_class){
        targetSampleList[str_detect(targetSampleList, paste0(".[", single_class, "]"))]
      }))
    } else {
      matches <- targetSampleList[str_detect(targetSampleList, paste0(".[", class, "]"))]
    }
    return(matches)
  })
  
  setMethod("retrieveTargetSampleListSegments", signature(map = "ReferencedCopyNumberMap", class = "character", getAbsolute="logical"), function(map, class, getAbsolute){ 
    targetSampleList <- retrieveTargetSampleList(map, class)
    segmentList <- lapply(targetSampleList, function(sample){
      segments <- NA
      # TODO: Warn on missing/empty segments
      if(getAbsolute == TRUE){
        segments <- map@map[[sample]]@absoluteSegments
      } else {
        segments <- map@map[[sample]]@chromosomalSegments
      }
        
      return(segments)
    })
    targetSampleSegments <- do.call(rbind, segmentList)
  })
  
}

print("Defining Generics")
defineGenerics()
print("Defining FACETS interface class")
defineFacetsClass()
print("Defining standard data class")
defineStandardClass()
print("Complete")

################################################################################################################################################

testFacetsClasses <- function(){
  referenceHashmap <- new("ReferencedFacetsOutputMap", reference="hN30")
  outputContainer <- new("FacetsOutputContainer", xx=saveSnps, fit=saveFit)
  addEntry(map=referenceHashmap, target="hT40", outputContainer=outputContainer)
  removeEntry(map=referenceHashmap, target="hT40")
  referenceHashmap@map$hT40
  print(referenceHashmap@reference)
  addXX(map=referenceHashmap, target="hT41", xx = saveSnps)
  addFit(map=referenceHashmap, target="hT41", fit = saveFit)
}

testStandardClass <- function(){
  referenceStandardMap <- new("ReferencedCopyNumberMap", reference="hN30")
  profile <- new("CopyNumberProfile", ratio=saveSnps, segments=saveFit)
  addEntry(map=referenceStandardMap, target="hT40", profile=profile)
  addSegments(map=referenceStandardMap, target="hT41", segments=saveFit)
  addRatio(map=referenceStandardMap, target="hT41", ratio=saveSnps)
}