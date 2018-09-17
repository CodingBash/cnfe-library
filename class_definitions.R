

defineGenerics <- function(){
  setGeneric("addEntry", function(map, target, profile) {
    standardGeneric("addEntry")
  })
  setGeneric("addRatio", function(map, target, ratio) {
    standardGeneric("addRatio")
  })
  setGeneric("addSegments", function(map, target, segments) {
    standardGeneric("addSegments")
  })
  setGeneric("removeEntry", function(map, target) {
    standardGeneric("removeEntry")
  })
  setGeneric("retrieveTargetSampleList", function(map) {
    standardGeneric("retrieveTargetSampleList")
  })
  setGeneric("retrieveTargetSampleList", function(map, class) {
    standardGeneric("retrieveTargetSampleList")
  })
  setGeneric("retrieveTargetSampleListSegments", function(map, class) {
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
  
  
  setClass("CopyNumberProfile", representation(segments="data.frame", ratio="data.frame"))
  setClass("ReferencedCopyNumberMap", representation(map="environment", reference="character"),
           prototype(map=new.env(), reference=NA_character_))
  setMethod("addEntry", signature(map = "ReferencedCopyNumberMap", target="character", profile="CopyNumberProfile"), function(map, target, profile){ 
    map@map[[target]] <- profile
  })
  setMethod("addRatio", signature(map = "ReferencedCopyNumberMap", target="character", ratio="data.frame"), function(map, target, ratio){ 
    if(is.null(map@map[[target]])) {
      profile <- new("CopyNumberProfile", ratio=ratio)
      map@map[[target]] <- profile
    } else {
      profile <- map@map[[target]]
      profile@ratio <- ratio
      map@map[[target]] <- profile
    }
  })
  setMethod("addSegments", signature(map = "ReferencedCopyNumberMap", target="character", segments="data.frame"), function(map, target, segments){ 
    if(is.null(map@map[[target]])) {
      profile <- new("CopyNumberProfile", segments=segments)
      map@map[[target]] <- profile
    } else {
      profile <- map@map[[target]]
      profile@segments <- segments
      map@map[[target]] <- profile
    }
  })
  setMethod("removeEntry", signature(map = "ReferencedCopyNumberMap", target="character"), function(map, target){ 
    remove(list=c(target), envir=map@map)
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
  
  setMethod("retrieveTargetSampleListSegments", signature(map = "ReferencedCopyNumberMap", class = "character"), function(map, class){ 
    targetSampleList <- retrieveTargetSampleList(map, class)
    segmentList <- lapply(targetSampleList, function(sample){
      return(map@map[[sample]]@segments)
    })
    targetSampleSegments <- do.call(rbind, segmentList)
  })
  
}

defineGenerics()
defineFacetsClass()
defineStandardClass()

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

testStandardClass(){
  referenceStandardMap <- new("ReferencedCopyNumberMap", reference="hN30")
  profile <- new("CopyNumberProfile", ratio=saveSnps, segments=saveFit)
  addEntry(map=referenceStandardMap, target="hT40", profile=profile)
  addSegments(map=referenceStandardMap, target="hT41", segments=saveFit)
  addRatio(map=referenceStandardMap, target="hT41", ratio=saveSnps)
}