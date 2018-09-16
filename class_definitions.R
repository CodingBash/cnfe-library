library(hashmap)

setClass("FacetsOutputContainer", representation(xx="data.frame", fit="data.frame"),
         prototype(xx=NULL, fit=NULL))
setClass("ReferencedFacetsOutputMap", representation(map="environment", reference="character"),
         prototype(map=new.env(), reference=NA_character_))
setGeneric("addEntry", function(map, target, outputContainer) {
  standardGeneric("addEntry")
})
setGeneric("removeEntry", function(map, target) {
  standardGeneric("removeEntry")
})

setMethod("addEntry", signature(map = "ReferencedFacetsOutputMap", target="character", outputContainer="FacetsOutputContainer"), function(map, target, outputContainer){ 
    map@map[[target]] <- outputContainer
})

setMethod("removeEntry", signature(map = "ReferencedFacetsOutputMap", target="character"), function(map, target){ 
  remove(list=c(target), envir=map@map)
})

################################################################################################################################################

testClasses <- function(){
  referenceHashmap <- new("ReferencedFacetsOutputMap", reference="hN30")
  outputContainer <- new("FacetsOutputContainer", xx=saveSnps, fit=saveFit)
  addEntry(map=referenceHashmap, target="hT40", outputContainer=outputContainer)
  removeEntry(map=referenceHashmap, target="hT40")
  referenceHashmap@map$hT40
  print(referenceHashmap@reference)
}