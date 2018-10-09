library(CORE)

core <- function(inputSegments, distrib = "Rparallel", maxmark = 100, nshuffle = 500, njobs = 4){
  # TODO: Does CORE accepted chromosome or absolute coordinates. Or any?
  outputCOREobj <- runCORE(inputSegments, distrib = distrib, maxmark = maxmark, nshuffle = nshuffle, njobs = njobs)
  # TODO: Rescaling of CORE table?
  COREtable <- retrieveCORETable(outputCOREobj)
  # TODO: Reformat into standard output
  return(COREtable)
}


#
# Run CORE analysis
#
runCORE <- function(inputCORESegments, inputCOREBoundaries, distrib="Rparallel", maxmark=10, nshuffle=50, seedme, njobs=4) {
  
  myCOREobj<- NA
  
  if(missing(inputCOREBoundaries)){
    myCOREobj <- CORE(dataIn=inputCORESegments, maxmark=maxmark, nshuffle=0,
                      seedme=seedme)
  } else {
    myCOREobj <- CORE(dataIn=inputCORESegments, boundaries = inputCOREBoundaries, maxmark=maxmark, nshuffle=0,
                      seedme=seedme)
  }
  
  newCOREobj<-CORE(dataIn=myCOREobj,keep=c("maxmark","seedme","boundaries"),
                   nshuffle=nshuffle,distrib=distrib,njobs=njobs)
  return(newCOREobj)  
}

#
# Format the CORE table from the CORE result object
#
retrieveCORETable <- function(COREobj){
  COREtable <- data.frame(COREobj$coreTable)
  # TODO: Consistency of chr at beginning of chromosome
  COREtable$chrom <- paste("chr", COREtable$chrom, sep = "")
  return(COREtable)
}


core_selection <- function(){
  coreTable <- NA # Input coreTable
  
  # . . . selection based on conditional
  
  coreSelection <- NA # Output coreSelection
}
