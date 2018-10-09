setwd("C:/Users/bbece/Documents/Git-Projects/Git-Research-Projects/cnfe-library")
# TODO: Organize imports into files which use them
source("./class_definitions.R")
source("./utility_functions.R")
source("./cn_calculation.R")
source("./cn_tool_adapter.R")
source("./cnprep_processing.R")
source("./segment_slicing.R")
source("./isolate_segments.R")

# targetStrategy ("SLICE" | "ISOLATE")
run_pipeline <- function(genome = "hg19", targetStrategy = "SLICE"){
  # Calculate Copy Number Profiles
  facetsCopyNumberResults = calculate_cn_profile_test(tool="FACETS") 
  
  # Transform Copy Number Profiles to Standard Format
  standardCopyNumberMapList <- facets_adapter(facetsCopyNumberResults, retrieveChromosomeSizes("hg19"))
  
  # Run CNprep on Copy Number Profiles
  segtableResults <- runCNprep(standardCopyNumberMapList, parallel = FALSE, mclust_model = "E", minjoin = 0, ntrial = 10)
  
  # TODO: Standardize target results between slicing or isolation (design improvement)
  amplificationSegments <- NA
  deletionSegments <- NA
  isolatedSegments <- NA
  if(targetStrategy == "SLICE"){
    amplificationSegments <- runSlicing(segtableResults, probes = TRUE, amplification = TRUE)
    deletionSegments <- runSlicing(segtableResults, probes = TRUE, amplification = FALSE)
  } else if (targetStrategy == "ISOLATE"){
    amplification_criteria <- function(segment){
      result <- segment[["marginalprob"]]  < 0.001 & segment[["mediandev"]] > 0
      return(result)
    }
    deletion_criteria <- function(segment){
      result <- segment[["marginalprob"]]  < 0.001 & segment[["mediandev"]] > 0
      return(result)
    }
    amplificationSegments <- runIsolation(segtableResults, criteria_function = amplification_criteria, probes = TRUE)
    deletionSegments <- runIsolation(segtableResults, critieraFunction = deletion_criteria, probes = TRUE)
  } else {
    # TODO: Parameter verification should be put at beginning of method
    stop(paste0("targetStrategy not recognized. Accepted inputs: 'SLICE' | 'ISOLATE'. You entered: ", targetStrategy))
  }
  
  # Isolate segments by structure
  
  # Run CORE on 
}