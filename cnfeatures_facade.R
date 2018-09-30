setwd("C:/Users/bbece/Documents/Git-Projects/Git-Research-Projects/cnfe-library")
source("./class_definitions.R")
source("./utility_functions.R")
source("./cn_calculation.R")
source("./cn_tool_adapter.R")
source("./cnprep_processing.R")
source("./segment_slicing.R")

# targetStrategy ("SLICE" | "ISOLATE")
run_pipeline <- function(genome = "hg19", targetStrategy = "SLICE"){
  # Calculate Copy Number Profiles
  facetsCopyNumberResults = calculate_cn_profile_test(tool="FACETS") 
  
  # Transform Copy Number Profiles to Standard Format
  standardCopyNumberMapList <- facets_adapter(facetsCopyNumberResults, retrieveChromosomeSizes("hg19"))
  
  # Run CNprep on Copy Number Profiles
  segtableResults <- runCNprep(standardCopyNumberMapList, parallel = FALSE, mclust_model = "E", minjoin = 0, ntrial = 10)
  
  
  if(targetStrategy == "SLICE"){
    runSlicing(segtableResults)
  } else if (targetStrategy == "ISOLATE"){
    
  } else {
    # TODO: Parameter verification should be put at beginning of method
    stop(paste0("targetStrategy not recognized. Accepted inputs: 'SLICE' | 'ISOLATE'. You entered: ", targetStrategy))
  }
  
  # Isolate segments by structure
  
  # Run CORE on 
}