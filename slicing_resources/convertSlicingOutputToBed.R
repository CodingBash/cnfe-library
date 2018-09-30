# setwd("~/Git-Projects/Git-Research-Projects/CNprep-event-slicing/")
source("helperFunctions.R")

table <- read.table("output/organoidSlices.txt", sep = "\t", header = TRUE)
chromosomeSizes <- readRDS("resources/meta/chromosomeSizes.rds")

ids <- unique(table$profID)

for(id in ids){
  probeBed <- table[table$profID == id,c(6,2,3,11)]
  probeBed[[1]] <- paste0("chr", probeBed[[1]])
  probeBedA <- probeBed[probeBed[[4]] == 1, ]
  probeBedD <- probeBed[probeBed[[4]] == 0, ]
  
  genomeBed <- table[table$profID == id,c(6,9,10,11)]
  genomeBed <- absoluteToChromosomeBPConversion(genomeBed, chromosomeSizes)
  genomeBed[[1]] <- paste0("chr", genomeBed[[1]])
  
  genomeBedA <- genomeBed[genomeBed[[4]] == 1, ]
  genomeBedD <- genomeBed[genomeBed[[4]] == 0, ]
  
  write.table(probeBed, file = paste0("output/bed/", id, "_slicingProbes.bed"), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(probeBedA, file = paste0("output/bed/", id, "_slicingProbesA.bed"), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(probeBedD, file = paste0("output/bed/", id, "_slicingProbesD.bed"), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  
  write.table(genomeBed, file = paste0("output/bed/", id, "_slicingGenome.bed"), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(genomeBedA, file = paste0("output/bed/", id, "_slicingGenomeA.bed"), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(genomeBedD, file = paste0("output/bed/", id, "_slicingGenomeD.bed"), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
}
