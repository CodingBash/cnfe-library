#
# Given a genome, generate the chromosome sizes
# TODO: Remove "chr" and "X" and "Y". This will make chromosomal <-> absolute conversion easier
#
generateChromosomeSizes <- function(genome){
  seqlengths(genome) <- seqlengths(Hsapiens)
  
  # Create chromosome vector
  chrom_vec <- c(NA)
  chrom_vec.index <- 1
  for(i in append(seq(1,22, by=1), c("X", "Y"))){
    chrom_vec[chrom_vec.index] <- paste("chr", i, sep = "")  
    chrom_vec.index <- chrom_vec.index + 1
  }
  chromosomeSizes <- data.frame(stringsAsFactors = FALSE)
  
  for(chrom_i in chrom_vec){
    df = data.frame(chrom = chrom_i, size = seqlengths(genome)[chrom_i])
    chromosomeSizes <- rbind(chromosomeSizes, df)
  }
  return(chromosomeSizes)
}

#
# Convert a single row from chromosomal units to absolute units
#
genomicConversionForSingleEntry <- function(chrom, start, end, chromosomeSizes, operator){
  chrom_r <- as.numeric(chrom)
  if (is.na(chrom_r) || length(chrom_r) == 0){
    stop("Chromosome entry was NA or empty")
  } 
  
  total_units <- 0
  if(chrom_r == 1){
    # DO NOTHING since chrom 1 is already in correct unit
  } else if(chrom_r %in% seq(2,22)){
    for(i in seq(1, as.numeric(chrom_r) - 1)){
      total_units <- total_units + chromosomeSizes[paste("chr", i, sep = ""), ]$size
    }  
  } else if (chrom_r == "X" | chrom_r == 23) {
    for(i in seq(1, 22)){
      total_units <- total_units + chromosomeSizes[paste("chr", i, sep = ""), ]$size
    }  
  }  else if (chrom_r == "Y" | chrom_r == 24) {
    for(i in seq(1, 22)){
      total_units <- total_units + chromosomeSizes[paste("chr", i, sep = ""), ]$size
    }  
    total_units <- total_units + chromosomeSizes["chrX", ]$size
  } else {
    stop(paste0("Could not convert for chromosome: ", chrom_r))
  }
  converted_start <- operator(start, total_units)
  converted_end <- operator(end, total_units)
  returnme <- data.frame(start = converted_start, end = converted_end)
  return(returnme)
}

#
# Take a input of multiple segments with the chromosomeSizes, and convert
# segment maploc from chrom.location to absolute.location in bp units
#
genomicConversion <- function(inputSegments, chromosomeSizes, chromosomeToAbsolute){
  operator <- if(chromosomeToAbsolute == TRUE) match.fun(FUN = "+") else match.fun(FUN = "-")

  for(row.index in seq(1, nrow(inputSegments))){
    # Pre-process chromosome input
    chrom <- inputSegments[row.index, ][[1]]
    if(is.character(chrom) == TRUE){
      if(substr(chrom, 1, 3) == "chr"){
        chrom <- substr(chrom, 4, nchar(chrom))  
      } else if (!(chrom %in% c("X", "Y"))){
        chrom <- as.numeric(chrom)  
      }
    }
    
    # Convert row
    convertedRow <- genomicConversionForSingleEntry(chrom = chrom, start = inputSegments[row.index, ][[2]], end = inputSegments[row.index, ][[3]], chromosomeSizes = chromosomeSizes, operator=operator)
    
    #
    # Update start and end maploc with new absolute location
    #
    inputSegments[row.index, ][[2]] <- convertedRow$start
    inputSegments[row.index, ][[3]] <- convertedRow$end
  }
  return(inputSegments)
}