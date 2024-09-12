# Creating function to generate GRanges object summarising coverage information
## This function takes as input a .tsv file generated from the alignment_stat_gen.sh script
## Specficially it takes a list of coverage values (1 entry per line) for 100bp windows across a specified genomic region
## It outputs GRanges object that can be visualised with GVIZ

coverage_GRange<-function(cov_file, chr = chrom, strt = start, nd = end) {
  # Read in coverage stats file
  cov_scores <- scan(file = cov_file)
  
  # Generate IRanges object to be used for ranges in GRanges object
  cov_window <- seq(strt, nd-1000, 1000)
  cov_width <- rep(1000, length(cov_window))
  cov_irange <- IRanges(c(cov_window), width=c(cov_width))
  
  # Create GRanges object 
  gr  <- GRanges(
    seqnames = chr,
    ranges = cov_irange, # Using IRanges object created above
    strand = "*", # Strand unimportant, using *
    Cov = cov_scores # Creating the metadata field to be plotted
  )
  
  # Return create granges obejct
  return(gr)
}


# Creating function to generate GRanges object summarising mapping quality information
## This function takes as input a .tsv file generated from the alignment_stat_gen.sh script
## Specficially it takes a list of average mapping quality values (1 entry per line) for 10000bp windows across a specified genomic region
## Rather tha showoing average mapping quality, this function checks to see if the average value for a window is above or below 57
## Greater than 57 is assigned 1, below is assigned 0
## It outputs GRanges object that can be visualised with GVIZ

# Creating function to generate GRanges object summarising coverage information
mapq_GRange<-function(mapq_file, chr = chrom, strt = start, nd = end) {
  # Read in coverage stats file
  mapq_scores <- scan(file = mapq_file)
  
  # Generate IRanges object to be used for ranges in GRanges object
  mapq_window <- seq(strt, nd-10000, 10000)
  mapq_width <- rep(10000, length(mapq_window))
  mapq_irange <- IRanges(c(mapq_window), width=c(mapq_width))
  accessability <- ifelse(mapq_scores < 57, 0, 1)
  
  # Create GRanges object 
  gr  <- GRanges(
    seqnames = chr, # All windows are on chr14
    ranges = mapq_irange, # Using IRanges object created above
    strand = "*", # Strand unimportant, using *
    MAPQ = accessability # Creating the metadata field to be plotted
  )
  
  # Return create granges obejct
  return(gr)
}


# Creating function to generate GRanges object summarising LOD scores
## This function takes as input a number of parametric linkage tables to be combined (merlin output), and an associated map file
## This function only works for a single chromomse (provided as the chrom argument)
## It outputs GRanges object that can be visualised with GVIZ

LOD_grange<-function(..., map, chrom) {
  
  # Read provided parametric data tables into a list
  path_list <- list(...)
  
  if (length(path_list) > 1 ){
    
    # Create list containing data frames housing linkage information
    table_list <- NULL
    for (n in 1:length(path_list)){
      table_list[[n]] <- paste("linkage_table", n, sep = "_")
    }
    
    # Read linkage information into data.frames stored in table_list
    for (n in 1:length(path_list)){
      dat <- read_tsv(path_list[[n]])
      table_list[[n]] <- dat
    }
    
    # create variable for number of linkage tables being combined
    table_num <- length(path_list)
    
    # Create data.frame to store and combine linkage information
    cumulative_LOD <- NULL
    cumulative_LOD <- rbind(cumulative_LOD, table_list[[1]][c("CHR","POS","LABEL","MODEL","LOD")])
    cumulative_LOD$LOD <- as.numeric(as.character(cumulative_LOD$LOD))
    cumulative_names <- c("CHR","POS","LABEL","MODEL","LOD_1")
    for (n in 2:table_num){
      cumulative_names[[4+n]] <- paste("LOD", n, sep = "_")
    }
  
    # Read LOD score data into cumulative data.frame
    for (n in 2:table_num) {
      table_list[[n]]$LOD <- as.numeric(as.character(table_list[[n]]$LOD))
      cumulative_LOD <- cbind(cumulative_LOD, table_list[[n]]["LOD"])
    }
    colnames(cumulative_LOD) <- cumulative_names
    
    # This line removes some spurious lines introduced into the data.frame
    cumulative_LOD <- cumulative_LOD[!grepl("linkage", cumulative_LOD$CHR),]
  
    # Convert any values below 0 to NA
    for (l in 1:table_num) {
      cumulative_LOD[[4+l]] <- replace(cumulative_LOD[[4+l]], which(cumulative_LOD[[4+l]] < 0), NA)
    }
  
    # At loci where all LOD scores are positive sum all LOD scores and put value in new "Total" column
    cumulative_LOD <- cumulative_LOD %>% 
      mutate(LOD = dplyr::select(., -1:-table_num) %>% rowSums(na.rm = TRUE))
    cumulative_LOD$LOD[!complete.cases(cumulative_LOD)] <- NA
    
    # Save cumulative_LOD as lod_table
    lod_table <- cumulative_LOD
  } else if (length(path_list) == 1) {
    lod_table <- read_tsv(...)
    lod_table <- lod_table[ ,1:5]
    lod_table$LOD <- replace(lod_table$LOD, which(lod_table$LOD < 0), NA)
  }
  
  # Add additional columns to data.frame for marker information
  lod_table$LABEL2 <- NA
  lod_table$MARKER <- NA
  lod_table$MARKER_POS <- NA
  lod_table$MARKER_COORD <- NA
  for (n in 1:nrow(lod_table)) {
    lod_table$LABEL2[n] <- lod_table$LABEL[n+1] 
  }
  lod_table$LABEL2[nrow(lod_table)] <- 999
    
  # Read in map file for information on SNPs used
  snp_map <- read_tsv(map)
  
  # Using the map file from merlin to identify markers used in the linkage analysis and assign them to corresponding genomic interval
  for (n in 1:nrow(snp_map)) {
    pos <- snp_map$POSITION[n]
    for (r in 1:nrow(lod_table)) {
      if ( pos >= lod_table$LABEL[r] & pos < lod_table$LABEL2[r]) {
        lod_table$MARKER[r] <- snp_map$MARKER[n]
        lod_table$MARKER_POS[r] <- snp_map$POSITION[n]
      }
    }
  }
  
  # Create biomart object containing SNP information to get GRCh38 coordinates for used SNPs
  snpmart = useEnsembl(biomart = "snp", dataset="hsapiens_snp")
  rs_list <- snp_map$MARKER
  snp_pos <- getBM(attributes = c("refsnp_id", "chrom_start"),
                   filters = "snp_filter",
                   values = rs_list,
                   mart = snpmart)
  
  # Add genomic position of markers to cumulative_LOD data.frame
  for (n in 1:nrow(snp_pos)){
    snp_start <- snp_pos$chrom_start[n]
    snp_name <- snp_pos$refsnp_id[n]
    for (r in 1:nrow(lod_table)) {
      if ( snp_name == lod_table$MARKER[r] & !is.na(lod_table$MARKER[r])) {
        lod_table$MARKER_COORD[r] <- snp_start
      }
    }
  }
  
  # Generate data frame to create GRange object, containing LOD scores and corresponding genomic position.
  LOD_grange_df <- data.frame(lod_table$LOD, lod_table$MARKER_COORD)
  colnames(LOD_grange_df) <- c("LOD", "Marker_Coord")
  LOD_grange_df["LOD"][is.na(LOD_grange_df["LOD"])] <- -3
  
  # Remove windows with no associated SNPs
  LOD_grange_df <- na.omit(LOD_grange_df)
  
  # Create objects for input into GRange object
  IRange_start <- LOD_grange_df$Marker_Coord
  LODscores <- LOD_grange_df$LOD
  LOD_IRange <- IRanges(start=IRange_start, width = 1)
  
  gr <- GRanges(
    seqnames = chrom,
    ranges = LOD_IRange,
    strand = "*",
    LOD = LODscores
  )
  
  return(gr)
}



# Creating function to read in data from Merlin linkage analysis to generated a plot summarising multiple iterations of linkage analysis in a single plot.
## This function takes as input a number of parametric linkage tables to be combined (Merlin output).

WGLinkage_plot_prep <-function(file_path) {

  # List all parametric.tbl files
  parametric_list <- list.files(path=file_path, 
                                pattern="-parametric.tbl",
                                full.names=TRUE)
  
  # Read each of the listed files into data frame
  linkage_table <- NULL
  for (f in parametric_list) {
    print(f)
    dat <- read_tsv(f)
    linkage_table <- rbind(linkage_table, dat)
  }
  
  # Sort dataframe based on chromosome
  linkage_table <- linkage_table %>% arrange(CHR)
  return(linkage_table)
  
}

# Creating function to read in data from a single chromosome from Merlin data linkage analysis to generate a plot showing summed LOD scores for a single chromosome.
## This function takes as input any number of linkage analyses, along with a chromosome of interest.

singleCHR_combined_linkage_plot_prep <- function(..., chrom = 1) {
  
  # Generate list of input paths
  path_list <- list(...)
  
  # Initialise list to store iterations of linkage data
  table_list <- NULL
  
  # Populate list with data.frame names, 1 per linkage iteration
  for (n in 1:length(path_list)){
    table_list[[n]] <- paste("linkage_table", n, sep = "_")
  }
  
  # Store linkage data within data.frames stored in list
  for (n in 1:length(path_list)){
    # Read contents of '-parametic.tbl' from supplied directories into data.frames
    parametric_table <- list.files(path=path_list[[n]], 
                                   pattern="-parametric.tbl",
                                   full.names=TRUE)
    
    # Store contents of "-parametric.tbl" in data.frames stored in above made list
    for (f in parametric_table) {
      dat <- read_tsv(f)
      #      print(dat)
      table_list[[n]] <- rbind(table_list[[n]], dat)
    }
    
    # Sort data.frame based on chromosome
    table_list[[n]] <- table_list[[n]] %>% arrange(CHR)
  }
  
  # Initialise data.frame with LOD scores for all iterations
  cumulative_LOD <- NULL
  
  # Store genomic position data and LOD scores from first linkage iteration
  cumulative_LOD <- rbind(cumulative_LOD, 
                          table_list[[1]][c("CHR","POS","LABEL","MODEL","LOD")])
  
  # Convert LOD from char to num to allow plotting
  cumulative_LOD$LOD <- as.numeric(as.character(cumulative_LOD$LOD))
  
  # Generate list to rename final cumulative_LOD data.frame to give LOD columns different names
  cumulative_names <- c("CHR","POS","LABEL","MODEL","LOD")
  for (n in 2:length(path_list)){
    cumulative_names[[4+n]] <- paste("LOD", n, sep = "_")
  }
  
  # Convert all remaining LOD scores to num and add to cumulative_LOD data.frame 
  for (n in 2:length(path_list)) {
    table_list[[n]]$LOD <- as.numeric(as.character(table_list[[n]]$LOD))
    cumulative_LOD <- cbind(cumulative_LOD, table_list[[n]]["LOD"])
  }
  
  # Rename columns using previously made list
  colnames(cumulative_LOD) <- cumulative_names
  
  # Somewhere in this code "linkage_table_1" gets included in data.frame, this removes it
  cumulative_LOD <- cumulative_LOD[!grepl("linkage",cumulative_LOD$CHR),]
  
  # Add all LOD columns together, place output in new "total" column
  cumulative_LOD <-cumulative_LOD %>% 
    mutate(Total = select(., -1:-4) %>% rowSums(na.rm = TRUE))
  
  # Convert CHR column to num, and sort 'cumulative_LOD' numerically (chromosome order)
  cumulative_LOD$CHR <- as.numeric(as.character(cumulative_LOD$CHR))
  cumulative_LOD <- cumulative_LOD %>% arrange(cumulative_LOD$CHR)
  
  # Downsample to specific chromosome
  singleChr <- subset(cumulative_LOD, CHR == chrom)
  
  # Setting LOD scores below 0 to NA
  singleChr$LOD <- replace(singleChr$LOD, which(singleChr$LOD < 0), NA)
  singleChr$LOD_2 <- replace(singleChr$LOD_2, which(singleChr$LOD_2 < 0), NA)
  singleChr$LOD_3 <- replace(singleChr$LOD_3, which(singleChr$LOD_3 < 0), NA)
  singleChr$LOD_4 <- replace(singleChr$LOD_4, which(singleChr$LOD_4 < 0), NA)
  
  # Downsample columns of data table to contain only relevant fields
  keep <- c("CHR","POS","LABEL","MODEL","LOD","LOD_2","LOD_3","LOD_4")
  
  singleChr <- singleChr[keep]
  
  singleChr <- singleChr %>% 
    mutate(Total = select(., -1:-4) %>% rowSums(na.rm = FALSE))
  
  singleChr["Total"][is.na(singleChr["Total"])] <- -3
  
  return(singleChr)
}
           