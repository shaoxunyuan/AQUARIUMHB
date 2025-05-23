#' Generate Reference Isoform Data
#' 
#' This function processes circRNA isoform data from multiple sources, 
#' combines them, and generates a comprehensive reference isoform dataset.
#' 
#' @param datapathfile Path to the text file containing sample paths
#' @param outputfile Path to the output file where results will be saved
#' 
#' @return A data frame containing the combined and processed reference isoforms
#' 
#' @importFrom data.table fread setnames
#' @import dplyr
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' MakeReferenceIsoform("DataPathFile.txt", "ReferenceIsoformFinal.txt")
#' }
MakeReferenceIsoform <- function(datapathfile, outputfile) {
  message("Starting reference isoform generation...")
  
  # Load external datasets
  message("Loading external datasets...")
  loadLongIsoformFiles()
  FLcircAS <- as.data.frame(FLcircAS)
  IsoCirc <- as.data.frame(IsoCirc)
  
  # Function to generate isoform IDs from coordinate data
  generate_isoform_id <- function(stout_list_full) {
    message("Generating isoform IDs...")
    
    # Convert coordinate string to start|end format
    convert_coordinates <- function(coord_str) {
      if (is.na(coord_str) || coord_str == "") return("NA|NA")
      
      coords <- strsplit(coord_str, ",")[[1]]
      start_coords <- sapply(strsplit(coords, "-"), function(x) x[1])
      end_coords <- sapply(strsplit(coords, "-"), function(x) x[2])
      start_str <- paste(start_coords, collapse = ",")
      end_str <- paste(end_coords, collapse = ",")
      paste(start_str, end_str, sep = "|")
    }
    
    # Apply coordinate conversion and create isoform IDs
    stout_list_full$isoformID <- paste0(
      "chr", 
      stout_list_full$chr, 
      "|", 
      sapply(stout_list_full$isoform_cirexon, convert_coordinates), 
      "|", 
      stout_list_full$strand
    )
    return(stout_list_full)
  }
  
  # Function to parse transcript string and extract exon information
  parse_transcript <- function(transcript_str) {
    # Split string
    parts <- unlist(strsplit(transcript_str, "\\|"))
    
    # Get exon start and end positions (assume parts[2]=starts, parts[3]=ends)
    starts <- as.numeric(unlist(strsplit(parts[2], ",")))
    ends <- as.numeric(unlist(strsplit(parts[3], ",")))
    
    # Calculate exon count
    exon_count <- length(starts)
    
    # Calculate exon lengths
    exon_length <- ends - starts + 1  # Exon length calculation
    exon_total_length <- sum(exon_length)  # Total exon length
    
    # Return results
    return(list(
      exon_count = exon_count, 
      exon_length = exon_length, 
      exon_total_length = exon_total_length
    ))
  }

  # Function to summarize isoform data by grouping and merging sources
  summarize_isoform <- function(isoform_data) {
    message("Summarizing isoform data...")
    
    isoform_list <- split(isoform_data, isoform_data$isoformID)
    
    result_list <- lapply(isoform_list, function(one_isoform) {
      one_isoform <- one_isoform[order(one_isoform$ReferenceSource, decreasing = FALSE), ]
      reference_sources <- unique(one_isoform$ReferenceSource)
      
      # Safely combine reference sources
      ReferenceSource = paste0(grep("Full", reference_sources, value = TRUE), ",", 
							   grep("FLcircAS", reference_sources, value = TRUE), ",", 
							   grep("IsoCirc", reference_sources, value = TRUE))
      
      one_isoform$ReferenceSource <- ReferenceSource
      return(dplyr::distinct(one_isoform))
    })
    
    return(do.call(rbind, result_list))
  }
  
  # Read and process all stout.list files
  message("Reading and processing input files...")
  
  # Read path file
  DataFilePath <- data.table::fread(datapathfile, data.table = FALSE)
  
  stout.list.all.list <- list()
  
  for (i in seq_along(DataFilePath$SamplePath)) {
    message(sprintf(
      "Processing file %d of %d: %s", 
      i, 
      nrow(DataFilePath), 
      DataFilePath$SamplePath[i]
    ))
    
    # Build full path for stout.list
    stout.list.path <- file.path(DataFilePath$SamplePath[i], "vis/stout.list")
    
    # Read file (note: header=FALSE as per original code)
    stout.list <- data.table::fread(
      stout.list.path, 
      data.table = FALSE, 
      header = FALSE
    )
    
    # Set column names
    colnames(stout.list) <- c(
      "Image_ID", "bsj", "chr", "start", "end", "total_exp",
      "isoform_number", "isoform_exp", "isoform_length", 
      "isoform_state", "strand", "gene_id", "isoform_cirexon"
    )
    
    # Select required columns
    stout.list <- stout.list[, c(
      "chr", "bsj", "start", "end", 
      "isoform_state", "strand", "isoform_cirexon"
    )]
    
    stout.list.all.list[[i]] <- stout.list
  }
  
  # Combine all data
  message("Combining and processing data...")
  stout.list.all <- do.call(rbind, stout.list.all.list)
  rownames(stout.list.all) <- NULL
  
  # Filter for full isoforms and remove duplicates
  stout.list.all <- stout.list.all[stout.list.all$isoform_state == "Full", ]
  stout.list.all <- dplyr::distinct(stout.list.all)
  
  # Generate isoform IDs
  stout.list.all <- generate_isoform_id(stout.list.all)
  
  # Parse transcript information
  message("Parsing transcript information...")
  stout.list.all <- stout.list.all %>%
    rowwise() %>%
    mutate(
      parse_results = list(parse_transcript(isoformID)),
      exon_count = parse_results$exon_count,
      exon_length = paste(parse_results$exon_length, collapse = ","),
      exon_total_length = parse_results$exon_total_length
    ) %>%
    dplyr::select(-parse_results)
  
  # Select final columns and add reference source
  stout.list.all <- stout.list.all[, c(
    "chr", "bsj", "start", "end", "isoformID", "strand", 
    "exon_count", "exon_length", "exon_total_length"
  )]
  stout.list.all$ReferenceSource <- "Full"
  
  # Combine with external datasets
  message("Combining with external datasets...")
  ReferenceIsoform <- rbind(
    FLcircAS[, names(stout.list.all)], 
    IsoCirc[, names(stout.list.all)], 
    stout.list.all
  )
  
  # Count occurrences of each isoform
  message("Counting isoform occurrences...")
  isoform_count <- dplyr::count(ReferenceIsoform, isoformID) %>% arrange(desc(n))
  
  # Split into groups by occurrence count
  isoform_count3 <- isoform_count[isoform_count$n == 3, ]
  isoform_count2 <- isoform_count[isoform_count$n == 2, ]
  isoform_count1 <- isoform_count[isoform_count$n == 1, ]
  
  isoform3 <- ReferenceIsoform[ReferenceIsoform$isoformID %in% isoform_count3$isoformID, ]
  isoform2 <- ReferenceIsoform[ReferenceIsoform$isoformID %in% isoform_count2$isoformID, ]
  isoform1 <- ReferenceIsoform[ReferenceIsoform$isoformID %in% isoform_count1$isoformID, ]
  
  # Summarize each group
  message("Summarizing isoform groups...")
  isoform3_summary <- summarize_isoform(isoform3)
  isoform2_summary <- summarize_isoform(isoform2)
  
  # Combine final results
  ReferenceIsoform_Long_And_Full <- rbind(isoform1, isoform2_summary, isoform3_summary)
  
  # Write to file
  message(sprintf("Writing results to %s...", outputfile))
  write.table(
    ReferenceIsoform_Long_And_Full, 
    file = outputfile, 
    sep = "\t", 
    quote = FALSE, 
    col.names = TRUE, 
    append = FALSE, 
    row.names = FALSE
  )
  
  message("Reference isoform generation completed successfully!")
}