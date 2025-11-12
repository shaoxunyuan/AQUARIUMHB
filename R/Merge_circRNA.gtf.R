#' Merge and Process Circular RNA GTF Annotation Files
#'
#' Combines multiple circRNA GTF annotation files, standardizes column names,
#' processes isoform states, and generates a consolidated dataframe of circRNA isoforms.
#'
#' @param SamplePath Dataframe containing sample paths (output from fread("DataPathFile.txt")).
#' @param output_file Character, path for the output file (default: "final_circRNA.df.txt").
#'
#' @return Dataframe with merged circRNA isoform annotations.
#'
#' @importFrom rtracklayer import
#' @importFrom dplyr distinct arrange select
#' @importFrom data.table fread
#'
#' @examples
#' \dontrun{
#' SamplePath <- loadSamplePathFile()
#' merged_isoforms <- Merge_circRNA.gtf(SamplePath, output_file = "final_circRNA.datatable.txt")
#' }
#'
#' @export
Merge_circRNA.gtf <- function(SamplePath, output_file = "final_circRNA.datatable.txt") {
  
  # Define isoform state levels for ordering
  isoform_state_level <- c("Full", "breakinRef_ref", "breakinRef_gtf", "breakoutRef_gtf", "onlyinRef_ref", "onlyoutRef_gtf")
  
  # Collect all GTF file paths
  filelist <- character()
  for (i in 1:nrow(SamplePath)) {
    files <- list.files(paste0(SamplePath$FullPath[i], "/quant/"), full.names = TRUE)
    filelist <- append(filelist, files)
  }
  
  # Import all GTF files
  gtf_list <- list()
  for (i in seq_along(filelist)) {
    message("Processing (", i, "/", length(filelist), "): ", filelist[i])
    gtf_df <- import(filelist[i]) %>% as.data.frame()
    gtf_list[[i]] <- gtf_df
  }
  
  # Merge all GTF dataframes
  merged_gtf <- do.call(rbind, gtf_list)
  rownames(merged_gtf) <- NULL
  
  # Standardize column names
  colnames(merged_gtf)[grep("seqnames", colnames(merged_gtf))] <- "chr"
  colnames(merged_gtf)[grep("start", colnames(merged_gtf))] <- "exon_start"
  colnames(merged_gtf)[grep("end", colnames(merged_gtf))] <- "exon_end"
  colnames(merged_gtf)[grep("width", colnames(merged_gtf))] <- "exon_width"
  
  # Select and sort key columns
  merged_gtf <- merged_gtf %>%
    select(chr, strand, bsj, transcript_id, isoform_state, ReferenceSource) %>%
    arrange(chr, transcript_id) %>%
    distinct()
  
  # Summarize processing statistics
  message("Total records after merging: ", nrow(merged_gtf))
  message("Unique transcript IDs: ", length(unique(merged_gtf$transcript_id)))
  
  # Split by transcript ID
  gtf_by_transcript <- split(merged_gtf, merged_gtf$transcript_id)
  
  # Process each transcript group
  process_transcript <- function(transcript_df) {
    # Convert isoform_state to ordered factor
    ordered_states <- factor(transcript_df$isoform_state, levels = isoform_state_level, ordered = TRUE)
    
    # Collapse ordered states into a comma-separated string
    transcript_df$isoformstate <- paste(ordered_states, collapse = ",")
    
    # Remove original isoform_state column and deduplicate
    transcript_df %>% select(-isoform_state) %>% distinct()
  }
  
  # Apply processing to all transcript groups
  processed_list <- lapply(gtf_by_transcript, process_transcript)
  final_df <- do.call(rbind, processed_list)
  rownames(final_df) <- NULL
  
  # Write output file
  write.table(final_df, file = output_file, row.names = FALSE, quote = FALSE, sep = "\t")
}