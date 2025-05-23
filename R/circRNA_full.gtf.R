#' Generate Full-Length circRNA GTF Files
#' 
#' This function processes circRNA isoform data to generate GTF files for full-length circRNAs, 
#' incorporating reference sources and formatting data according to GTF specifications.
#' 
#' @title Generate Full-Length circRNA GTF Files
#' @param datapathfile Path to the text file containing sample paths (2 columns: SampleID, SamplePath)
#' @param referencefile Path to the reference isoform dataset file (output from MakeReferenceIsoform)
#' @return Invisible NULL (writes GTF files to specified directories)
#' @importFrom data.table fread fwrite
#' @importFrom dplyr filter select mutate
#' @importFrom plyr mapvalues
#' @export
#' @examples
#' \dontrun{
#' circRNA_full.gtf(
#'   datapathfile = "PRJNA429023/DataPathFile.txt",
#'   referencefile = "ReferenceIsoformFinal.txt"
#' )
#' }
circRNA_full.gtf <- function(datapathfile, referencefile) {
  # Load reference isoform data
  message("Loading reference isoform dataset...")
  ReferenceSet <- fread(referencefile, data.table = FALSE) %>%
    dplyr::select(isoformID, ReferenceSource)
  
  # Load sample path data
  message("Reading sample path configuration...")
  SampleData <- fread(datapathfile, data.table = FALSE) %>%
    dplyr::rename(SampleID = SampleID, SamplePath = SamplePath)
  
  # Validate input files
  if (!file.exists(datapathfile)) stop(sprintf("File not found: %s", datapathfile))
  if (!file.exists(referencefile)) stop(sprintf("File not found: %s", referencefile))
  
  # Process each sample
  for (i in seq_len(nrow(SampleData))) {
    SampleID <- SampleData$SampleID[i]
    SamplePath <- SampleData$SamplePath[i]
    
    message(sprintf("\nProcessing sample %d/%d: %s", i, nrow(SampleData), SampleID))
    
    # Create quantification directory
    dirquant <- file.path(SamplePath, "quant")
    if (!dir.exists(dirquant)) {
      message(sprintf("Creating directory: %s", dirquant))
      dir.create(dirquant, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Path to stout.list
    stout_list_path <- file.path(SamplePath, "vis", "stout.list")
    if (!file.exists(stout_list_path)) {
      message(sprintf("Warning: stout.list not found for %s", SampleID))
      next
    }
    
    # Read and process stout.list
    message("Reading and parsing stout.list...")
    stout.list <- fread(
      stout_list_path, 
      data.table = FALSE, 
      sep = "\t", 
      header = FALSE,
      col.names = c(
        "Image_ID", "bsj", "chr", "start", "end", "total_exp",
        "isoform_number", "isoform_exp", "isoform_length", 
        "isoform_state", "strand", "gene_id", "isoform_cirexon"
      )
    ) %>%
      filter(isoform_state == "Full") %>%
      select(chr, start, end, strand, bsj, isoform_cirexon)
    
    # Skip if no full isoforms found
    if (nrow(stout.list) == 0) {
      message("No full-length isoforms found. Skipping GTF generation.")
      next
    }
    
    # Generate GTF entries
    message("Generating GTF entries...")
    gtf_entries <- list()
    
    for (row_idx in seq_len(nrow(stout.list))) {
      row_data <- stout.list[row_idx, ]
      
      # Parse exon coordinates
      exons <- strsplit(row_data$isoform_cirexon, ",")[[1]] %>%
        lapply(function(x) strsplit(x, "-")[[1]]) %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        setNames(c("exon_start", "exon_end")) %>%
        mutate(
          exon_start = as.numeric(exon_start),
          exon_end = as.numeric(exon_end)
        )
      
      # Create isoform ID
      exon_starts <- paste(exons$exon_start, collapse = ",")
      exon_ends <- paste(exons$exon_end, collapse = ",")
      isoform_id <- paste0("chr", row_data$chr, "|", exon_starts, "|", exon_ends, "|", row_data$strand)
      
      # Map reference source
      reference_source <- plyr::mapvalues(isoform_id,ReferenceSet$isoformID,ReferenceSet$ReferenceSource,warn.missing = FALSE)
      
      # Build GTF attributes
      attributes <- paste0(
        'bsj "', row_data$bsj, '"; ',
        'transcript_id "', isoform_id, '"; ',
        'isoform_state "Full"; ',
        'ReferenceSource "', reference_source, '"'
      )
      
      # Create GTF row
      gtf_row <- data.frame(
        chr = row_data$chr,
        source = "ciri",
        type = "exon",
        start = exons$exon_start,
        end = exons$exon_end,
        score = ".",
        strand = row_data$strand,
        phase = ".",
        attributes = attributes,
        stringsAsFactors = FALSE
      )
      
      gtf_entries[[row_idx]] <- gtf_row
    }
    
    # Combine GTF entries
    gtf_output <- do.call(rbind, gtf_entries)
    
    # Write to GTF file
    output_path <- file.path(dirquant, "circRNA_full.gtf")
    message(sprintf("Writing GTF file to: %s", output_path))
    fwrite(
      gtf_output,
      file = output_path,
      sep = "\t",
      quote = FALSE,
      col.names = FALSE,
      row.names = FALSE
    )
    
    message("Processing complete for current sample.")
  }
  
  message("\nAll samples processed successfully.")
  invisible(NULL)
}