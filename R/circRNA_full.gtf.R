#' Generate circRNA Full GTF File
#'
#' This function processes visualization data from circRNA detection tools 
#' to generate a GTF file containing exon information for circular RNAs 
#' with 'Full' isoform state.
#'
#' @param SamplePath Data of input, containing sample information, must have columns:	SampleName	FullPath (level to sample directory).
#' @param ReferenceSet Referenceset data of all possible full-lenghth isoforms.
#'
#' @return Generates a circRNA_full.gtf file in the quantification directory for each sample.
#' @export
#'
#' @examples
#' \dontrun{
#' circRNA_full.gtf(SamplePath = samplepath, 
#'                  ReferenceSet = ReferenceSet)
#' }
circRNA_full.gtf <- function(SamplePath = samplepath, 
                             ReferenceSet = ReferenceSet) {
  # Load required packages
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  if (!requireNamespace("plyr", quietly = TRUE)) {
    stop("Package 'plyr' is required but not installed.")
  }
  
  # Read reference file
  message("Reading reference isoform data...")
  nrow(ReferenceSet)
  
  # Read data path file
  message("Reading data path file...")
  nrow(SamplePath)
  
  # Process each sample
  for (i in 1:nrow(SamplePath)) {
    SampleName <- SamplePath$SampleName[i]
    message(paste0("Processing sample: ", SampleName))
    
    dirquant <- paste0(SamplePath$FullPath[i], "quant/")
    if (!dir.exists(dirquant)) {
      message(paste0("Creating directory: ", dirquant))
      dir.create(dirquant, recursive = TRUE)
    }
    
    # Build path to stout.list file
    stout.list.path <- file.path(SamplePath$FullPath[i], "vis/stout.list")
    
    # Check if stout.list file exists
    if (!file.exists(stout.list.path)) {
      warning(paste0("stout.list file not found for sample ", SampleName, ". Skipping this sample."))
      next
    }
    
    message(paste0("Reading stout.list for sample ", SampleName))
    stout.list <- data.table::fread(stout.list.path, data.table = FALSE, sep = "\t", header = FALSE)
    
    # Rename columns with unique names
    colnames(stout.list) <- c("Image_ID", "bsj", "chr", "start", "end", "total_exp",
                             "isoform_number", "isoform_exp", "isoform_length", 
                             "isoform_state", "strand", "gene_id", "isoform_cirexon")
    
    # Filter for 'Full' isoforms
    message("Filtering for 'Full' isoforms...")
    type_full <- stout.list[stout.list$isoform_state == "Full", 
                           c('chr', 'start', 'end', 'strand', 'bsj', 'isoform_state', 'isoform_cirexon')]
    
    # Check if any 'Full' isoforms found
    if (nrow(type_full) == 0) {
      warning(paste0("No 'Full' isoforms found for sample ", SampleName))
      next
    }
    
    # Process each 'Full' isoform to generate GTF entries
    message(paste0("Processing ", nrow(type_full), " 'Full' isoforms..."))
    gtf_Full.list <- list() 
    
    for(index in 1:nrow(type_full)){
      onerow <- type_full[index,]
      chr <- onerow$chr
      start <- onerow$start
      end <- onerow$end
      strand <- onerow$strand
      isoform_cirexon <- onerow$isoform_cirexon
      
      exon <- strsplit(strsplit(isoform_cirexon, ",")[[1]], "-")
      exon <- data.frame(t(data.frame(exon)))
      rownames(exon) <- NULL
      colnames(exon) <- c("start","end")
      exon$start <- as.numeric(exon$start)
      exon$end <- as.numeric(exon$end) 
      
      exonstart <- paste(exon$start, collapse = ",")
      exonend <- paste(exon$end, collapse = ",")
      isoformID <- paste0("chr", chr, "|", exonstart, "|", exonend, "|", strand)
      bsj <- paste0(chr, ":", start, "|", end)
      
      results <- data.frame(chr = chr, ciri = "ciri", type = "exon", start = exon$start, end = exon$end,
                           attr1 = ".", strand = strand, attr2 = ".", bsj = bsj, isoformID = isoformID)
      gtf_Full.list[[index]] <- results           
    }
    
    # Combine all entries
    type_full.gtf <- do.call(rbind, gtf_Full.list)
    
    # Map reference sources
    message("Mapping reference sources...")
    type_full.gtf$ReferenceSource <- plyr::mapvalues(type_full.gtf$isoformID, 
                                                     ReferenceSet$isoformID, 
                                                     ReferenceSet$ReferenceSource, 
                                                     warn_missing = FALSE)
    
    # Create attribute column
    type_full.gtf$attr <- paste0('bsj "', type_full.gtf$bsj, '"; ',  
                                'transcript_id "', type_full.gtf$isoformID, '"; ',  
                                'isoform_state "', "Full", '"; ',  
                                'ReferenceSource "', type_full.gtf$ReferenceSource, '"; ')
    
    # Reorder columns to GTF format
    type_full.gtf <- type_full.gtf[,c("chr","ciri","type","start","end","attr1","strand","attr2","attr")]
    
    # Write output to file
    output_file <- paste0(dirquant, "circRNA_full.gtf")
    message(paste0("Writing output to: ", output_file))

    convert_sci_in_text <- function(x) {
      m <- gregexpr("[+-]?(?:\\d+(?:\\.\\d+)?|\\.\\d+)[eE][+-]?\\d+", x, perl = TRUE)
      hits <- regmatches(x, m)
    
      repl <- lapply(hits, function(v) {
        if (length(v) == 0) return(v)
        format(as.numeric(v), scientific = FALSE, trim = TRUE)
      })
    
      regmatches(x, m) <- repl
      x
    }
    
    type_full.gtf$attr <- convert_sci_in_text(type_full.gtf$attr)

    
    options(scipen = 999)
    write.table(type_full.gtf, file = output_file, sep = "\t", quote = FALSE, 
                col.names = FALSE, append = FALSE, row.names = FALSE)
    
    message(paste0("Completed processing for sample ", SampleName))
  }
  
  message("All samples processed successfully!")
  return(invisible(NULL))
}
