#' Generate GTF File for circRNA Only Isoforms
#'
#' This function processes circRNA isoforms that are not present in the visualization data
#' (stout.list) but exist in the CIRI report. These isoforms are classified as "only" isoforms
#' and are processed using either reference data or genome annotation (GTF).
#'
#' @param SamplePath Data of input, containing sample information, must have columns:	SampleName	FullPath (level to sample directory).
#' @param ReferenceSet Referenceset data of all possible full-lenghth isoforms.
#'
#' @return Generates a circRNA_only.gtf file in the quantification directory for each sample.
#' @export
#'
#' @examples
#' \dontrun{
#' circRNA_only.gtf(SamplePath = samplepath, 
#'                  ReferenceSet = ReferenceSet)
#' }
circRNA_only.gtf <- function(SamplePath = samplepath, 
                             ReferenceSet = ReferenceSet) {
  
  # Helper function: Supply missing exon information using GTF annotation
  only_supplyfrom_gtf <- function(bsj){
    # Extract chromosome, start, and end positions from BSJ
    chr <- unlist(strsplit(bsj, split = "[:|]"))[1]
    start <- as.numeric(unlist(strsplit(bsj, split = "[:|]"))[2])
    end <- as.numeric(unlist(strsplit(bsj, split = "[:|]"))[3])
    
    # Retrieve exon information from GTF annotation
    supply_gtf_exon <- gtf_exontable[gtf_exontable$seqnames == chr & 
                                     gtf_exontable$start >= as.numeric(start) & 
                                     gtf_exontable$end <= as.numeric(end), ]
    
    if(nrow(supply_gtf_exon) > 0){
      # Format exon information
      supply_exon <- supply_gtf_exon[, c("seqnames", "start", "end")]
      names(supply_exon) <- c("chrom", "exonStart", "exonEnd")
      supply_exon <- SimRVSequences::combine_exons(supply_exon)
      supply_exon <- as.data.frame(supply_exon)
      
      # Generate isoform ID
      isoformID <- paste0("chr", chr, "|",
                         paste(supply_exon$exonStart, collapse = ","), "|",
                         paste(supply_exon$exonEnd, collapse = ","), "|", 
                         unique(supply_gtf_exon$strand))
      
      # Create results data frame
      results <- data.frame(chr = chr, start = start, end = end, 
                           strand = unique(supply_gtf_exon$strand), 
                           bsj = bsj, isoformID = isoformID,
                           isoform_state = "onlyoutRef_gtf", 
                           ReferenceSource = "gtf")
      
      # Validate that start and end positions match the exon structure
      if (start == supply_exon$exonStart[1] & end == supply_exon$exonEnd[nrow(supply_exon)]) {
        return(results)
      } else {
        return(data.frame())
      }
    } else {
      return(data.frame())
    }
  }
  
  # Helper function: Convert exon table to GTF format
  exontable_to_gtf <- function(exon_df){
    # Initialize list to store GTF entries
    gtf.list <- list()
    
    # Process each row of the exon data frame
    for(index in 1:nrow(exon_df)){
      onerow <- exon_df[index, ]
      chr <- onerow$chr
      ciri <- "ciri"
      type <- "exon" 
      start <- onerow$start
      end <- onerow$end
      attr1 <- "."
      strand <- onerow$strand
      attr2 <- "."
      attr <- paste0('bsj "', onerow$bsj, '"; ', 
                    'transcript_id "', onerow$isoformID, '"; ',
                    'isoform_state "', onerow$isoform_state, '"; ',
                    'ReferenceSource "', onerow$ReferenceSource, '"; ')
      
      # Create GTF entry
      gtfresults <- data.frame(chr = chr, source = "ciri", type = "exon", 
                              start = start, end = end, attr1 = ".", 
                              strand = strand, attr2 = ".", attr = attr)
      gtf.list[[index]] <- gtfresults
    }
    
    # Combine all GTF entries
    gtf <- do.call(rbind, gtf.list)
    return(gtf)
  }
  
  # Load required annotation files
  message("Loading annotation files...")
  loadAnnotationFiles()
  gtf_exontable <- Homo_sapiens.GRCh38.94.chr.gtf_exontable
  
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
    
    # Create quantification directory if not exists
    dirquant <- paste0(SamplePath$FullPath[i], "quant/")
    if (!dir.exists(dirquant)) {
      message(paste0("Creating directory: ", dirquant))
      dir.create(dirquant, recursive = TRUE)
    }

    # Read visualization data
    stout.list.path <- file.path(SamplePath$FullPath[i], "vis/stout.list")
    message(paste0("Reading stout.list for sample ", SampleName))
    stout.list <- data.table::fread(stout.list.path, data.table = FALSE, sep = "\t", header = FALSE)
    colnames(stout.list) <- c("Image_ID", "bsj", "chr", "start", "end", "total_exp",
                             "isoform_number", "isoform_exp", "isoform_length", 
                             "isoform_state", "strand", "gene_id", "isoform_cirexon") 
    stout.list$chr <- gsub("chr", "", stout.list$chr)
    
    # Read CIRI report
    ciri.report.path <- file.path(SamplePath$FullPath[i], "full/ciri.report")
    message(paste0("Reading CIRI report for sample ", SampleName))
    ciri.report <- data.table::fread(ciri.report.path, data.table = FALSE, sep = "\t", header = TRUE)

    # Get different types of BSJs
    message("Classifying circRNA isoforms...")
    type_only <- setdiff(ciri.report$circRNA_ID, stout.list$bsj)
    type_onlyinRef <- intersect(ReferenceSet$bsj, type_only)
    type_onlyoutRef <- setdiff(type_only, ReferenceSet$bsj)
    
    message(paste0("Found ", length(type_only), " 'only' isoforms: ", 
                  length(type_onlyinRef), " in reference, ", 
                  length(type_onlyoutRef), " not in reference"))

    # Process isoforms from reference data
    message("Processing isoforms from reference data...")
    type_onlyinRef.df <- ReferenceSet[ReferenceSet$bsj %in% type_onlyinRef, ]
    type_onlyinRef.df.list <- split(type_onlyinRef.df, type_onlyinRef.df$bsj)
    type_onlyinRef.list <- list()
    
    # Process each BSJ in reference data (show progress every 100 items)
    n_ref_bsj <- length(type_onlyinRef.df.list)
    for (index in 1:n_ref_bsj) {
      if (index %% 100 == 0 || index == n_ref_bsj) {
        percent <- round(index / n_ref_bsj * 100, 1)
        message(paste0("Sample ", SampleName, ": Processed ", index, "/", n_ref_bsj, 
                      " reference BSJs (", percent, "%)"))
      }
      
      onelist <- type_onlyinRef.df.list[[index]]
      
      # Prioritize Full/Blood isoforms, then longest isoform
      if (length(grep("Full|Blood", onelist$ReferenceType)) > 0) {
        select_isoform_for_bsj <- onelist[grep("Full|Blood", onelist$ReferenceSource), ]
        select_isoform_for_bsj <- select_isoform_for_bsj[which.max(select_isoform_for_bsj$exon_total_length), ]
      } else {
        select_isoform_for_bsj <- onelist[which.max(onelist$exon_total_length), ]
      }
      
      type_onlyinRef.list[[index]] <- select_isoform_for_bsj
    }
    
    # Combine results
    type_onlyinRef.df <- do.call(rbind, type_onlyinRef.list)
    type_onlyinRef.df$isoform_state <- "onlyinRef_ref"
    type_onlyinRef.df <- type_onlyinRef.df[, c("chr", "start", "end", "strand", "bsj", "isoformID", "isoform_state", "ReferenceSource")]

    # Process isoforms from GTF annotation
    message("Processing isoforms from GTF annotation...")
    if (length(type_onlyoutRef) > 0) {
      type_onlyoutRef.list <- list()
      n_non_ref_bsj <- length(type_onlyoutRef)
      
      # Process each BSJ using GTF annotation (show progress every 100 items)
      for (index in 1:n_non_ref_bsj) {
        if (index %% 100 == 0 || index == n_non_ref_bsj) {
          percent <- round(index / n_non_ref_bsj * 100, 1)
          message(paste0("Sample ", SampleName, ": Processed ", index, "/", n_non_ref_bsj, 
                        " non-reference BSJs (", percent, "%)"))
        }
        
        bsj <- type_onlyoutRef[index]
        results <- only_supplyfrom_gtf(bsj)
        type_onlyoutRef.list[[index]] <- results
      }
      
      if (length(type_onlyoutRef.list) > 0) {
        type_onlyoutRef.df <- do.call(rbind, type_onlyoutRef.list)
        type_onlyoutRef.df <- type_onlyoutRef.df[!duplicated(type_onlyoutRef.df), ]
      } else {
        type_onlyoutRef.df <- data.frame()
      }
    } else {
      message("No non-reference 'only' isoforms found.")
      type_onlyoutRef.df <- data.frame()
    }

    # Combine all results
    message("Combining all processed isoforms...")
    type_only.df <- rbind(type_onlyinRef.df, type_onlyoutRef.df)
    
    # Convert exon information to GTF format
    message("Converting exon information to GTF format...")
    type_only_exon.list <- list()
    n_isoforms <- nrow(type_only.df)
    
    # Process each isoform (show progress every 100 items)
    for (index in 1:n_isoforms) {
      if (index %% 100 == 0 || index == n_isoforms) {
        percent <- round(index / n_isoforms * 100, 1)
        message(paste0("Sample ", SampleName, ": Processed exon information for ", index, "/", n_isoforms, 
                      " isoforms (", percent, "%)"))
      }
      
      onerow <- type_only.df[index, ]
      exon_start <- as.numeric(unlist(strsplit(unlist(strsplit(onerow$isoformID, split = "\\|"))[[2]], split = ",")))
      exon_end <- as.numeric(unlist(strsplit(unlist(strsplit(onerow$isoformID, split = "\\|"))[[3]], split = ",")))
      
      exon <- data.frame(start = exon_start, end = exon_end)
      exon <- exon[order(exon$start, decreasing = FALSE), ]
      
      multirow <- onerow[rep(1, nrow(exon)), ]
      multirow$start <- exon$start
      multirow$end <- exon$end
      
      type_only_exon.list[[index]] <- multirow
    }
    
    # Combine all exon information
    type_only_exon <- do.call(rbind, type_only_exon.list)
    rownames(type_only_exon) <- NULL
    
    # Convert to GTF format
    type_only_gtf <- exontable_to_gtf(type_only_exon)
    type_only_gtf$chr <- gsub("chr", "", type_only_gtf$chr)
    
    # Write output to file
    output_file <- paste0(dirquant, "circRNA_only.gtf")
    message(paste0("Writing output to: ", output_file))
    options(scipen = 999)
    write.table(type_only_gtf, file = output_file, sep = "\t", quote = FALSE, 
                col.names = FALSE, append = FALSE, row.names = FALSE)
    
    message(paste0("Completed processing for sample ", SampleName))
  }
  
  message("All samples processed successfully!")
}
