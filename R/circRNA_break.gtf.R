#' Generate GTF File for circRNA Break Isoforms
#'
#' This function processes break isoforms of circRNA from visualization data
#' and generates GTF files containing exon information. Break isoforms are
#' those with internal information partially missing, represented by "0-0"
#' in the exon structure. The function attempts to complete these isoforms
#' using reference data or genome annotation.
#'
#' @param SamplePath Data of input, containing sample information, must have columns:	SampleName	FullPath (level to sample directory).
#' @param ReferenceSet Referenceset data of all possible full-lenghth isoforms.
#'
#' @return Generates a circRNA_break.gtf file in the quantification directory for each sample.
#' @export
#'
#' @examples
#' \dontrun{
#' circRNA_break.gtf(SamplePath = samplepath, 
#'                  ReferenceSet = ReferenceSet)
#' }
circRNA_break.gtf <- function(SamplePath = samplepath, 
                             ReferenceSet = ReferenceSet) {

    # Convert isoformID to exon structure string
    isoformID_to_exon = function(isoformID){
        exon_start <- strsplit(strsplit(isoformID, "[|]")[[1]][2],split = ",")[[1]]
        exon_end <- strsplit(strsplit(isoformID, "[|]")[[1]][3],split = ",")[[1]]
        exon <- paste(paste(exon_start, exon_end, sep = "-"),collapse = ",")
        exon
    }
    
    # Supply missing exon information for break isoform using GTF annotation
    break_supplyfrom_gtf <- function(onerow,isoform_state="breakinRef_gtf"){
                
        # Parse existing exon structure
        isoform_cirexon_ciri <- unlist(strsplit(onerow$isoform_cirexon, ","))
        exon <- strsplit(unlist(strsplit(isoform_cirexon_ciri, ",")), split = "-")
        exon <- data.frame(t(data.frame(exon)));
        rownames(exon) <- NULL
        colnames(exon) <- c("exonStart", "exonEnd")
        exon$chrom <- onerow$chr
        exon <- exon[, c("chrom", "exonStart", "exonEnd")];
        
        # Identify break point
        index0 <- which(isoform_cirexon_ciri == "0-0")
        isoform_cirexon_ciri <- strsplit(isoform_cirexon_ciri, '-')
        before0 <- as.integer(isoform_cirexon_ciri[[index0 - 1]][2]);
        after0 <- as.integer(isoform_cirexon_ciri[[index0 + 1]][1]);
        
        # Retrieve supplementary exons from GTF annotation
        supply_exon <- gtf_exontable[gtf_exontable$seqnames == onerow$chr & 
                                     gtf_exontable$start >= before0 & 
                                     gtf_exontable$end <= after0, c("seqnames", "start", "end")]
        supply_exon <- distinct(supply_exon)
        names(supply_exon) <- c("chrom", "exonStart", "exonEnd")
        
        # Combine existing and supplementary exons
        exon_break_supply <- rbind(exon, supply_exon)
        exon_break_supply$exonStart <- as.numeric(exon_break_supply$exonStart)
        exon_break_supply <- exon_break_supply[!exon_break_supply$exonStart == 0, ]
        exon_break_supply$exonEnd <- as.numeric(exon_break_supply$exonEnd)
        
        # Optimize exon structure
        exon_break_supply <- SimRVSequences::combine_exons(exon_break_supply)
        exon_break_supply$exonStart <- as.numeric(exon_break_supply$exonStart)
        exon_break_supply$exonEnd <- as.numeric(exon_break_supply$exonEnd)
        
        # Generate new isoformID and metadata
        exonstart <- paste(exon_break_supply$exonStart,collapse = ",");
        exonend <- paste(exon_break_supply$exonEnd,collapse = ",");
        chr=onerow$chr;start=onerow$start;end=onerow$end;strand=onerow$strand
        len=sum(as.numeric(exon_break_supply$exonEnd)-as.numeric(exon_break_supply$exonStart)+1)  
        BSJ_ID=paste0("chr",onerow$chr,"|",onerow$start,"|",onerow$end,"|",onerow$strand)
        exon_start=exonstart
        exon_end=exonend
        isoformID <- paste0("chr",chr,"|",exonstart,"|",exonend,"|",strand);
        chr_isoform_cirexon=paste0(onerow$chr,"@",onerow$isoform_cirexon);
        
        # Return results
        results <- data.frame(chr=chr,start=start,end=end,strand=strand,
                              bsj=onerow$bsj,isoformID=isoformID,
                              isoform_state=isoform_state,ReferenceSource="gtf")
        return(results)
    }
    
    # Convert exon table to GTF format
    exontable_to_gtf <- function(exon_df){
        message("Converting exon table to GTF format...")
        gtf.list <- list()
        for(index in 1:nrow(exon_df)){
            onerow <- exon_df[index,]
            chr = onerow$chr
            ciri = "ciri"
            type = "exon" 
            start = onerow$start
            end = onerow$end
            attr1 = "."
            strand = onerow$strand
            attr2 = "."
            attr = paste0('bsj "',onerow$bsj,'"; ','transcript_id "',onerow$isoformID,'"; ',
                          'isoform_state "',onerow$isoform_state,'"; ',
                          'ReferenceSource "',onerow$ReferenceSource,'"; ')
            gtfresults <- data.frame(chr=chr,source="ciri",type="exon",start=start,end=end,
                                     attr1=".",strand=strand,attr2=".",attr=attr)
            gtf.list[[index]] <- gtfresults
        }
        gtf <- do.call(rbind,gtf.list)
        gtf                                           
    }
    
    # Load required annotation files
    message("Loading annotation files...")
    loadAnnotationFiles()
    gtf_exontable = Homo_sapiens.GRCh38.94.chr.gtf_exontable

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
        
        # Filter break isoforms
        message("Filtering break isoforms...")
        type_break <- stout.list[stout.list$isoform_state == "Break", ]

        # Separate break isoforms into those in reference and those not
        message("Separating break isoforms into reference and non-reference groups...")
        type_breakinRef <- type_break[type_break$bsj %in% intersect(type_break$bsj, ReferenceSet$bsj), ]
        ReferenceSet_breakinRef <- ReferenceSet[ReferenceSet$bsj %in% type_breakinRef$bsj, ]
        type_breakoutRef_gtf <- type_break[type_break$bsj %in% setdiff(type_break$bsj, ReferenceSet$bsj), ]

        # Process break isoforms in reference (show progress every 100 items)
        message("Processing break isoforms in reference set...")
        type_breakinRef.list <- list()
        n_ref <- nrow(type_breakinRef)
        for(index in 1:n_ref) {
            if (index %% 100 == 0 || index == n_ref) {
                percent <- round(index / n_ref * 100, 1)
                message(paste0("Sample ", SampleName, ": Processed ", index, "/", n_ref, 
                              " reference break isoforms (", percent, "%)"))
            }
            onerow <- type_breakinRef[index, ]
            isoform_cirexon_ciri <- unlist(strsplit(onerow$isoform_cirexon, ","))
            break_index <- which(isoform_cirexon_ciri == "0-0")
            exon_before_break <- paste(isoform_cirexon_ciri[c(1:(break_index - 1))], collapse = ",")
            exon_after_break <- paste(isoform_cirexon_ciri[c(break_index + 1):length(isoform_cirexon_ciri)], collapse = ",")
            
            # Find matching reference isoforms
            ReferenceSet_bsj <- ReferenceSet_breakinRef[ReferenceSet_breakinRef$bsj == onerow$bsj, ]
            ReferenceSet_bsj <- ReferenceSet_bsj[order(ReferenceSet_bsj$exon_total_length, decreasing = TRUE), ]
            select_isoform_for_bsj <- data.frame()
          
            for(isoformindex in 1:nrow(ReferenceSet_bsj)) {
                select_referece <- ReferenceSet_bsj[isoformindex, ]
                select_referece$chr_isoform_cirexon <- paste0(onerow$chr, "@", onerow$isoform_cirexon)
                
                # Get exon structure from reference isoform
                supply_exon <- isoformID_to_exon(select_referece$isoformID)
                
                # Check if reference exon structure matches before and after break
                is_subset_before <- grepl(paste0("^", exon_before_break), supply_exon, perl = TRUE)  
                is_subset_after <- grepl(paste0(exon_after_break, "$"), supply_exon, perl = TRUE) 
                
                if (is_subset_before & is_subset_after) {
                    select_isoform_for_bsj <- rbind(select_isoform_for_bsj, select_referece)
                }
            }

            # Select best matching reference isoform or use GTF annotation
            if(nrow(select_isoform_for_bsj) > 0) {
                select_isoform_for_bsj$isoform_state <- "breakinRef_ref"
                
                # Prioritize Full/Blood isoforms, then longest isoform
                if(length(grep("Full|Blood", select_isoform_for_bsj$ReferenceSource)) > 0) {
                    select_isoform_for_bsj <- select_isoform_for_bsj[grep("Full|Blood", select_isoform_for_bsj$ReferenceSource), ] 
                    select_isoform_for_bsj <- select_isoform_for_bsj[which.max(select_isoform_for_bsj$exon_total_length), ]
                } else {
                    select_isoform_for_bsj <- select_isoform_for_bsj[which.max(select_isoform_for_bsj$exon_total_length), ]
                }
                select_isoform_for_bsj <- select_isoform_for_bsj[, c("chr", "start", "end", "strand", "bsj", "isoformID",
                                                                       "isoform_state", "ReferenceSource")]
            } else {
                select_isoform_for_bsj <- break_supplyfrom_gtf(onerow, "breakinRef_gtf")
            }
            type_breakinRef.list[[index]] <- select_isoform_for_bsj
        }

        type_breakinRef.df <- do.call(rbind, type_breakinRef.list)

        # Process break isoforms not in reference using GTF annotation (show progress every 100 items)
        message("Processing break isoforms not in reference set...")
        n_non_ref <- nrow(type_breakoutRef_gtf)
        if(n_non_ref > 0) {
            type_breakoutRef.list <- list()
            for(index in 1:n_non_ref) {
                if (index %% 100 == 0 || index == n_non_ref) {
                    percent <- round(index / n_non_ref * 100, 1)
                    message(paste0("Sample ", SampleName, ": Processed ", index, "/", n_non_ref, 
                                  " non-reference break isoforms (", percent, "%)"))
                }
                onerow <- type_breakoutRef_gtf[index, ]
                type_breakoutRef.list[[index]] <- break_supplyfrom_gtf(onerow, "breakoutRef_gtf")
            }
            type_breakoutRef.df <- do.call(rbind, type_breakoutRef.list)
        } else {
            message("No non-reference break isoforms found.")
            type_breakoutRef.df <- data.frame()
        }

        # Combine results
        message("Combining processed break isoforms...")
        type_break.df <- rbind(type_breakinRef.df, type_breakoutRef.df)

        # Process exon information from isoformID
        message("Processing exon information from isoformID...")
        type_break_exon.list <- list()
        n_isoforms <- nrow(type_break.df)
        for(index in 1:n_isoforms) {
            if (index %% 100 == 0 || index == n_isoforms) {
                percent <- round(index / n_isoforms * 100, 1)
                message(paste0("Sample ", SampleName, ": Processed exon information for ", index, "/", n_isoforms, 
                              " isoforms (", percent, "%)"))
            }
            onerow <- type_break.df[index, ]
            exon_start <- strsplit(onerow$isoformID, split = "[|]")[[1]][2]
            exon_end <- strsplit(onerow$isoformID, split = "[|]")[[1]][3]
            exon_start <- unlist(strsplit(exon_start, split = ","))
            exon_end <- unlist(strsplit(exon_end, split = ","))
            
            exon <- data.frame(start = exon_start, end = exon_end) # Exon composition data frame
            exon <- exon[order(exon$start, decreasing = FALSE), ]
            
            multirow <- onerow[rep(1, nrow(exon)), ]
            multirow$start <- exon$start
            multirow$end <- exon$end
            multirow$start <- as.numeric(multirow$start)
            multirow$end <- as.numeric(multirow$end)
            
            type_break_exon.list[[index]] <- multirow
        }

        type_break_exon <- do.call(rbind, type_break_exon.list)
        rownames(type_break_exon) <- NULL

        # Generate GTF format
        message("Generating GTF format...")
        type_break_gtf <- exontable_to_gtf(type_break_exon)
        type_break_gtf$chr = gsub("chr","",type_break_gtf$chr)

        # Write output to file
        output_file <- paste0(dirquant, "circRNA_break.gtf")
        message(paste0("Writing output to: ", output_file))
        write.table(type_break_gtf, file = output_file, sep = "\t", quote = FALSE, col.names = FALSE, append = FALSE, row.names = FALSE)
        
        message(paste0("Completed processing for sample ", SampleName))
    }
    
    message("All samples processed successfully!")
}
