#' Generate GTF File for circRNA Break Isoforms (Industrial Strength)
#'
#' @import data.table
#' @import GenomicRanges
#' @import progress
#' @export
circRNA_break.gtf <- function(SamplePath = samplepath, 
                              ReferenceSet = ReferenceSet) {
  
  # --- 环境与科学计数法设置 ---
  old_scipen <- options(scipen = 999)$scipen
  on.exit(options(scipen = old_scipen))
  
  # 依赖检查
  if (!requireNamespace("progress", quietly = TRUE)) stop("Package 'progress' is required for industrial bars.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) stop("Package 'GenomicRanges' is required.")

  # 1. 加载注释文件
  message("\033[32m>>> Step 1: Loading annotation files...\033[39m")
  loadAnnotationFiles() 
  gtf_exontable_dt <- data.table::as.data.table(Homo_sapiens.GRCh38.94.chr.gtf_exontable)
  gtf_gr <- GenomicRanges::GRanges(
    seqnames = gtf_exontable_dt$seqnames,
    ranges = IRanges::IRanges(start = gtf_exontable_dt$start, end = gtf_exontable_dt$end)
  )

  # 2. 预处理参考集
  message("\033[32m>>> Step 2: Pre-indexing ReferenceSet...\033[39m")
  ReferenceSet_dt <- data.table::as.data.table(ReferenceSet)
  ref_bsj_list <- unique(ReferenceSet_dt$bsj)

  # 内部补全函数 (保持之前的优化版本)
  break_supplyfrom_gtf <- function(onerow, gtf_gr, isoform_state = "breakinRef_gtf") {
    parts <- unlist(strsplit(as.character(onerow$isoform_cirexon), ","))
    idx0 <- which(parts == "0-0")
    before0 <- as.integer(strsplit(parts[idx0 - 1], "-")[[1]][2])
    after0 <- as.numeric(strsplit(parts[idx0 + 1], "-")[[1]][1])
    
    query_gr <- GenomicRanges::GRanges(seqnames = onerow$chr, ranges = IRanges::IRanges(start = before0, end = after0))
    hits <- GenomicRanges::findOverlaps(gtf_gr, query_gr, type = "within")
    hit_indices <- S4Vectors::queryHits(hits)
    
    existing_exons <- data.table::data.table(
      chrom = onerow$chr,
      start = as.numeric(sapply(strsplit(parts[parts != "0-0"], "-"), `[`, 1)),
      end = as.numeric(sapply(strsplit(parts[parts != "0-0"], "-"), `[`, 2))
    )
    supply_exons <- data.table::as.data.table(gtf_gr[hit_indices])[, .(chrom = as.character(seqnames), start, end)]
    
    combined_gr <- GenomicRanges::reduce(GenomicRanges::GRanges(
      c(existing_exons$chrom, supply_exons$chrom),
      IRanges::IRanges(c(existing_exons$start, supply_exons$start), c(existing_exons$end, supply_exons$end))
    ))
    final_exons <- as.data.frame(combined_gr)
    
    return(data.frame(
      chr = onerow$chr, start = onerow$start, end = onerow$end, strand = onerow$strand,
      bsj = onerow$bsj, 
      isoformID = paste0("chr", onerow$chr, "|", paste(final_exons$start, collapse = ","), "|", 
                         paste(final_exons$end, collapse = ","), "|", onerow$strand),
      isoform_state = isoform_state, ReferenceSource = "gtf", stringsAsFactors = FALSE
    ))
  }

  SamplePath_dt <- data.table::as.data.table(SamplePath)
  
  # 3. 样本循环
  for (i in 1:nrow(SamplePath_dt)) {
    SampleName <- SamplePath_dt$SampleName[i]
    message(sprintf("\n\033[34m[*] Processing Sample [%d/%d]: %s\033[39m", i, nrow(SamplePath_dt), SampleName))
    
    dirquant <- paste0(SamplePath_dt$FullPath[i], "quant/")
    if (!dir.exists(dirquant)) dir.create(dirquant, recursive = TRUE)

    stout_path <- file.path(SamplePath_dt$FullPath[i], "vis/stout.list")
    if (!file.exists(stout_path)) next

    stout.list <- data.table::fread(stout_path, header = FALSE)
    data.table::setnames(stout.list, c("Image_ID", "bsj", "chr", "start", "end", "total_exp",
                                     "isoform_number", "isoform_exp", "isoform_length", 
                                     "isoform_state", "strand", "gene_id", "isoform_cirexon"))
    
    type_break <- stout.list[isoform_state == "Break"]
    if (nrow(type_break) == 0) next

    type_breakinRef <- type_break[bsj %in% ref_bsj_list]
    type_breakoutRef <- type_break[!(bsj %in% ref_bsj_list)]
    
    res_list <- list()

    # --- 工业级进度条设置 ---
    n_ref <- nrow(type_breakinRef)
    if (n_ref > 0) {
      # 定义进度条样式
      pb <- progress::progress_bar$new(
        format = "  Matching Reference [:bar] :percent | ETA: :eta | Rate: :rate",
        total = n_ref, clear = FALSE, width = 80
      )
      
      for (j in 1:n_ref) {
        pb$tick() # 推进进度条
        
        row_curr <- type_breakinRef[j, ]
        parts <- unlist(strsplit(as.character(row_curr$isoform_cirexon), ","))
        idx0 <- which(parts == "0-0")
        pre_str <- paste(parts[1:(idx0-1)], collapse = ",")
        post_str <- paste(parts[(idx0+1):length(parts)], collapse = ",")
        
        ref_subset <- ReferenceSet_dt[bsj == row_curr$bsj]
        match_idx <- grepl(paste0("\\|", pre_str), ref_subset$isoformID) & 
                     grepl(paste0(post_str, "\\|"), ref_subset$isoformID)
        
        if (any(match_idx)) {
          best <- ref_subset[match_idx][which.max(exon_total_length)]

          ref_strand <- sub("^.*\\|([+-])$", "\\1", best$isoformID)
          strand_use <- ifelse(row_curr$strand %in% c("+","-"),
                              row_curr$strand,
                              ref_strand)

          res_list[[length(res_list) + 1]] <- data.frame(
            chr = row_curr$chr,
            start = row_curr$start,
            end = row_curr$end,
            strand = strand_use,
            bsj = row_curr$bsj,
            isoformID = best$isoformID,
            isoform_state = "breakinRef_ref",
            ReferenceSource = best$ReferenceSource,
            stringsAsFactors = FALSE
          )
        }
        else {
          res_list[[length(res_list) + 1]] <- break_supplyfrom_gtf(row_curr, gtf_gr, "breakinRef_gtf")
        }
      }
    }

    # 处理不在参考集中的 Break (通常较快，直接处理)
    if (nrow(type_breakoutRef) > 0) {
      message("  --- Supplying non-reference isoforms from GTF...")
      out_res <- lapply(1:nrow(type_breakoutRef), function(k) break_supplyfrom_gtf(type_breakoutRef[k, ], gtf_gr, "breakoutRef_gtf"))
      res_list <- c(res_list, out_res)
    }

    all_processed <- data.table::rbindlist(res_list)
    if (nrow(all_processed) == 0) next
    
    # 展开 GTF 行
    all_processed[, `:=`(starts_str = data.table::tstrsplit(isoformID, "\\|")[[2]], 
                         ends_str = data.table::tstrsplit(isoformID, "\\|")[[3]])]
    
    gtf_final <- all_processed[, .(
      start = unlist(strsplit(as.character(starts_str), ",")),
      end = unlist(strsplit(as.character(ends_str), ","))
    ), by = .(chr, bsj, isoformID, isoform_state, ReferenceSource, strand)]
    
    gtf_final[, strand := ifelse(strand %in% c("+","-"), strand, ".")]
    gtf_table <- data.table::data.table(
      V1 = gsub("chr", "", gtf_final$chr), V2 = "ciri", V3 = "exon",
      V4 = gtf_final$start, V5 = gtf_final$end, V6 = ".", V7 = gtf_final$strand, V8 = ".",
      V9 = paste0('bsj "', gtf_final$bsj, '"; transcript_id "', gtf_final$isoformID, 
                  '"; isoform_state "', gtf_final$isoform_state, '"; ReferenceSource "', 
                  gtf_final$ReferenceSource, '";')
    )

    data.table::fwrite(gtf_table, file = paste0(dirquant, "circRNA_break.gtf"), sep = "\t", quote = FALSE, col.names = FALSE)
    message(paste0("  \033[32m✔\033[39m Saved ", nrow(gtf_table), " records."))
  }
  
  message("\n\033[32m>>> [SUCCESS] All samples processed.\033[39m")
  return(invisible(NULL))
}