#' Generate GTF File for circRNA Break Isoforms (Optimized & Secure)
#'
#' @param SamplePath Data of input, containing sample information, must have columns: SampleName, FullPath.
#' @param ReferenceSet Referenceset data of all possible full-length isoforms.
#'
#' @import data.table
#' @import GenomicRanges
#' @export
circRNA_break.gtf <- function(SamplePath, ReferenceSet) {
  
  # 防止科学计数法输出 (针对整个 R 环境 session)
  old_scipen <- options(scipen = 999)$scipen
  on.exit(options(scipen = old_scipen)) # 函数结束时自动恢复原始设置
  
  message("Loading annotation files and indexing...")
  if (!exists("Homo_sapiens.GRCh38.94.chr.gtf_exontable")) {
    loadAnnotationFiles()
  }
  
  # 转换为 data.table 并建立索引提高检索速度
  gtf_exontable <- as.data.table(Homo_sapiens.GRCh38.94.chr.gtf_exontable)
  gtf_gr <- GRanges(
    seqnames = gtf_exontable$seqnames,
    ranges = IRanges(start = gtf_exontable$start, end = gtf_exontable$end)
  )

  # 内部函数：补全 Break 异构体 (优化区间匹配)
  break_supplyfrom_gtf <- function(onerow, gtf_gr, isoform_state = "breakinRef_gtf") {
    parts <- unlist(strsplit(onerow$isoform_cirexon, ","))
    idx0 <- which(parts == "0-0")
    
    before0 <- as.integer(strsplit(parts[idx0 - 1], "-")[[1]][2])
    after0 <- as.integer(strsplit(parts[idx0 + 1], "-")[[1]][1])
    
    # 极速区间重叠查找
    query_gr <- GRanges(seqnames = onerow$chr, ranges = IRanges(start = before0, end = after0))
    hits <- findOverlaps(gtf_gr, query_gr, type = "within")
    
    # 构建当前已有外显子表
    existing_exons <- data.table(
      chrom = onerow$chr,
      start = as.numeric(sapply(strsplit(parts[parts != "0-0"], "-"), `[`, 1)),
      end = as.numeric(sapply(strsplit(parts[parts != "0-0"], "-"), `[`, 2))
    )
    
    supply_exons <- as.data.table(gtf_gr[queryHits(hits)])[, .(chrom = as.character(seqnames), start, end)]
    
    # 合并、排序并使用 GenomicRanges::reduce 自动处理重叠和相邻区间
    combined_gr <- reduce(GRanges(c(existing_exons$chrom, supply_exons$chrom),
                                  IRanges(c(existing_exons$start, supply_exons$start),
                                          c(existing_exons$end, supply_exons$end))))
    final_df <- as.data.frame(combined_gr)
    
    exonstart <- paste(final_df$start, collapse = ",")
    exonend <- paste(final_df$end, collapse = ",")
    
    return(data.frame(
      chr=onerow$chr, start=onerow$start, end=onerow$end, strand=onerow$strand,
      bsj=onerow$bsj, 
      isoformID=paste0("chr", onerow$chr, "|", exonstart, "|", exonend, "|", onerow$strand),
      isoform_state=isoform_state, ReferenceSource="gtf", stringsAsFactors = FALSE
    ))
  }

  # 处理样本
  SamplePath <- as.data.table(SamplePath)
  ReferenceSet <- as.data.table(ReferenceSet)

  for (i in 1:nrow(SamplePath)) {
    SampleName <- SamplePath$SampleName[i]
    dirquant <- paste0(SamplePath$FullPath[i], "/quant/")
    if (!dir.exists(dirquant)) dir.create(dirquant, recursive = TRUE)
    
    stout_path <- file.path(SamplePath$FullPath[i], "vis/stout.list")
    if (!file.exists(stout_path)) next
    
    message(paste0("Processing sample: ", SampleName))
    stout.list <- fread(stout_path, header = FALSE)
    setnames(stout.list, c("Image_ID", "bsj", "chr", "start", "end", "total_exp",
                           "isoform_number", "isoform_exp", "isoform_length", 
                           "isoform_state", "strand", "gene_id", "isoform_cirexon"))
    
    type_break <- stout.list[isoform_state == "Break"]
    if(nrow(type_break) == 0) next
    
    # 分类处理
    ref_bsjs <- unique(ReferenceSet$bsj)
    type_breakinRef <- type_break[bsj %in% ref_bsjs]
    type_breakoutRef <- type_break[!(bsj %in% ref_bsjs)]
    
    res_list <- list()
    
    # 1. 处理在参考集中的 (使用正则向量化匹配思路)
    if (nrow(type_breakinRef) > 0) {
      for (j in 1:nrow(type_breakinRef)) {
        row_curr <- type_breakinRef[j, ]
        parts <- unlist(strsplit(row_curr$isoform_cirexon, ","))
        idx0 <- which(parts == "0-0")
        pre <- paste(parts[1:(idx0-1)], collapse=",")
        post <- paste(parts[(idx0+1):length(parts)], collapse=",")
        
        ref_subset <- ReferenceSet[bsj == row_curr$bsj]
        # 匹配逻辑：匹配开头和结尾
        match_vec <- grepl(paste0("^", pre), ref_subset$isoformID) & grepl(paste0(post, "$"), ref_subset$isoformID)
        
        if (any(match_vec)) {
          best <- ref_subset[match_vec][which.max(exon_total_length)]
          res_list[[length(res_list)+1]] <- data.frame(
            chr=row_curr$chr, start=row_curr$start, end=row_curr$end, strand=row_curr$strand,
            bsj=row_curr$bsj, isoformID=best$isoformID, 
            isoform_state="breakinRef_ref", ReferenceSource=best$ReferenceSource
          )
        } else {
          res_list[[length(res_list)+1]] <- break_supplyfrom_gtf(row_curr, gtf_gr, "breakinRef_gtf")
        }
      }
    }
    
    # 2. 处理不在参考集中的
    if (nrow(type_breakoutRef) > 0) {
      out_res <- lapply(1:nrow(type_breakoutRef), function(k) break_supplyfrom_gtf(type_breakoutRef[k,], gtf_gr, "breakoutRef_gtf"))
      res_list <- c(res_list, out_res)
    }
    
    # 汇总结果并向量化生成 GTF 结构
    all_processed <- rbindlist(res_list)
    
    # 向量化展开所有外显子
    # 第一步：把 isoformID 分解为开始列和结束列
    all_processed[, `:=`(starts_str = tstrsplit(isoformID, "\\|")[[2]], 
                         ends_str = tstrsplit(isoformID, "\\|")[[3]])]
    
    # 第二步：展开逗号分隔符 (极速向量化)
    gtf_final <- all_processed[, .(
      start = unlist(strsplit(as.character(starts_str), ",")),
      end = unlist(strsplit(as.character(ends_str), ","))
    ), by = .(chr, bsj, isoformID, isoform_state, ReferenceSource, strand)]
    
    # 构造 GTF 格式列
    gtf_table <- data.table(
      V1 = gsub("chr", "", gtf_final$chr),
      V2 = "ciri",
      V3 = "exon",
      V4 = gtf_final$start,
      V5 = gtf_final$end,
      V6 = ".",
      V7 = gtf_final$strand,
      V8 = ".",
      V9 = paste0('bsj "', gtf_final$bsj, '"; ',
                  'transcript_id "', gtf_final$isoformID, '"; ',
                  'isoform_state "', gtf_final$isoform_state, '"; ',
                  'ReferenceSource "', gtf_final$ReferenceSource, '";')
    )
    
    # 输出文件
    output_file <- paste0(dirquant, "circRNA_break.gtf")
    fwrite(gtf_table, file = output_file, sep = "\t", quote = FALSE, col.names = FALSE)
  }
  
  message("All samples processed successfully!")
}