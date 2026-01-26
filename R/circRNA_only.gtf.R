#' Generate GTF File for circRNA Only Isoforms (Optimized)
#'
#' @param SamplePath Data of input, containing sample information, must have columns: SampleName, FullPath.
#' @param ReferenceSet Referenceset data of all possible full-length isoforms.
#'
#' @import data.table
#' @import GenomicRanges
#' @export
circRNA_only.gtf <- function(SamplePath, ReferenceSet) {
  
  # 设置科学计数法惩罚，并在函数退出时恢复
  old_scipen <- options(scipen = 999)$scipen
  on.exit(options(scipen = old_scipen))
  
  # 1. 预加载与预处理注释文件
  message("Loading annotation files and indexing...")
  if (!exists("Homo_sapiens.GRCh38.94.chr.gtf_exontable")) {
    loadAnnotationFiles()
  }
  
  # 转换为 data.table 提升查询速度
  gtf_exontable <- as.data.table(Homo_sapiens.GRCh38.94.chr.gtf_exontable)
  # 构建 GRanges 对象用于极速区间 overlap 计算
  gtf_gr <- GRanges(
    seqnames = gtf_exontable$seqnames,
    ranges = IRanges(start = gtf_exontable$start, end = gtf_exontable$end)
  )

  # 内部函数：从 GTF 补全只有 BSJ 坐标的异构体 (优化版)
  only_supplyfrom_gtf <- function(bsj, gtf_gr) {
    # 解析 BSJ 坐标
    parts <- unlist(strsplit(bsj, "[:|]"))
    chr_val <- parts[1]
    start_val <- as.numeric(parts[2])
    end_val <- as.numeric(parts[3])
    
    # 查找完全位于 BSJ 范围内的外显子
    query_gr <- GRanges(seqnames = chr_val, ranges = IRanges(start = start_val, end = end_val))
    hits <- findOverlaps(gtf_gr, query_gr, type = "within")
    
    if (length(hits) > 0) {
      # 提取匹配的外显子并合并相邻/重叠区间 (reduce)
      matched_exons <- gtf_gr[queryHits(hits)]
      merged_gr <- reduce(matched_exons)
      final_exons <- as.data.frame(merged_gr)
      
      # 验证首尾坐标是否完全匹配 (对应原代码逻辑)
      if (start_val == final_exons$start[1] && end_val == final_exons$end[nrow(final_exons)]) {
        
        exon_starts <- paste(final_exons$start, collapse = ",")
        exon_ends <- paste(final_exons$end, collapse = ",")
        strand_val <- as.character(runValue(strand(matched_exons))[1])
        
        isoformID <- paste0("chr", chr_val, "|", exon_starts, "|", exon_ends, "|", strand_val)
        
        return(data.frame(
          chr = chr_val, start = start_val, end = end_val, strand = strand_val,
          bsj = bsj, isoformID = isoformID,
          isoform_state = "onlyoutRef_gtf", ReferenceSource = "gtf",
          stringsAsFactors = FALSE
        ))
      }
    }
    return(NULL) # 不匹配则返回 NULL，方便后续 rbindlist 自动过滤
  }

  # 2. 遍历样本处理
  SamplePath <- as.data.table(SamplePath)
  ReferenceSet <- as.data.table(ReferenceSet)

  for (i in 1:nrow(SamplePath)) {
    SampleName <- SamplePath$SampleName[i]
    dirquant <- paste0(SamplePath$FullPath[i], "/quant/")
    if (!dir.exists(dirquant)) dir.create(dirquant, recursive = TRUE)
    
    # 文件路径准备
    stout_path <- file.path(SamplePath$FullPath[i], "vis/stout.list")
    ciri_path <- file.path(SamplePath$FullPath[i], "full/ciri.report")
    
    if (!file.exists(stout_path) || !file.exists(ciri_path)) {
      message(paste0("Skip: ", SampleName, " (Missing input files)"))
      next
    }
    
    message(paste0("Processing sample: ", SampleName))
    
    # 读取数据
    stout_list <- fread(stout_path, header = FALSE)
    setnames(stout_list, c("Image_ID", "bsj", "chr", "start", "end", "total_exp",
                           "isoform_number", "isoform_exp", "isoform_length", 
                           "isoform_state", "strand", "gene_id", "isoform_cirexon"))
    
    ciri_report <- fread(ciri_path, header = TRUE)
    
    # 识别 "only" 异构体 (CIRI report 中有但 stout.list 中没有)
    type_only_ids <- setdiff(ciri_report$circRNA_ID, stout_list$bsj)
    
    if (length(type_only_ids) == 0) {
      message(paste0("No 'only' isoforms for sample ", SampleName))
      next
    }
    
    # 区分参考集内和参考集外
    type_onlyinRef_ids <- intersect(ReferenceSet$bsj, type_only_ids)
    type_onlyoutRef_ids <- setdiff(type_only_ids, ReferenceSet$bsj)
    
    res_list <- list()
    
    # --- 处理参考集内的 Only ---
    if (length(type_onlyinRef_ids) > 0) {
      ref_hits <- ReferenceSet[bsj %dt_in% type_onlyinRef_ids]
      
      # 对每个 BSJ 选择最优异构体 (利用 data.table 的分组排序功能，极快)
      # 优先级：Full/Blood > 长度最大
      ref_hits[, is_priority := grepl("Full|Blood", ReferenceSource)]
      best_ref <- ref_hits[order(bsj, -is_priority, -exon_total_length), 
                           .SD[1], by = bsj]
      
      best_ref[, isoform_state := "onlyinRef_ref"]
      # 统一列名以匹配后续合并
      res_list[[1]] <- best_ref[, .(chr, start, end, strand, bsj, isoformID, isoform_state, ReferenceSource)]
    }
    
    # --- 处理参考集外的 Only (使用 GTF 补全) ---
    if (length(type_onlyoutRef_ids) > 0) {
      message(paste0("Supplying ", length(type_onlyoutRef_ids), " isoforms from GTF..."))
      # 使用 lapply 替代 for 循环
      out_res <- lapply(type_onlyoutRef_ids, function(b) only_supplyfrom_gtf(b, gtf_gr))
      res_list[[2]] <- rbindlist(out_res)
    }
    
    # 合并所有分类结果
    all_only <- rbindlist(res_list, use.names = TRUE)
    
    if (nrow(all_only) == 0) {
      message(paste0("No valid isoforms generated for ", SampleName))
      next
    }
    
    # --- 核心优化：生成 GTF 数据结构 (完全向量化) ---
    # 1. 快速展开外显子
    all_only[, `:=`(starts_str = tstrsplit(isoformID, "\\|")[[2]], 
                    ends_str = tstrsplit(isoformID, "\\|")[[3]])]
    
    gtf_final <- all_only[, .(
      start = unlist(strsplit(as.character(starts_str), ",")),
      end = unlist(strsplit(as.character(ends_str), ","))
    ), by = .(chr, bsj, isoformID, isoform_state, ReferenceSource, strand)]
    
    # 2. 向量化构造 GTF 列
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
    
    # 写出结果 (使用 fwrite，速度远超 write.table)
    output_file <- paste0(dirquant, "circRNA_only.gtf")
    fwrite(gtf_table, file = output_file, sep = "\t", quote = FALSE, col.names = FALSE)
    
    message(paste0("Completed processing for sample ", SampleName))
  }
  
  message("All samples processed successfully!")
}

# 辅助操作符补充（如果环境中没有 %dt_in%）
`%dt_in%` <- function(x, y) x %in% y