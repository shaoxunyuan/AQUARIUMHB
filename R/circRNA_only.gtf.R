#' Generate GTF File for circRNA Only Isoforms (Optimized & Industrial)
#'
#' @param SamplePath Data of input, containing sample information, must have columns: SampleName, FullPath.
#' @param ReferenceSet Referenceset data of all possible full-length isoforms.
#'
#' @import data.table
#' @import GenomicRanges
#' @import progress
#' @export
circRNA_only.gtf <- function(SamplePath = samplepath, 
                             ReferenceSet = ReferenceSet) {
  
  # --- 环境与科学计数法设置 ---
  old_scipen <- options(scipen = 999)$scipen
  on.exit(options(scipen = old_scipen))
  
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) stop("Package 'GenomicRanges' is required.")
  if (!requireNamespace("progress", quietly = TRUE)) stop("Package 'progress' is required.")

  # 1. 按照习惯加载注释文件并建立索引
  message("\033[32m>>> Step 1: Loading annotation files and building index...\033[39m")
  loadAnnotationFiles()
  gtf_exontable_dt <- data.table::as.data.table(Homo_sapiens.GRCh38.94.chr.gtf_exontable)
  
  # 构建 GRanges 索引用于极速查找
  gtf_gr <- GenomicRanges::GRanges(
    seqnames = gtf_exontable_dt$seqnames,
    ranges = IRanges::IRanges(start = gtf_exontable_dt$start, end = gtf_exontable_dt$end)
  )

  # 2. 预处理参考集 (循环外执行，提升效率)
  message("\033[32m>>> Step 2: Pre-indexing ReferenceSet...\033[39m")
  ReferenceSet_dt <- data.table::as.data.table(ReferenceSet)
  # 提取用于快速匹配的 BSJ 映射表
  ref_map <- unique(ReferenceSet_dt[, .(bsj, isoformID, ReferenceSource, exon_total_length)])
  # 标记优先级 (Full/Blood 优先)
  ref_map[, is_priority := grepl("Full|Blood", ReferenceSource)]
  data.table::setkey(ref_map, bsj)

  # 内部辅助函数：从 GTF 补全 Only 异构体 (优化版)
  only_supplyfrom_gtf <- function(bsj_id, gtf_gr) {
    # 解析 BSJ 坐标
    parts <- unlist(strsplit(bsj_id, "[:|]"))
    chr_val <- parts[1]; start_val <- as.numeric(parts[2]); end_val <- as.numeric(parts[3])
    
    # 极速检索完全落在 BSJ 内的外显子
    query_gr <- GenomicRanges::GRanges(seqnames = chr_val, ranges = IRanges::IRanges(start = start_val, end = end_val))
    hits <- GenomicRanges::findOverlaps(gtf_gr, query_gr, type = "within")
    hit_indices <- S4Vectors::queryHits(hits)
    
    if (length(hit_indices) > 0) {
      # 提取并合并外显子
      matched_exons <- gtf_gr[hit_indices]
      merged_gr <- GenomicRanges::reduce(matched_exons)
      final_ex_df <- as.data.frame(merged_gr)
      
      # 验证首尾是否对齐 (原代码逻辑)
      if (start_val == final_ex_df$start[1] && end_val == final_ex_df$end[nrow(final_ex_df)]) {
        strand_val <- as.character(runValue(strand(matched_exons))[1])
        isoID <- paste0("chr", chr_val, "|", paste(final_ex_df$start, collapse = ","), "|", 
                        paste(final_ex_df$end, collapse = ","), "|", strand_val)
        
        return(data.frame(
          chr = chr_val, start = start_val, end = end_val, strand = strand_val,
          bsj = bsj_id, isoformID = isoID,
          isoform_state = "onlyoutRef_gtf", ReferenceSource = "gtf", stringsAsFactors = FALSE
        ))
      }
    }
    return(NULL)
  }

  SamplePath_dt <- data.table::as.data.table(SamplePath)
  
  # 3. 样本循环
  for (i in 1:nrow(SamplePath_dt)) {
    SampleName <- SamplePath_dt$SampleName[i]
    message(sprintf("\n\033[34m[*] [%d/%d] Processing Sample: %s\033[39m", i, nrow(SamplePath_dt), SampleName))
    
    dirquant <- paste0(SamplePath_dt$FullPath[i], "quant/")
    if (!dir.exists(dirquant)) dir.create(dirquant, recursive = TRUE)

    # 读取 stout.list 和 ciri.report
    stout_path <- file.path(SamplePath_dt$FullPath[i], "vis/stout.list")
    ciri_path <- file.path(SamplePath_dt$FullPath[i], "full/ciri.report")
    if (!file.exists(stout_path) || !file.exists(ciri_path)) next

    stout_list <- data.table::fread(stout_path, header = FALSE)
    data.table::setnames(stout_list, c("Image_ID", "bsj", "chr", "start", "end", "total_exp",
                                     "isoform_number", "isoform_exp", "isoform_length", 
                                     "isoform_state", "strand", "gene_id", "isoform_cirexon"))
    
    ciri_report <- data.table::fread(ciri_path, header = TRUE)

    # 识别 Only 异构体
    type_only_ids <- setdiff(ciri_report$circRNA_ID, stout_list$bsj)
    if (length(type_only_ids) == 0) next

    # 分类
    type_onlyinRef_ids <- intersect(unique(ref_map$bsj), type_only_ids)
    type_onlyoutRef_ids <- setdiff(type_only_ids, type_onlyinRef_ids)
    
    res_list <- list()

    # 处理在参考集中的 Only
    if (length(type_onlyinRef_ids) > 0) {
      message(paste0("  --- Selecting ", length(type_onlyinRef_ids), " isoforms from ReferenceSet..."))
      # 使用 data.table 极速筛选和排序
      best_ref <- ref_map[bsj %in% type_onlyinRef_ids]
      best_ref <- best_ref[order(bsj, -is_priority, -exon_total_length), .SD[1], by = bsj]
      
      best_ref[, isoform_state := "onlyinRef_ref"]
      res_list[[1]] <- best_ref[, .(chr = sapply(strsplit(bsj, ":"), `[`, 1), 
                                     start = as.numeric(sapply(strsplit(sapply(strsplit(bsj, ":"), `[`, 2), "\\|"), `[`, 1)),
                                     end = as.numeric(sapply(strsplit(bsj, "\\|"), `[`, 2)),
                                     strand = sapply(strsplit(isoformID, "\\|"), `[`, 4),
                                     bsj, isoformID, isoform_state, ReferenceSource)]
    }

    # 处理参考集外的 Only (集成进度条)
    if (length(type_onlyoutRef_ids) > 0) {
      pb <- progress::progress_bar$new(
        format = "  GTF Search [:bar] :percent | ETA: :eta",
        total = length(type_onlyoutRef_ids), clear = FALSE, width = 60
      )
      
      out_res <- lapply(type_onlyoutRef_ids, function(b) {
        pb$tick()
        only_supplyfrom_gtf(b, gtf_gr)
      })
      res_list[[2]] <- data.table::rbindlist(out_res)
    }

    all_only <- data.table::rbindlist(res_list, use.names = TRUE)
    if (nrow(all_only) == 0) next
    
    # 向量化展开外显子行
    all_only[, `:=`(starts_str = data.table::tstrsplit(isoformID, "\\|")[[2]], 
                    ends_str = data.table::tstrsplit(isoformID, "\\|")[[3]])]
    
    gtf_final <- all_only[, .(
      start = unlist(strsplit(as.character(starts_str), ",")),
      end = unlist(strsplit(as.character(ends_str), ","))
    ), by = .(chr, bsj, isoformID, isoform_state, ReferenceSource, strand)]
    
    # 构造 GTF
    gtf_table <- data.table::data.table(
      V1 = gsub("chr", "", gtf_final$chr), V2 = "ciri", V3 = "exon",
      V4 = gtf_final$start, V5 = gtf_final$end, V6 = ".", V7 = gtf_final$strand, V8 = ".",
      V9 = paste0('bsj "', gtf_final$bsj, '"; transcript_id "', gtf_final$isoformID, 
                  '"; isoform_state "', gtf_final$isoform_state, '"; ReferenceSource "', 
                  gtf_final$ReferenceSource, '";')
    )

    data.table::fwrite(gtf_table, file = paste0(dirquant, "circRNA_only.gtf"), sep = "\t", quote = FALSE, col.names = FALSE)
    message(paste0("  \033[32m✔\033[39m Saved ", nrow(gtf_table), " exon records."))
  }
  
  message("\n\033[32m>>> [SUCCESS] All 'only' isoforms processed.\033[39m")
  return(invisible(NULL))
}