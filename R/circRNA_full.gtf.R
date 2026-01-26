#' Generate circRNA Full GTF File (Industrial Strength)
#'
#' @param SamplePath Data of input, containing sample information, must have columns: SampleName, FullPath.
#' @param ReferenceSet Referenceset data of all possible full-length isoforms.
#'
#' @import data.table
#' @import progress
#' @export
circRNA_full.gtf <- function(SamplePath = samplepath, 
                             ReferenceSet = ReferenceSet) {
  
  # --- 环境设置与依赖检查 ---
  old_scipen <- options(scipen = 999)$scipen
  on.exit(options(scipen = old_scipen))
  
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
  if (!requireNamespace("progress", quietly = TRUE)) stop("Package 'progress' is required for industrial bars.")

  # --- Step 1: 循环外处理预加载的 ReferenceSet ---
  message("\033[32m>>> Step 1: Pre-indexing ReferenceSet (one-time operation)...\033[39m")
  ReferenceSet_dt <- data.table::as.data.table(ReferenceSet)
  
  # 建立高速 Join 索引
  ref_map <- unique(ReferenceSet_dt[, .(isoformID, ReferenceSource)])
  data.table::setkey(ref_map, isoformID)
  
  SamplePath_dt <- data.table::as.data.table(SamplePath)
  num_samples <- nrow(SamplePath_dt)
  message(sprintf("--- ReferenceSet indexed. Ready to process %d samples.", num_samples))

  # --- Step 2: 样本循环 ---
  for (i in 1:num_samples) {
    SampleName <- SamplePath_dt$SampleName[i]
    # 使用 ANSI 颜色高亮样本标题
    message(sprintf("\n\033[34m[*] [%d/%d] Processing Sample: %s\033[39m", i, num_samples, SampleName))
    
    dirquant <- paste0(SamplePath_dt$FullPath[i], "quant/")
    if (!dir.exists(dirquant)) dir.create(dirquant, recursive = TRUE)
    
    stout.list.path <- file.path(SamplePath_dt$FullPath[i], "vis/stout.list")
    if (!file.exists(stout.list.path)) {
      message("\033[31m  ! Skipping: stout.list not found.\033[39m")
      next
    }
    
    stout.list <- data.table::fread(stout.list.path, header = FALSE)
    data.table::setnames(stout.list, c("Image_ID", "bsj", "chr", "start", "end", "total_exp",
                                     "isoform_number", "isoform_exp", "isoform_length", 
                                     "isoform_state", "strand", "gene_id", "isoform_cirexon"))
    
    type_full <- stout.list[isoform_state == "Full"]
    num_full <- nrow(type_full)
    
    if (num_full == 0) {
      message("  --- No 'Full' isoforms found.")
      next
    }

    # --- 工业级进度条：针对 ID 生成阶段 ---
    # 因为 sapply 循环是该函数在单样本处理中最耗时的部分
    pb <- progress::progress_bar$new(
      format = "  Generating IDs [:bar] :percent | ETA: :eta",
      total = num_full, clear = FALSE, width = 60
    )

    type_full[, `:=`(
      bsj_str = paste0(chr, ":", start, "|", end),
      isoformID = sapply(1:.N, function(j){
        pb$tick() # 推进进度条
        exons <- unlist(strsplit(as.character(isoform_cirexon[j]), ","))
        coords <- data.table::tstrsplit(exons, "-")
        paste0("chr", chr[j], "|", 
               paste(coords[[1]], collapse = ","), "|", 
               paste(coords[[2]], collapse = ","), "|", 
               strand[j])
      })
    )]

    # --- 展开与 Join ---
    message("  --- Expanding coordinates and joining reference...")
    gtf_data <- type_full[, {
      exons <- unlist(strsplit(as.character(isoform_cirexon), ","))
      coords <- data.table::tstrsplit(exons, "-")
      list(
        start = as.numeric(coords[[1]]),
        end = as.numeric(coords[[2]])
      )
    }, by = .(chr, strand, bsj = bsj_str, isoformID)]

    # 极速匹配
    gtf_data <- ref_map[gtf_data, on = "isoformID"]

    # 构造属性列
    gtf_data[, attr := paste0('bsj "', bsj, '"; ',  
                              'transcript_id "', isoformID, '"; ',  
                              'isoform_state "Full"; ',  
                              'ReferenceSource "', ReferenceSource, '";')]

    # 构造最终表格
    final_gtf <- gtf_data[, .(
      chr = chr, source = "ciri", type = "exon",
      start = start, end = end, score = ".",
      strand = strand, frame = ".", attributes = attr
    )]

    # 写出文件
    output_file <- paste0(dirquant, "circRNA_full.gtf")
    data.table::fwrite(final_gtf, file = output_file, sep = "\t", quote = FALSE, 
                       col.names = FALSE, row.names = FALSE)
    
    message(sprintf("  \033[32m✔\033[39m Success: %d exon rows written.", nrow(final_gtf)))
  }
  
  message("\n\033[32m>>> [COMPLETE] All samples processed successfully!\033[39m")
  return(invisible(NULL))
}