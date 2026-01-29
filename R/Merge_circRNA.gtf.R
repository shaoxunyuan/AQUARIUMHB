#' Merge and Process Circular RNA GTF Annotation Files
#'
#' Combines multiple circRNA GTF annotation files, standardizes column names,
#' processes isoform states, and generates a consolidated dataframe of circRNA isoforms.
#'
#' @param SamplePath Dataframe containing sample paths (output from fread("DataPathFile.txt")).
#' @param output_file Character, path for the output file (default: "final_circRNA.df.txt").
#' @param parallel Boolean, whether to use parallel processing (default: TRUE).
#' @param n_cores Integer, number of cores to use for parallel processing (default: NULL, auto-detect).
#'
#' @return Dataframe with merged circRNA isoform annotations.
#'
#' @importFrom rtracklayer import
#' @importFrom dplyr distinct arrange select mutate map map_dbl map_int map_chr
#' @importFrom data.table fread setDT rbindlist
#' @import purrr
#' @importFrom parallel detectCores mclapply
#'
#' @examples
#' \dontrun{
#' SamplePath <- loadSamplePathFile()
#' merged_isoforms <- Merge_circRNA.gtf(SamplePath, output_file = "final_circRNA.datatable.txt")
#' }
#'
#' @export
Merge_circRNA.gtf <- function(SamplePath, output_file = "final_circRNA.datatable.txt", parallel = TRUE, n_cores = NULL) {
  
  # library(rtracklayer);library(dplyr);library(data.table);library(purrr);
  # 定义 isoform_state 的优先级顺序（用于后续排序）
  isoform_state_level <- c("Full", "breakinRef_ref", "breakinRef_gtf", "breakoutRef_gtf", "onlyinRef_ref", "onlyoutRef_gtf")
  
  # 收集所有样本路径下的 GTF 文件路径
  filelist <- unlist(lapply(1:nrow(SamplePath), function(i) {
    list.files(
      paste0(SamplePath$FullPath[i], "/quant/"),
      pattern = "^(circRNA_full|circRNA_break|circRNA_only)\\.gtf$",
      full.names = TRUE
    )
  }))
  
  # 自动检测核心数
  if (is.null(n_cores)) {
    # 确保 parallel 包可用
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("parallel package not available, using sequential processing")
      n_cores <- 1
    } else {
      n_cores <- max(1, parallel::detectCores() - 1)
    }
  }
  
  # 高效读取 GTF 文件
  message("Reading GTF files...")
  if (parallel && n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    gtf_list <- parallel::mclapply(filelist, function(file) {
      tryCatch({
        import(file) %>% as.data.frame()
      }, error = function(e) {
        message("Error reading ", file, ": ", e$message)
        return(NULL)
      })
    }, mc.cores = n_cores)
  } else {
    gtf_list <- lapply(filelist, function(file) {
      tryCatch({
        import(file) %>% as.data.frame()
      }, error = function(e) {
        message("Error reading ", file, ": ", e$message)
        return(NULL)
      })
    })
  }
  
  # 过滤掉 NULL 结果
  gtf_list <- Filter(Negate(is.null), gtf_list)
  
  # 合并所有 GTF 数据框为一个大表（使用 data.table 提高速度）
  message("Merging GTF data...")
  merged_gtf <- rbindlist(gtf_list, fill = TRUE)
  rownames(merged_gtf) <- NULL
  
  # 标准化列名：将 rtracklayer 导入的列名统一为自定义名称
  colnames(merged_gtf)[grep("seqnames", colnames(merged_gtf))] <- "Chr"
  colnames(merged_gtf)[grep("start", colnames(merged_gtf))] <- "ExonStart"
  colnames(merged_gtf)[grep("end", colnames(merged_gtf))] <- "ExonEnd"
  colnames(merged_gtf)[grep("width", colnames(merged_gtf))] <- "ExonWidth"
  
  # 选择关键列，并按染色体和转录本 ID 排序去重
  merged_gtf <- merged_gtf %>%
    select(Chr, strand, bsj, transcript_id, isoform_state, ReferenceSource) %>%
    arrange(Chr, transcript_id) %>%
    distinct()
  
  # 输出合并后的统计信息
  message("Total records after merging: ", nrow(merged_gtf))
  message("Unique transcript IDs: ", length(unique(merged_gtf$transcript_id)))
  
  # 转换为 data.table 以提高分组处理速度
  setDT(merged_gtf)
  
  # 高效处理转录本状态 - 使用 data.table 分组操作
  message("Processing transcript states...")
  results_df <- merged_gtf[, {
    # 对每个转录本组，按优先级排序并拼接状态
    ordered_states <- factor(isoform_state, levels = isoform_state_level, ordered = TRUE)
    unique_states <- unique(ordered_states)
    sorted_states <- sort(unique_states)
    list(
      Chr = unique(Chr)[1],
      Strand = unique(strand)[1],
      bsj = unique(bsj)[1],
      isoformID = unique(transcript_id)[1],
      IsoformState = paste(sorted_states, collapse = ","),
      ReferenceSource = paste(unique(ReferenceSource), collapse = ",")
    )
  }, by = transcript_id]
  
  # 定义向量化的 isoformID 解析函数
  isoform_sequence_info <- function(isoformID) {
    # 向量化处理 isoformID 解析
    parts_list <- strsplit(isoformID, "|", fixed = TRUE)
    
    chr <- sapply(parts_list, function(x) x[1])
    starts_list <- lapply(parts_list, function(x) as.numeric(strsplit(x[2], ",")[[1]]))
    ends_list <- lapply(parts_list, function(x) as.numeric(strsplit(x[3], ",")[[1]]))
    Strand <- sapply(parts_list, function(x) x[4])
    
    # 计算外显子信息
    ExonCount <- sapply(starts_list, length)
    ExonCoordinate <- mapply(function(s, e) paste(paste(s, e, sep = "-"), collapse = ","), starts_list, ends_list)
    exon_lengths_list <- mapply(function(s, e) e - s + 1, starts_list, ends_list, SIMPLIFY = FALSE)
    IsoformLength <- sapply(exon_lengths_list, sum)
    ExonLength <- sapply(exon_lengths_list, function(x) paste(x, collapse = ","))
    
    # 计算转录本起止位置
    transcript_start <- mapply(function(s, e, str) {
      if (str == "+") min(s) else max(e)
    }, starts_list, ends_list, Strand)
    
    transcript_end <- mapply(function(s, e, str) {
      if (str == "+") max(e) else min(s)
    }, starts_list, ends_list, Strand)
    
    # 计算其他信息
    bsj <- paste0(gsub("chr", "", chr), ":", transcript_start, "|", transcript_end)
    BSJ_ID <- paste(chr, transcript_start, transcript_end, Strand, sep = "|")
    GenomicLength <- transcript_end - transcript_start + 1
    
    # 返回结果列表
    list(
      BSJ_ID = BSJ_ID,
      chr = chr,
      circRNAStart = transcript_start,
      circRNAEnd = transcript_end,
      Strand = Strand,
      GenomicLength = GenomicLength,
      IsoformLength = IsoformLength,
      ExonCount = ExonCount,
      ExonCoordinate = ExonCoordinate,
      ExonLength = ExonLength,
      bsj = bsj
    )
  }
  
  # 使用向量化函数处理所有 isoformID
  message("Parsing isoform IDs...")
  seq_info <- isoform_sequence_info(results_df$isoformID)
  
  # 将解析结果添加到数据框
  results_df_with_seqinfo <- results_df[, c("bsj", "isoformID", "Chr", "Strand", "IsoformState", "ReferenceSource")]
  results_df_with_seqinfo$BSJ_ID <- seq_info$BSJ_ID
  results_df_with_seqinfo$circRNAStart <- seq_info$circRNAStart
  results_df_with_seqinfo$circRNAEnd <- seq_info$circRNAEnd
  results_df_with_seqinfo$GenomicLength <- seq_info$GenomicLength
  results_df_with_seqinfo$IsoformLength <- seq_info$IsoformLength
  results_df_with_seqinfo$ExonCount <- seq_info$ExonCount
  results_df_with_seqinfo$ExonCoordinate <- seq_info$ExonCoordinate
  results_df_with_seqinfo$ExonLength <- seq_info$ExonLength

  # 调整最终输出列的顺序
  setcolorder(results_df_with_seqinfo, c("bsj","BSJ_ID","isoformID","Chr","circRNAStart","circRNAEnd","Strand",
                                          "GenomicLength","IsoformLength","ExonCount","ExonCoordinate","ExonLength",
                                          "IsoformState","ReferenceSource"))

  message("Finished Recording CircRNA Isoform Sequence Information.")

  # 转换为 data.table 以提高后续操作速度
  setDT(results_df_with_seqinfo)
  
  # 开始根据 IsoformState 字段推断 IsoformType 类型
  message("Start Recording CircRNA IsoformState and IsoformType Information.")
  
  # 使用 data.table 高效操作添加 IsoformType 列
  results_df_with_seqinfo[, IsoformType := NA_character_]
  results_df_with_seqinfo[grepl("Full", IsoformState), IsoformType := "Full"]
  results_df_with_seqinfo[!grepl("Full", IsoformState) & grepl("break", IsoformState, ignore.case = TRUE), IsoformType := "Break"]
  results_df_with_seqinfo[!grepl("Full", IsoformState) & !grepl("break", IsoformState, ignore.case = TRUE) & grepl("only", IsoformState, ignore.case = TRUE), IsoformType := "BSJOnly"]

  message("Finished Recording CircRNA IsoformState and IsoformType Information.")

  # 开始处理 ReferenceSource 字段，标准化其格式并推断 ReferenceType
  message("Start Recording CircRNA ReferenceSource and ReferenceType Information.")

  # 定义向量化函数：对 ReferenceSource 字符串按特定顺序重组
  ReCombineReferenceSource <- function(x) {
    # 向量化处理
    sapply(x, function(item) {
      if (is.na(item) || item == "") return(item)
      parts <- unlist(strsplit(item, split = ","))
      
      full_part    <- if ("Full" %in% parts) "Full" else character(0)
      flcircas_part <- parts[grep("FLcircAS", parts, fixed = FALSE)]
      isocirc_part  <- parts[grep("IsoCirc", parts, fixed = FALSE)]
      gtf_part     <- if ("gtf" %in% parts) "gtf" else character(0)
      
      all_parts <- c(
        if (length(full_part) > 0) paste(full_part, collapse = ","),
        if (length(flcircas_part) > 0) paste(flcircas_part, collapse = ","),
        if (length(isocirc_part) > 0) paste(isocirc_part, collapse = ","),
        if (length(gtf_part) > 0) paste(gtf_part, collapse = ",")
      )
      
      if (length(all_parts) == 0) return("")
      paste(all_parts, collapse = ";")
    })
  }

  # 应用重组函数到 ReferenceSource 列
  results_df_with_seqinfo[, ReferenceSource := ReCombineReferenceSource(ReferenceSource)]

  # 初始化 ReferenceType 列
  results_df_with_seqinfo[, ReferenceType := NA_character_]

  # 根据 ReferenceSource 内容分类 ReferenceType（使用 data.table 高效操作）
  results_df_with_seqinfo[grepl("Full", ReferenceSource) & grepl("Blood", ReferenceSource), ReferenceType := "Sample/Blood"]
  results_df_with_seqinfo[grepl("Full", ReferenceSource) & !grepl("Blood", ReferenceSource) & (grepl("FLcircAS", ReferenceSource) | grepl("IsoCirc", ReferenceSource)), ReferenceType := "Sample/non-Blood"]
  results_df_with_seqinfo[grepl("\\bFull\\b", ReferenceSource) & !grepl("Blood|FLcircAS|IsoCirc", ReferenceSource), ReferenceType := "Sample"]
  results_df_with_seqinfo[!grepl("Full", ReferenceSource) & grepl("Blood", ReferenceSource), ReferenceType := "Blood"]
  results_df_with_seqinfo[!grepl("Full|Blood", ReferenceSource) & (grepl("FLcircAS", ReferenceSource) | grepl("IsoCirc", ReferenceSource)), ReferenceType := "non-Blood"]
  results_df_with_seqinfo[grepl("\\bgtf\\b", ReferenceSource), ReferenceType := "gtf"]

  message("Finished Recording CircRNA ReferenceSource and ReferenceType Information.")

  # 将最终结果写入指定输出文件（使用 data.table 的 fwrite 提高速度）
  message("Writing output file...")
  data.table::fwrite(results_df_with_seqinfo, file = output_file, row.names = FALSE, quote = FALSE, sep = "\t")

  # 返回结果
  return(as.data.frame(results_df_with_seqinfo))
}


