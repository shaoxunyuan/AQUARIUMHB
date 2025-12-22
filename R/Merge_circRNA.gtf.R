#' Merge and Process Circular RNA GTF Annotation Files
#'
#' Combines multiple circRNA GTF annotation files, standardizes column names,
#' processes isoform states, and generates a consolidated dataframe of circRNA isoforms.
#'
#' @param SamplePath Dataframe containing sample paths (output from fread("DataPathFile.txt")).
#' @param output_file Character, path for the output file (default: "final_circRNA.df.txt").
#'
#' @return Dataframe with merged circRNA isoform annotations.
#'
#' @importFrom rtracklayer import
#' @importFrom dplyr distinct arrange select
#' @importFrom data.table fread
#' @import purrr
#'
#' @examples
#' \dontrun{
#' SamplePath <- loadSamplePathFile()
#' merged_isoforms <- Merge_circRNA.gtf(SamplePath, output_file = "final_circRNA.datatable.txt")
#' }
#'
#' @export
Merge_circRNA.gtf <- function(SamplePath, output_file = "final_circRNA.datatable.txt") {
  
  # library(rtracklayer);library(dplyr);library(data.table);library(purrr);
  # 定义 isoform_state 的优先级顺序（用于后续排序）
  isoform_state_level <- c("Full", "breakinRef_ref", "breakinRef_gtf", "breakoutRef_gtf", "onlyinRef_ref", "onlyoutRef_gtf")
  
  # 收集所有样本路径下的 GTF 文件路径
  filelist <- character()
  for (i in 1:nrow(SamplePath)) {
	      files <- list.files(
	      paste0(SamplePath$FullPath[i], "/quant/"),
	      pattern = "^(circRNA_full|circRNA_break|circRNA_only)\\.gtf$",
	      full.names = TRUE
      )
  filelist <- append(filelist, files)
  }
  
  # 逐个读取 GTF 文件并转换为 data.frame，存入列表
  gtf_list <- list()
  for (i in seq_along(filelist)) {
    message("Processing (", i, "/", length(filelist), "): ", filelist[i])
    gtf_df <- import(filelist[i]) %>% as.data.frame()
    gtf_list[[i]] <- gtf_df
  }
  
  # 合并所有 GTF 数据框为一个大表
  merged_gtf <- do.call(rbind, gtf_list)
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
  
  # 按 transcript_id 分组，准备后续处理
  gtf_by_transcript <- split(merged_gtf, merged_gtf$transcript_id)
  
  # 定义处理单个转录本组的函数：
  # 将 isoform_state 转为有序因子后拼接成字符串，并移除原始列
  process_transcript <- function(transcript_df) {
    # 转换为有序因子以保留状态优先级
    ordered_states <- factor(transcript_df$isoform_state, levels = isoform_state_level, ordered = TRUE)
    
    # 拼接所有状态为逗号分隔字符串
    transcript_df$isoformstate <- paste(ordered_states, collapse = ",")
    
    # 移除原始 isoform_state 列并去重
    transcript_df %>% select(-isoform_state) %>% distinct()
  }
  
  # 对每个转录本组应用处理函数
  processed_list <- lapply(gtf_by_transcript, process_transcript)
  results_df <- do.call(rbind, processed_list)
  rownames(results_df) <- NULL
  
  # 重命名关键列为更规范的名称
  colnames(results_df)[which(colnames(results_df) == "transcript_id")] <- "isoformID"
  colnames(results_df)[which(colnames(results_df) == "isoformstate")] <- "IsoformState"
  colnames(results_df)[which(colnames(results_df) == "strand")] <- "Strand"
  
  # 定义从 isoformID 解析序列信息的函数
  isoform_sequence_info <- function(isoformID) {
        parts <- strsplit(isoformID, "|", fixed = TRUE)[[1]]
		chr    <- parts[1]
		starts <- as.numeric(strsplit(parts[2], ",")[[1]])
		ends   <- as.numeric(strsplit(parts[3], ",")[[1]])
		Strand <- parts[4]   # 本來就沒逗號，不要 strsplit！
        
        if (length(starts) != length(ends)) {
          stop("外显子起点和终点数量不匹配")
        }
        
        ExonCount        <- length(starts)                     # 外显子数量
        ExonCoordinate   <- paste0(paste0(starts, "-", ends), collapse = ",")
        exon_lengths     <- ends - starts + 1                  # 单个外显子长度（包含首尾）
        total_length     <- sum(exon_lengths)                  # 总外显子长度
        
        # 根据链方向确定转录本起止位置
        if (parts[4] == "+") {
          transcript_start <- min(starts)
          transcript_end   <- max(ends)
        } else if (parts[4] == "-") {
          transcript_start <- max(ends)
          transcript_end   <- min(starts)
        } else {
          stop("链方向必须是'+'或'-'")
        }
        
        bsj     <- paste0(gsub("chr", "", chr), ":", transcript_start, "|", transcript_end)
        BSJ_ID  <- paste0(chr, "|", transcript_start, "|", transcript_end, "|", Strand)
        GenomicLength <- transcript_end - transcript_start + 1
        IsoformLength <- sum(exon_lengths)
        
        # 返回解析结果（list 格式）
        list(
          BSJ_ID         = BSJ_ID,
          isoformID      = isoformID,
          Chr            = chr,
          circRNAStart   = transcript_start,
          circRNAEnd     = transcript_end,
          Strand         = Strand,
          GenomicLength  = GenomicLength,
          IsoformLength  = IsoformLength,
          ExonCount      = length(starts),
          ExonCoordinate = ExonCoordinate,
          ExonLength     = paste(exon_lengths, collapse = ",")
        )
  }
  
  # 使用 purrr::map 对每个 isoformID 调用解析函数，并展开为新列
  results_df <- results_df %>% mutate(isoformID = as.character(isoformID)) 
  results_df_with_seqinfo <- results_df %>%
      mutate(transcript_info = map(isoformID, isoform_sequence_info)) %>%
      mutate(
        circRNAStart               = map_dbl(transcript_info, ~ .x$circRNAStart),
        circRNAEnd                 = map_dbl(transcript_info, ~ .x$circRNAEnd),
        ExonCount          = map_int(transcript_info, ~ .x$ExonCount),
        ExonLength         = map_chr(transcript_info, ~ .x$ExonLength),
        IsoformLength   = map_dbl(transcript_info, ~ .x$IsoformLength),
        Chr                 = map_chr(transcript_info, ~ .x$Chr),
        Strand              = map_chr(transcript_info, ~ .x$Strand),
        GenomicLength      = map_dbl(transcript_info, ~ .x$GenomicLength),
        BSJ_ID              = map_chr(transcript_info, ~ .x$BSJ_ID),
        ExonCoordinate     = map_chr(transcript_info, ~ .x$ExonCoordinate)
      ) %>%
      select(-transcript_info)

  # 调整最终输出列的顺序
  results_df_with_seqinfo = results_df_with_seqinfo[,c("bsj","BSJ_ID","isoformID","Chr","circRNAStart","circRNAEnd","Strand",
                                                       "GenomicLength","IsoformLength","ExonCount","ExonCoordinate","ExonLength",
                                                       "IsoformState","ReferenceSource")]

  message("Finished Recording CircRNA Isoform Sequence Information.")

  
  # 开始根据 IsoformState 字段推断 IsoformType 类型
  message("Start Recording CircRNA IsoformState and IsoformType Information.")
  results_df_with_IsoformType = results_df_with_seqinfo

  results_df_with_IsoformType$IsoformType = NA

  # 步骤1：若包含 "Full"，则类型为 "Full"
  results_df_with_IsoformType$IsoformType[grepl("Full", results_df_with_IsoformType$IsoformState)] <- "Full"

  # 步骤2：不含 "Full" 但含 "break"（忽略大小写），则类型为 "Break"
  mask_break <- !grepl("Full", results_df_with_IsoformType$IsoformState) & grepl("break", results_df_with_IsoformType$IsoformState, ignore.case = TRUE)
  results_df_with_IsoformType$IsoformType[mask_break] <- "Break"

  # 步骤3：不含 "Full" 和 "break" 但含 "only"，则类型为 "BSJOnly"
  mask_only <- !grepl("Full", results_df_with_IsoformType$IsoformState) & 
              !grepl("break", results_df_with_IsoformType$IsoformState, ignore.case = TRUE) & 
              grepl("only", results_df_with_IsoformType$IsoformState, ignore.case = TRUE)
  results_df_with_IsoformType$IsoformType[mask_only] <- "BSJOnly"

  message("Finished Recording CircRNA IsoformState and IsoformType Information.")

  # 开始处理 ReferenceSource 字段，标准化其格式并推断 ReferenceType
  message("Start Recording CircRNA ReferenceSource and ReferenceType Information.")

  # 定义函数：对 ReferenceSource 字符串按特定顺序重组（Full > FLcircAS > IsoCirc > gtf）
  ReCombineReferenceSource <- function(x) {
    if (is.na(x) || x == "") return(x)  # 处理缺失或空值
    parts <- unlist(strsplit(x, split = ","))
    
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
  }

  # 应用重组函数到 ReferenceSource 列
  results_df_with_ReferenceSource = results_df_with_IsoformType
  results_df_with_ReferenceSource$ReferenceSource <- sapply(results_df_with_ReferenceSource$ReferenceSource, ReCombineReferenceSource)

  # 转换为 data.table 以便高效赋值
  setDT(results_df_with_ReferenceSource)

  # 初始化 ReferenceType 列
  results_df_with_ReferenceSource$ReferenceType = NA
  results_df_with_ReferenceSource[, ReferenceType := as.character(ReferenceType)]

  # 根据 ReferenceSource 内容分类 ReferenceType
  # idx1: 同时含 Full 和 Blood → "Sample/Blood"
  idx1 <- grepl("Full", results_df_with_ReferenceSource$ReferenceSource) & grepl("Blood", results_df_with_ReferenceSource$ReferenceSource)
  results_df_with_ReferenceSource[idx1, c( "ReferenceType")] <- list( "Sample/Blood")

  # idx2: 含 Full、不含 Blood，但含 FLcircAS 或 IsoCirc → "Sample/non-Blood"
  idx2 <- grepl("Full", results_df_with_ReferenceSource$ReferenceSource) & !grepl("Blood", results_df_with_ReferenceSource$ReferenceSource) & (grepl("FLcircAS", results_df_with_ReferenceSource$ReferenceSource) | grepl("IsoCirc", results_df_with_ReferenceSource$ReferenceSource))
  results_df_with_ReferenceSource[idx2, c( "ReferenceType")] <- list( "Sample/non-Blood")

  # idx3: 仅含 Full（无 Blood/FLcircAS/IsoCirc）→ "Sample"
  idx3 <- grepl("\\bFull\\b", results_df_with_ReferenceSource$ReferenceSource) & !grepl("Blood|FLcircAS|IsoCirc", results_df_with_ReferenceSource$ReferenceSource)
  results_df_with_ReferenceSource[idx3, c( "ReferenceType")] <- list( "Sample")

  # idx4: 不含 Full，但含 Blood → "Blood"
  idx4 <- !grepl("Full", results_df_with_ReferenceSource$ReferenceSource) & grepl("Blood", results_df_with_ReferenceSource$ReferenceSource)
  results_df_with_ReferenceSource[idx4, c( "ReferenceType")] <- list( "Blood")

  # idx5: 不含 Full/Blood，但含 FLcircAS 或 IsoCirc → "non-Blood"
  idx5 <- !grepl("Full|Blood", results_df_with_ReferenceSource$ReferenceSource) & (grepl("FLcircAS", results_df_with_ReferenceSource$ReferenceSource) | grepl("IsoCirc", results_df_with_ReferenceSource$ReferenceSource))
  results_df_with_ReferenceSource[idx5, c( "ReferenceType")] <- list( "non-Blood")

  # idx6: 仅含 gtf → "gtf"
  idx6 <- grepl("\\bgtf\\b", results_df_with_ReferenceSource$ReferenceSource)
  results_df_with_ReferenceSource[idx6, c( "ReferenceType")] <- list( "gtf")

  message("Finished Recording CircRNA ReferenceSource and ReferenceType Information.")

  
  # 将最终结果写入指定输出文件（制表符分隔，无行名和引号）
  write.table(results_df_with_ReferenceSource, file = output_file, row.names = FALSE, quote = FALSE, sep = "\t")

}
