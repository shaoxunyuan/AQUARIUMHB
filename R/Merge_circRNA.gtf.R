#' Merge and Process Circular RNA GTF Annotation Files
#'
#' Combines multiple circRNA GTF annotation files, standardizes column names,
#' processes isoform states, and generates a consolidated dataframe of 
#' circRNA isoforms.
#'
#' @param sample_path Dataframe containing sample paths (output from 
#'   fread("DataPathFile.txt")).
#' @param output_file Character, path for the output file (default: 
#'   "final_circRNA.datatable.txt").
#'
#' @return Dataframe with merged circRNA isoform annotations.
#'
#' @importFrom data.table fread fwrite rbindlist setDT unique
#' @importFrom stringi stri_match_first_regex stri_split_fixed
#'
#' @examples
#' \dontrun{
#' sample_path <- loadSamplePathFile()
#' merged_isoforms <- Merge_circRNA.gtf(
#'   sample_path,
#'   output_file = "final_circRNA.datatable.txt"
#' )
#' }
#'
#' @export
Merge_circRNA.gtf <- function(
  sample_path,
  output_file = "final_circRNA.datatable.txt"
) {
  if (!is.data.frame(sample_path)) stop("sample_path must be a dataframe")
  if (!"FullPath" %in% colnames(sample_path)) stop("sample_path must contain a 'FullPath' column")
  if (!is.character(output_file) || length(output_file) != 1) stop("output_file must be a character string")

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  if (!requireNamespace("stringi", quietly = TRUE)) {
    stop("Package 'stringi' is required.")
  }

  isoform_state_level <- c(
    "Full", "breakinRef_ref", "breakinRef_gtf",
    "breakoutRef_gtf", "onlyinRef_ref", "onlyoutRef_gtf"
  )

  # 收集所有样本路径下的 GTF 文件路径（避免 for+append 的 O(n^2) 扩容）
  file_list <- unlist(lapply(sample_path$FullPath, function(p) {
    list.files(
      file.path(p, "quant"),
      pattern = "^(circRNA_full|circRNA_break|circRNA_only)\\.gtf$",
      full.names = TRUE
    )
  }), use.names = FALSE)

  if (!length(file_list)) {
    stop("No circRNA gtf files found under <sample>/quant/")
  }

  read_one_gtf_fast <- function(f) {
    dt <- data.table::fread(
      f,
      sep = "\t",
      header = FALSE,
      quote = "",
      col.names = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "attr"),
      select = c("V1", "V7", "attr"),
      showProgress = FALSE
    )
    data.table::setnames(dt, c("V1", "V7"), c("Chr_raw", "Strand"))

    m_bsj <- stringi::stri_match_first_regex(dt$attr, '\\bbsj "([^"]+)"')
    m_tid <- stringi::stri_match_first_regex(dt$attr, '\\btranscript_id "([^"]+)"')
    m_iso <- stringi::stri_match_first_regex(dt$attr, '\\bisoform_state "([^"]+)"')
    m_ref <- stringi::stri_match_first_regex(dt$attr, '\\bReferenceSource "([^"]+)"')

    dt[, `:=`(
      bsj = m_bsj[, 2],
      isoformID = m_tid[, 2], # transcript_id == isoformID
      isoform_state = m_iso[, 2],
      ReferenceSource = m_ref[, 2]
    )]
    dt[, c("Chr_raw", "attr") := NULL]
    dt
  }

  # 逐个读取（I/O + 解析是大头，保持 message 进度输出）
  dt_list <- vector("list", length(file_list))
  for (i in seq_along(file_list)) {
    message("Processing (", i, "/", length(file_list), "): ", file_list[i])
    dt_list[[i]] <- read_one_gtf_fast(file_list[[i]])
  }

  DT <- data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  rm(dt_list)
  invisible(gc())

  # 只保留关键列并去重（对应原脚本 select + distinct）
  DT <- DT[, .(Strand, bsj, isoformID, isoform_state, ReferenceSource)]
  DT <- data.table::unique(DT)

  message("Total records after merging: ", nrow(DT))
  message("Unique isoform IDs: ", data.table::uniqueN(DT$isoformID))

  # 按 isoformID 聚合拼接 IsoformState（替代 split+lapply）
  DT[, isoform_state := factor(isoform_state, levels = isoform_state_level, ordered = TRUE)]
  DTm <- DT[
    ,
    .(
      Strand = Strand[1],
      bsj = bsj[1],
      ReferenceSource = ReferenceSource[1],
      IsoformState = paste(as.character(sort(unique(isoform_state))), collapse = ",")
    ),
    by = .(isoformID)
  ]

  # 解析 isoformID（向量化拆分外层；内层 exon 坐标仍为 list，但避免 per-row 多次 strsplit/map）
  p <- data.table::tstrsplit(DTm$isoformID, "|", fixed = TRUE)
  chr_from_id <- p[[1]]
  starts_s <- p[[2]]
  ends_s <- p[[3]]
  strand_from_id <- p[[4]]

  starts_l <- stringi::stri_split_fixed(starts_s, ",", simplify = FALSE)
  ends_l <- stringi::stri_split_fixed(ends_s, ",", simplify = FALSE)
  exon_count <- lengths(starts_l)

  min_start <- vapply(starts_l, function(x) min(as.numeric(x)), numeric(1))
  max_end <- vapply(ends_l, function(x) max(as.numeric(x)), numeric(1))

  circRNAStart <- data.table::fifelse(strand_from_id == "+", min_start, max_end)
  circRNAEnd <- data.table::fifelse(strand_from_id == "+", max_end, min_start)

  exon_len_l <- Map(function(st, en) as.numeric(en) - as.numeric(st) + 1, starts_l, ends_l)
  isoform_len <- vapply(exon_len_l, sum, numeric(1))
  exon_len_s <- vapply(exon_len_l, function(x) paste(x, collapse = ","), character(1))

  exon_coord <- vapply(
    Map(function(st, en) paste0(st, "-", en), starts_l, ends_l),
    function(x) paste(x, collapse = ","),
    character(1)
  )

  genomic_len <- circRNAEnd - circRNAStart + 1
  BSJ_ID <- paste0(chr_from_id, "|", circRNAStart, "|", circRNAEnd, "|", strand_from_id)

  DTm[, `:=`(
    BSJ_ID = BSJ_ID,
    Chr = chr_from_id,
    circRNAStart = circRNAStart,
    circRNAEnd = circRNAEnd,
    Strand = strand_from_id,
    GenomicLength = genomic_len,
    IsoformLength = isoform_len,
    ExonCount = exon_count,
    ExonCoordinate = exon_coord,
    ExonLength = exon_len_s
  )]

  message("Finished Recording CircRNA Isoform Sequence Information.")

  # IsoformType（需要包含在最终输出）
  DTm[, IsoformType := data.table::fifelse(
    grepl("Full", IsoformState),
    "Full",
    data.table::fifelse(
      grepl("break", IsoformState, ignore.case = TRUE),
      "Break",
      data.table::fifelse(
        grepl("only", IsoformState, ignore.case = TRUE),
        "BSJOnly",
        NA_character_
      )
    )
  )]

  # ReferenceSource 重组（按唯一值缓存，避免每行 strsplit）
  recombine_reference_source <- function(x) {
    if (is.na(x) || x == "") return(x)
    parts <- unlist(strsplit(x, ",", fixed = TRUE))

    full_part <- if ("Full" %in% parts) "Full" else character(0)
    flcircas_part <- parts[grep("FLcircAS", parts)]
    isocirc_part <- parts[grep("IsoCirc", parts)]
    gtf_part <- if ("gtf" %in% parts) "gtf" else character(0)

    all_parts <- c(
      if (length(full_part)) paste(full_part, collapse = ","),
      if (length(flcircas_part)) paste(flcircas_part, collapse = ","),
      if (length(isocirc_part)) paste(isocirc_part, collapse = ","),
      if (length(gtf_part)) paste(gtf_part, collapse = ",")
    )
    if (!length(all_parts)) return("")
    paste(all_parts, collapse = ";")
  }

  u <- unique(DTm$ReferenceSource)
  mp <- stats::setNames(vapply(u, recombine_reference_source, character(1)), u)
  DTm[, ReferenceSource := mp[ReferenceSource]]

  # ReferenceType 分类（保持原逻辑）
  DTm[, ReferenceType := NA_character_]
  DTm[grepl("Full", ReferenceSource) & grepl("Blood", ReferenceSource), ReferenceType := "Sample/Blood"]
  DTm[
    grepl("Full", ReferenceSource) &
      !grepl("Blood", ReferenceSource) &
      (grepl("FLcircAS", ReferenceSource) | grepl("IsoCirc", ReferenceSource)),
    ReferenceType := "Sample/non-Blood"
  ]
  DTm[
    grepl("\\bFull\\b", ReferenceSource) & !grepl("Blood|FLcircAS|IsoCirc", ReferenceSource),
    ReferenceType := "Sample"
  ]
  DTm[!grepl("Full", ReferenceSource) & grepl("Blood", ReferenceSource), ReferenceType := "Blood"]
  DTm[
    !grepl("Full|Blood", ReferenceSource) &
      (grepl("FLcircAS", ReferenceSource) | grepl("IsoCirc", ReferenceSource)),
    ReferenceType := "non-Blood"
  ]
  DTm[grepl("\\bgtf\\b", ReferenceSource), ReferenceType := "gtf"]

  # 最终列顺序（包含 IsoformType）
  out <- DTm[, .(
    bsj, BSJ_ID, isoformID, Chr, circRNAStart, circRNAEnd,
    Strand, GenomicLength, IsoformLength, ExonCount, ExonCoordinate,
    ExonLength, IsoformState, IsoformType, ReferenceSource, ReferenceType
  )]

  data.table::fwrite(out, file = output_file, sep = "\t", quote = FALSE)
  out
}


