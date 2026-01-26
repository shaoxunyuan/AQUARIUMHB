#' Generate circRNA Full GTF File (Optimized)
#'
#' @param SamplePath Data of input, containing sample information, must have columns: SampleName, FullPath.
#' @param ReferenceSet Referenceset data of all possible full-length isoforms.
#'
#' @import data.table
#' @export
circRNA_full.gtf <- function(SamplePath, ReferenceSet) {
  
  # 防止科学计数法输出
  old_scipen <- options(scipen = 999)$scipen
  on.exit(options(scipen = old_scipen))
  
  # 确保必要的包已加载
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  # 预处理参考集：为了极速匹配，建立 key
  ReferenceSet <- as.data.table(ReferenceSet)
  setkey(ReferenceSet, isoformID)

  SamplePath <- as.data.table(SamplePath)

  # Process each sample
  for (i in 1:nrow(SamplePath)) {
    SampleName <- SamplePath$SampleName[i]
    message(paste0("Processing sample: ", SampleName))
    
    dirquant <- paste0(SamplePath$FullPath[i], "quant/")
    if (!dir.exists(dirquant)) {
      message(paste0("Creating directory: ", dirquant))
      dir.create(dirquant, recursive = TRUE)
    }
    
    stout.list.path <- file.path(SamplePath$FullPath[i], "vis/stout.list")
    
    if (!file.exists(stout.list.path)) {
      warning(paste0("stout.list file not found for sample ", SampleName, ". Skipping."))
      next
    }
    
    # 使用 fread 快速读取
    message(paste0("Reading stout.list for sample ", SampleName))
    stout.list <- fread(stout.list.path, header = FALSE)
    
    setnames(stout.list, c("Image_ID", "bsj", "chr", "start", "end", "total_exp",
                           "isoform_number", "isoform_exp", "isoform_length", 
                           "isoform_state", "strand", "gene_id", "isoform_cirexon"))
    
    # 筛选 'Full' 异构体
    type_full <- stout.list[isoform_state == "Full", 
                            .(chr, start, end, strand, bsj, isoform_state, isoform_cirexon)]
    
    if (nrow(type_full) == 0) {
      warning(paste0("No 'Full' isoforms found for sample ", SampleName))
      next
    }
    
    message(paste0("Processing ", nrow(type_full), " 'Full' isoforms..."))

    # --- 核心优化：向量化解析外显子坐标 ---
    
    # 1. 预处理 BSJ 格式 (对应原代码逻辑 chr:start|end)
    type_full[, bsj_formatted := paste0(chr, ":", start, "|", end)]
    
    # 2. 向量化展开外显子 (利用 data.table 的快拆功能)
    # 首先解析出 isoformID 需要的格式，并保留 BSJ
    gtf_expanded <- type_full[, {
      # 拆分外显子字符串 "start-end,start-end"
      exon_pairs <- unlist(strsplit(isoform_cirexon, ","))
      # 进一步拆分 start 和 end
      exon_coords <- tstrsplit(exon_pairs, "-")
      e_starts <- as.numeric(exon_coords[[1]])
      e_ends <- as.numeric(exon_coords[[2]])
      
      # 构造用于匹配 ReferenceSet 的 isoformID 字符串
      iso_id <- paste0("chr", chr, "|", paste(e_starts, collapse = ","), "|", paste(e_ends, collapse = ","), "|", strand)
      
      # 返回展开后的列表，data.table 会自动处理行对应关系
      list(
        start = e_starts,
        end = e_ends,
        isoformID = iso_id,
        bsj_id = bsj_formatted
      )
    }, by = .(chr, strand, bsj_id_raw = bsj)] # 这里的 by 保证了每一组(异构体)内部处理坐标

    # --- 核心优化：向量化映射 ReferenceSource ---
    
    # 利用 data.table 的 join 功能取代 plyr::mapvalues (在数万行数据下 join 快得多)
    # 只取 ReferenceSet 的唯一对应关系，防止 join 产生多余行
    ref_map <- unique(ReferenceSet[, .(isoformID, ReferenceSource)])
    setkey(ref_map, isoformID)
    
    gtf_expanded <- ref_map[gtf_expanded, on = "isoformID"]

    # --- 核心优化：向量化构造 GTF 属性列 ---
    
    # 构造属性字符串
    gtf_expanded[, attr := paste0('bsj "', bsj_id, '"; ',  
                                  'transcript_id "', isoformID, '"; ',  
                                  'isoform_state "Full"; ',  
                                  'ReferenceSource "', ReferenceSource, '";')]

    # 构造最终 GTF 表格
    type_full.gtf <- gtf_expanded[, .(
      chr = chr,
      ciri = "ciri",
      type = "exon",
      start = start,
      end = end,
      attr1 = ".",
      strand = strand,
      attr2 = ".",
      attr = attr
    )]

    # --- 输出处理 ---
    output_file <- paste0(dirquant, "circRNA_full.gtf")
    message(paste0("Writing output to: ", output_file))
    
    # 使用 fwrite 导出，它对数字处理非常快且默认不写科学计数法
    fwrite(type_full.gtf, file = output_file, sep = "\t", quote = FALSE, 
           col.names = FALSE, row.names = FALSE)
    
    message(paste0("Completed processing for sample ", SampleName))
  }
  
  message("All samples processed successfully!")
  return(invisible(NULL))
}