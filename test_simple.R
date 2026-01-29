# 简单测试 Merge_circRNA.gtf 函数

# 加载必要的库
library(data.table)
library(rtracklayer)
library(dplyr)
library(purrr)

# 加载函数
source("R/Merge_circRNA.gtf.R")

# 创建测试数据路径
sample_dirs <- list.dirs("inst/extdata", recursive = FALSE, full.names = TRUE)

# 只使用前3个样本进行测试
sample_dirs <- sample_dirs[1:3]

# 创建 SamplePath 数据框
SamplePath <- data.table(
  SampleID = basename(sample_dirs),
  FullPath = sample_dirs
)

cat("Testing with", nrow(SamplePath), "samples...\n")

# 测试函数
start_time <- Sys.time()
tryCatch({
  result <- Merge_circRNA.gtf(SamplePath, output_file = "test_output.txt", parallel = FALSE)
  cat("Success! Processed", nrow(result), "records.\n")
  cat("Time taken:", difftime(Sys.time(), start_time, units = "secs"), "seconds\n")
  
  # 查看结果前几行
  cat("\nFirst few rows of result:\n")
  print(head(result))
  
}, error = function(e) {
  cat("Error:", e$message, "\n")
  traceback()
})

# 清理测试文件
if (file.exists("test_output.txt")) {
  file.remove("test_output.txt")
  cat("\nTest file cleaned up.\n")
}