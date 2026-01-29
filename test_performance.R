#!/usr/bin/env Rscript

# 测试 Merge_circRNA.gtf 函数性能

# 加载必要的库
library(data.table)
library(rtracklayer)
library(dplyr)
library(purrr)

# 加载函数
source("R/Merge_circRNA.gtf.R")

# 创建测试数据路径
cat("Creating test sample paths...\n")
sample_dirs <- list.dirs("inst/extdata", recursive = FALSE, full.names = TRUE)

# 只使用前5个样本进行测试（避免测试时间过长）
if (length(sample_dirs) > 5) {
  sample_dirs <- sample_dirs[1:5]
}

# 创建 SamplePath 数据框
SamplePath <- data.table(
  SampleID = basename(sample_dirs),
  FullPath = sample_dirs
)

cat("Testing with", nrow(SamplePath), "samples...\n")

# 测试1：使用并行处理
cat("\n=== Test 1: With parallel processing ===\n")
start_time <- Sys.time()
tryCatch({
  result1 <- Merge_circRNA.gtf(SamplePath, output_file = "test_output_parallel.txt", parallel = TRUE)
  cat("Success! Processed", nrow(result1), "records.\n")
}, error = function(e) {
  cat("Error:", e$message, "\n")
  print(e)
})
end_time <- Sys.time()
parallel_time <- difftime(end_time, start_time, units = "secs")
cat("Time taken:", as.numeric(parallel_time), "seconds\n")

# 测试2：不使用并行处理
cat("\n=== Test 2: Without parallel processing ===\n")
start_time <- Sys.time()
tryCatch({
  result2 <- Merge_circRNA.gtf(SamplePath, output_file = "test_output_sequential.txt", parallel = FALSE)
  cat("Success! Processed", nrow(result2), "records.\n")
}, error = function(e) {
  cat("Error:", e$message, "\n")
  print(e)
})
end_time <- Sys.time()
sequential_time <- difftime(end_time, start_time, units = "secs")
cat("Time taken:", as.numeric(sequential_time), "seconds\n")

# 比较结果
cat("\n=== Performance Comparison ===\n")
cat("Parallel processing:", as.numeric(parallel_time), "seconds\n")
cat("Sequential processing:", as.numeric(sequential_time), "seconds\n")
if (exists("parallel_time") && exists("sequential_time")) {
  speedup <- as.numeric(sequential_time) / as.numeric(parallel_time)
  cat("Speedup:", round(speedup, 2), "x faster with parallel processing\n")
}

# 清理测试文件
if (file.exists("test_output_parallel.txt")) file.remove("test_output_parallel.txt")
if (file.exists("test_output_sequential.txt")) file.remove("test_output_sequential.txt")

cat("\nTest completed!\n")
