#' Load Sample Paths Automatically
#' 
#' Automatically find SamplePathFile.txt in extdata, read sample names, and generate their full paths.
#' 
#' @return A data frame with sample names and their full paths.
#' @export
#' 
#' @examples
#' ExampleSamplePath <- loadSamplePathFile()
loadSamplePathFile <- function() {
  spf_path <- system.file("extdata", "SamplePathFile.txt", package = "AQUARIUMHB")
  
  if (spf_path == "" || !file.exists(spf_path)) {
    stop("SamplePathFile.txt not found in extdata of AQUARIUMHB package.")
  }
  
  sample_names <- data.table::fread(spf_path, header = TRUE)$SamplePath
  
  ExampleSamplePath <- data.frame(
    SampleName = sample_names,
    FullPath = system.file("extdata", sample_names, package = "AQUARIUMHB"),
    stringsAsFactors = FALSE
  )
  
  return(ExampleSamplePath)
}