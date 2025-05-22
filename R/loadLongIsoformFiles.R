#' Load long isoform data files from package reference directory
#'
#' Loads two specified .RData files located in the package's inst/reference directory.
#'
#' @param file1 Name of the first .RData file (without path). Defaults to "FLcircAS.final.RData".
#' @param file2 Name of the second .RData file (without path). Defaults to "IsoCirc.final.RData".
#' @param package Name of the package. Defaults to "AQUARIUMHB".
#' @return Invisible NULL. The function loads data into the environment rather than returning it.
#' @export
#' @examples
#' # Load default files from AQUARIUMHB package
#' loadLongIsoformFiles()
#' 
#' # Load custom files from another package
#' loadLongIsoformFiles("custom1.RData", "custom2.RData", package = "MyPackage")
loadLongIsoformFiles <- function(file1 = "FLcircAS.final.RData", 
                                 file2 = "IsoCirc.final.RData",
                                 package = "AQUARIUMHB") {
  # Construct full paths using system.file()
  path1 <- system.file("reference", file1, package = package, mustWork = FALSE)
  path2 <- system.file("reference", file2, package = package, mustWork = FALSE)
  
  # Validate file existence
  if (!file.exists(path1)) {
    stop(paste("File", file1, "not found in package reference directory."))
  }
  
  if (!file.exists(path2)) {
    stop(paste("File", file2, "not found in package reference directory."))
  }
  
  # Load files into environment
  load(path1)
  load(path2)
  
  invisible(NULL)
}