#' Load long isoform data files from package reference directory
#'
#' Loads two specified .RData files located in the package's inst/reference directory.
#'
#' @param file1 Name of the first .RData file (without path). Defaults to "FLcircAS.final.RData".
#' @param file2 Name of the second .RData file (without path). Defaults to "IsoCirc.final.RData".
#' @param package Name of the package. Defaults to "AQUARIUMHB".
#' @param envir Environment to load the data into. Defaults to the global environment (.GlobalEnv).
#' @param verbose If TRUE, prints detailed information about the loading process.
#' @return Invisible NULL. The function loads data into the specified environment rather than returning it.
#' @export
#' @examples
#' # Load default files from AQUARIUMHB package
#' loadLongIsoformFiles()
#' 
#' # Load custom files with verbose output
#' loadLongIsoformFiles("custom1.RData", "custom2.RData", verbose = TRUE)
loadLongIsoformFiles <- function(file1 = "FLcircAS.final.RData", 
                                 file2 = "IsoCirc.final.RData",
                                 package = "AQUARIUMHB",
                                 envir = .GlobalEnv,
                                 verbose = FALSE) {
  # Construct full paths using system.file()
  path1 <- system.file("reference", file1, package = package, mustWork = FALSE)
  path2 <- system.file("reference", file2, package = package, mustWork = FALSE)
  
  # Validate package installation
  if (path1 == "" || path2 == "") {
    stop(paste("Package", package, "not found or files missing. Check package installation."))
  }
  
  # Validate file existence
  if (!file.exists(path1)) {
    stop(paste("File", file1, "not found in package reference directory."))
  }
  
  if (!file.exists(path2)) {
    stop(paste("File", file2, "not found in package reference directory."))
  }
  
  # Verbose output
  if (verbose) {
    message("=== Loading data files ===")
    message("Package:", package)
    message("File 1:", path1)
    message("File 2:", path2)
    message("File sizes:")
    message("- ", file1, ": ", format(file.size(path1), units = "auto"))
    message("- ", file2, ": ", format(file.size(path2), units = "auto"))
  }
  
  # Function to load and return object names
  load_and_report <- function(path, file_name) {
    if (verbose) message("\nLoading ", file_name, "...")
    
    # Capture objects loaded from the file
    obj_names <- load(path, envir = envir)
    
    if (verbose) {
      if (length(obj_names) == 0) {
        message("No objects found in ", file_name)
      } else {
        message("Loaded objects (", length(obj_names), "):")
        for (obj in obj_names) {
          obj_class <- class(get(obj, envir = envir))
          obj_size <- format(object.size(get(obj, envir = envir)), units = "auto")
          message("- ", obj, " (class: ", paste(obj_class, collapse = ", "), ", size: ", obj_size, ")")
        }
      }
    }
    
    return(obj_names)
  }
  
  # Load files and capture loaded object names
  loaded1 <- load_and_report(path1, file1)
  loaded2 <- load_and_report(path2, file2)
  
  # Final summary
  if (verbose) {
    total_objects <- length(loaded1) + length(loaded2)
    message("\n=== Summary ===")
    message("Total objects loaded: ", total_objects)
    message("Environment: ", environmentName(envir))
  }
  
  invisible(NULL)
}