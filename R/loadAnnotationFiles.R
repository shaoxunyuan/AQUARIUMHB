#' Load annotation data files from package reference directory
#'
#' Loads three specified data files located in the package's inst/reference directory,
#' including two .RData files and one .rds file.
#'
#' @param file1 Name of the first data file (without path). Defaults to "FLcircAS.final.RData".
#' @param file2 Name of the second data file (without path). Defaults to "IsoCirc.final.RData".
#' @param file3 Name of the third data file (without path). Defaults to "Homo_sapiens.GRCh38.94.chr.gtf_exontable.rds".
#' @param package Name of the package. Defaults to "AQUARIUMHB".
#' @param envir Environment to load the data into. Defaults to the global environment (.GlobalEnv).
#' @param verbose If TRUE, prints detailed information about the loading process.
#' @return Invisible NULL. The function loads data into the specified environment rather than returning it.
#' @export
#' @examples
#' # Load default files from AQUARIUMHB package
#' loadAnnotationFiles()
#' 
#' # Load custom files with verbose output
#' loadAnnotationFiles("custom1.RData", "custom2.RData", "custom3.rds", verbose = TRUE)
loadAnnotationFiles <- function(file1 = "FLcircAS.final.RData", 
                                file2 = "IsoCirc.final.RData",
                                file3 = "Homo_sapiens.GRCh38.94.chr.gtf_exontable.rds",
                                package = "AQUARIUMHB",
                                envir = .GlobalEnv,
                                verbose = FALSE) {
  # Construct full paths using system.file()
  path1 <- system.file("reference", file1, package = package, mustWork = FALSE)
  path2 <- system.file("reference", file2, package = package, mustWork = FALSE)
  path3 <- system.file("reference", file3, package = package, mustWork = FALSE)

  # Validate package installation
  if (path1 == "" || path2 == "" || path3 == "") {
    stop(paste("Package", package, "not found or files missing. Check package installation."))
  }

  # Validate file existence
  validate_file <- function(path, fname) {
    if (!file.exists(path)) {
      stop(paste("File", fname, "not found in package reference directory."))
    }
  }
  validate_file(path1, file1)
  validate_file(path2, file2)
  validate_file(path3, file3)

  # Verbose output
  if (verbose) {
    message("=== Loading data files ===")
    message("Package:", package)
    message("File 1:", path1)
    message("File 2:", path2)
    message("File 3:", path3)
    message("File sizes:")
    message("- ", file1, ": ", format(file.size(path1), units = "auto"))
    message("- ", file2, ": ", format(file.size(path2), units = "auto"))
    message("- ", file3, ": ", format(file.size(path3), units = "auto"))
  }

  # Generic load function with type detection
  load_data <- function(path, fname, verbose) {
    if (verbose) message("\nLoading ", fname, "...")
    
    # Detect file type and load accordingly
    if (grepl("\\.RData$", fname, ignore.case = TRUE)) {
      obj_names <- load(path, envir = envir)
    } else if (grepl("\\.rds$", fname, ignore.case = TRUE)) {
      obj <- readRDS(path)
      obj_names <- basename(sub("\\.rds$", "", fname))
      assign(obj_names, obj, envir = envir)
    } else {
      stop("Unsupported file type. Only .RData and .rds files are allowed.")
    }
    
    if (verbose) {
      if (length(obj_names) == 0) {
        message("No objects found in ", fname)
      } else {
        message("Loaded objects (", length(obj_names), "):")
        for (obj in obj_names) {
          if (exists(obj, envir = envir)) {
            obj_class <- class(get(obj, envir = envir))
            obj_size <- format(object.size(get(obj, envir = envir)), units = "auto")
            message("- ", obj, " (class: ", paste(obj_class, collapse = ", "), ", size: ", obj_size, ")")
          }
        }
      }
    }
    
    return(obj_names)
  }

  # Load all three files
  loaded1 <- load_data(path1, file1, verbose)
  loaded2 <- load_data(path2, file2, verbose)
  loaded3 <- load_data(path3, file3, verbose)

  # Final summary
  if (verbose) {
    total_objects <- length(loaded1) + length(loaded2) + length(loaded3)
    message("\n=== Summary ===")
    message("Total objects loaded: ", total_objects)
    message("Environment: ", environmentName(envir))
  }

  invisible(NULL)
}