% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadAnnotationFiles.R
\name{loadAnnotationFiles}
\alias{loadAnnotationFiles}
\title{Load annotation data files from package reference directory}
\usage{
loadAnnotationFiles(
  file1 = "FLcircAS.final.RData",
  file2 = "IsoCirc.final.RData",
  file3 = "Homo_sapiens.GRCh38.94.chr.gtf_exontable.rds",
  package = "AQUARIUMHB",
  envir = .GlobalEnv,
  verbose = FALSE
)
}
\arguments{
\item{file1}{Name of the first data file (without path). Defaults to "FLcircAS.final.RData".}

\item{file2}{Name of the second data file (without path). Defaults to "IsoCirc.final.RData".}

\item{file3}{Name of the third data file (without path). Defaults to "Homo_sapiens.GRCh38.94.chr.gtf_exontable.rds".}

\item{package}{Name of the package. Defaults to "AQUARIUMHB".}

\item{envir}{Environment to load the data into. Defaults to the global environment (.GlobalEnv).}

\item{verbose}{If TRUE, prints detailed information about the loading process.}
}
\value{
Invisible NULL. The function loads data into the specified environment rather than returning it.
}
\description{
Loads three specified data files located in the package's inst/reference directory,
including two .RData files and one .rds file.
}
\examples{
# Load default files from AQUARIUMHB package
loadAnnotationFiles()

# Load custom files with verbose output
loadAnnotationFiles("custom1.RData", "custom2.RData", "custom3.rds", verbose = TRUE)
}
