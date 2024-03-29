% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadCountsFromH5file.R
\name{loadCountsFromHSDS}
\alias{loadCountsFromHSDS}
\title{Load count matrix from HDF5-files.}
\usage{
loadCountsFromHSDS(
  es,
  url = "https://alserglab.wustl.edu/hsds/?domain=/counts"
)
}
\arguments{
\item{es, }{containing ExpressionSet loaded from GEO. Contains empty expression matrix.}

\item{url, }{containing url of the server and root domain.}
}
\value{
ExpressionSet with loaded count matrix
}
\description{
Load count matrix from HDF5-files.
}
\examples{
ess <- GEOquery::getGEO("GSE85653")
es <- ess[[1]]
url <- 'https://alserglab.wustl.edu/hsds/?domain=/counts'
es <- loadCountsFromHSDS(es, url)

}
