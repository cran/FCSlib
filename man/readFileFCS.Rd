% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readFileFcs.R
\name{readFileFCS}
\alias{readFileFCS}
\title{Read File FCS}
\usage{
readFileFCS(filename)
}
\arguments{
\item{filename}{the name of the file to read from.}
}
\value{
dataList   A list containing the data sets within the file.
}
\description{
Reads a FCS file and returns the data sets within the file.
}
\details{
Read a FCS file using the scan function and extract the data contained in the file.
}
\examples{
\donttest{
raw_fcs <- readFileFCS(FileName)
}
}
\author{
Raúl Pinto Cámara.
}
