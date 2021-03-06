% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readFileSPC.R
\name{readFileSPC}
\alias{readFileSPC}
\title{Read File SPC-140/150/130/830}
\usage{
readFileSPC(filename, nData = 1E8)
}
\arguments{
\item{filename}{the name of the file to read from.}

\item{nData}{parameter that defines the length of data to read.}
}
\value{
A list containing the Macrotime and the Microtime vectors.
}
\description{
Reads a SPC file and returns the Macrotime and Microtime.
}
\details{
Read a SPC file, with SPC-140/150/130/830 version, using the readBin function and extract the data contained in the file.
}
\note{
The nData parameter is used to overestimate the amount of data that the file can contain.
}
\examples{
\donttest{
spcData <- readFileSPC(FileName)
}
}
\references{
Becker, W., 2019. The Bh TCSPC Handbook. 8th ed. Berlin, Germany: Becker & Hickl GmbH, pp. 855-856.
}
\seealso{
\code{\link{asynACTCSPC}}
}
\author{
Raúl Pinto Cámara, José Damián Martínez Reyes.
}
