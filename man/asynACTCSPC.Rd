% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asynACTCSPC.R
\name{asynACTCSPC}
\alias{asynACTCSPC}
\title{Asynchronous Autocorrelation of Time-Correlated Single-Photon Counting}
\usage{
asynACTCSPC(macro, n = 5, B = 10)
}
\arguments{
\item{macro}{A numeric vector containig a Macrotime Data.}

\item{n}{numeric parameter that represents the number of layers of the cascade.}

\item{B}{numeric parameter that represents the number of values in every layer of the cascade.}
}
\value{
A numeric vector G containing either the autocorrelation for the input vector macro, with a length of n*B.
}
\description{
Calculates the auto-correlation of the Macrotime data, returning a correlation function.
}
\details{
This function creates list of tau's with a length of n*B, this list is used to perform the correlation of the data.
}
\examples{
\donttest{
spcData <- readFileSPC("atto532_atto655_m1.spc")
asynCorrData <- asynACTCSPC(macro = spcData$MacroTime)
}
}
\references{
wahl, M., Gregor, I., Patting, M. & Enderlein, J. Fast calculation of fluorescence correlation data with
asynchronous time-correlated single-photon counting. Opt. Express 11, 3583–3591 (2003).
}
\seealso{
\code{\link{readFileSPC}}
}
\author{
Raúl Pinto Cámara, José Damián Martínez Reyes.
}
