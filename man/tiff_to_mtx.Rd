% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tiff_to_mtx.R
\name{tiff_to_mtx}
\alias{tiff_to_mtx}
\title{Transformation of multiple-image TIFF files or arrays into a matrix}
\usage{
tiff_to_mtx(data, columns)
}
\arguments{
\item{data}{A character string indicating the name of a TIFF file or a 3D array}

\item{columns}{The number of columns of the resulting matrix}
}
\value{
A matrix
}
\description{
Transforms multiple-image TIFF files or 3D arrays into 2D matrices with a user-specified number of columns
}
\details{
Creates a matrix with a user-specified number of columns and a number of rows equal to the total amount of points in 'data' divided by 'columns'.
}
\examples{
\donttest{
### Please navigate to
### (https://github.com/FCSlib/FCSlib/tree/master/Sample\%20Data)
### to find this sample data

x <- readFileTiff("Example_file_name.tif")
class(x); dim(x)

x.m <- tiff_to_mtx(data = x, columns = 64)
class(x.m); dim(x.m)
}
}
\seealso{
\code{\link{binMatrix}}
}
\author{
Alejandro Linares
}
