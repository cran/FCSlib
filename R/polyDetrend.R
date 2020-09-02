#' @title Polynomial detrending of a time series
#' 
#' @description Performs the polynomial detrending algorithm over a vector
#' @aliases polyDetrend
#' @usage polyDetrend(f, acqTime, nIntervals, degree = 3, plot = TRUE)
#' @param f A vector
#' @param acqTime Point acquisition rate (in seconds)
#' @param nIntervals Number of intervals into which the vector will be grouped before a polynomial model fit 
#' @param degree The degree of the polynomial function
#' @param plot Boolean, set to TRUE (default) to plot the result
#' @details First, the binTimeSeries() function is used to obtain a binned version of 'f' of 'nIntervals' points.
#' A polynomial model of user-specified degree is then adjusted to the binned vector. The full time series is then
#' evaluated using the obtained model and the residuals are calculated. Finally, the maximum value of the fitted time
#' series is added to the residuals for trend correction.
#' 
#' @export
#' @importFrom graphics abline lines points
#' @importFrom stats lm
#' @return A vector
#' @author Alejandro Linares, Adan Guerrero, Haydee Hernández
#' 
#' @seealso \code{\link{expDetrend} \link{boxcarDetrend} \link{binTimeSeries}}
#' 
#' @examples
#' \donttest{
#' ### Please navigate to
#' ### (https://github.com/FCSlib/FCSlib/tree/master/Sample%20Data)
#' ### to find this sample data
#' 
#' x <- read.table("PB030.dat", header = F)
#' 
#' x.d <- polyDetrend(x[,2], acqTime = 4e-6, nIntervals = 100, degre = 5)
#' }

polyDetrend <- function(f, acqTime, nIntervals, degree = 3, plot = TRUE){
  a <- binTimeSeries(f, acqTime, nIntervals, plot = F)
  model <- lm(a$Counts~poly(a$Time, degree = degree))
  residuals <- model$residuals + model$fitted.values[1]
  t <- (1:length(f))*4e-6
  full.model <- lm(f~poly(t, degree = degree))
  full.residuals <- full.model$residuals + max(model$fitted.values)
  b <- binTimeSeries(full.residuals, acqTime, nIntervals, plot = F)
  
  if (plot){
    plot(a$Counts~a$Time, xlab = "Time (s)", ylab = "Mean photon counts", type = "l", main = "Binned single-point data")
    lines(model$fitted.values~a$Time, lwd = 2, col = "red")
    points(residuals~a$Time, col = "green", pch = 19)
    lines(b$Counts~b$Time, col = "blue")
    abline(h = model$fitted.values[1], lwd = 2)
  }
  
  return(full.residuals)
}
