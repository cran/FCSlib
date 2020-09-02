#' @title Exponential detrending of a time series
#' 
#' @description Performs the exponential detrending algorithm over a vector
#' @aliases expDetrend
#' @usage expDetrend(f, acqTime, nIntervals, pois = FALSE, plot = TRUE)
#' @param f A vector
#' @param acqTime Point acquisition rate (in seconds)
#' @param nIntervals Number of intervals into which the vector will be grouped before an exponential model fit 
#' @param pois Boolean, set to TRUE to add random, uncorrelated numbers sampled from a Poisson distribution
#' @param plot Boolean, set to TRUE (default) to plot the result
#' @details First, the binTimeSeries() function is used to obtain a binned version of 'f' of 'nIntervals' points.
#' An exponential model of the form (A0*e^(k*t) is then adjusted to the binned vector. The full time series is then
#' evaluated using the obtained model and the residuals are calculated. When 'pois = FALSE', A0 is directly added
#' to the full time series residuals for trend correction, in a quantity that will make the average counts remain
#' constant across the whole time series. On the other hand, when 'pois = TRUE', then A0 will be used as the lambda
#' parameter for a Poisson distribution from which integer numbers will be randomly sampled and added to the whole
#' series for trend correction. This procedure asures that no photon counts with decimals will be obtained, though,
#' it can add some noise to the data. 
#' 
#' @export
#' @importFrom graphics abline lines points
#' @importFrom stats lm rpois
#' @return A vector
#' @author Alejandro Linares, Adan Guerrero, Haydee Hernández
#' 
#' @seealso \code{\link{polyDetrend} \link{boxcarDetrend} \link{binTimeSeries}}
#' 
#' @examples
#' \donttest{
#' ### Please navigate to
#' ### (https://github.com/FCSlib/FCSlib/tree/master/Sample%20Data)
#' ### to find this sample data
#' 
#' x <- read.table("PB030.dat", header = F)
#' 
#' x.d <- expDetrend(x[,2], acqTime = 4e-6, nIntervals = 100, pois = F)
#' 
#' # Poisson sampling (this might take a little bit longer)
#' x.d.p <- expDetrend(x[,2], acqTime = 4e-6, nIntervals = 100, pois = T)
#' }

expDetrend <- function(f, acqTime, nIntervals, pois = FALSE, plot = TRUE){
  a <- binTimeSeries(f, acqTime, nIntervals, plot = F)
  em <- lm(log(a$Counts)~a$Time)
  A0 <- exp(em$coefficients[1])
  k <- em$coefficients[2]
  model <- A0*exp(k*a$Time)
  residuals <- a$Counts - model
  t <- (1:length(f))*acqTime
  full.model <- A0*exp(k*t)
  if (pois == F){
    full.residuals <- f - full.model + A0
  } else if (pois == T){
    full.residuals <- NULL
    for (i in 1:length(full.model)){
      full.residuals[i] <- f[i] + rpois(1, (A0 - full.model[i]))
    }
  }
  
  b <- binTimeSeries(full.residuals, acqTime, nIntervals, plot = F)
  if (plot){
    plot(a$Counts~a$Time, xlab = "Time (s)", ylab = "Mean photon counts", type = "l", main = "Binned single-point data")
    lines(model~a$Time, lwd = 2, col = "red")
    points(A0 + residuals~a$Time, col = "green", pch = 19)
    lines(b$Counts~b$Time, col = "blue")
    abline(h = model[1], lwd = 2)
  }

  return(full.residuals)
}
