#' @title Boxcar filter detrending of a time series
#' 
#' @description Performs the boxcar filter detrending algorithm over a vector
#' @aliases boxcarDetrend
#' @usage boxcarDetrend(f, w, acqTime, nIntervals, plot = TRUE)
#' @param f A vector
#' @param w Size of the time window of the moving average
#' @param acqTime Point acquisition rate (in seconds).
#' @param nIntervals Number of intervals into which the moving average vector will be grouped.
#' @param plot Boolean, set to TRUE (default) to plot the result
#' @details First, an amount of zeroes equal to (w-1) is added at the tail of 'f' to compensate for the moving average
#' effect when position (length(f) - w + 1) is reached. The moving average is then calculated and subtracted from the
#' original vector 'f' to obtain the residuals. The moving average vector is then binned and the first value of the
#' resulting vector is used as the lambda parameter for a Poisson distribution from which integer numbers will be
#' randomly sampled and added to the residuals vector for trend correction.
#' 
#' @export
#' @importFrom graphics abline lines points
#' @importFrom stats rpois
#' @return Detrended version of 'f'
#' @author Alejandro Linares, Adan Guerrero, Haydee Hernández
#' 
#' @seealso \code{\link{expDetrend} \link{polyDetrend} \link{binTimeSeries}}
#' 
#' @examples
#' \donttest{
#' ### Please navigate to
#' ### (https://github.com/FCSlib/FCSlib/tree/master/Sample%20Data)
#' ### to find this sample data
#' 
#' x <- read.table("PB030.dat", header = F)
#' 
#' x.d <- boxcarDetrend(x[,2], w = 100, acqTime = 4e-6, nIntervals = 100, plot = TRUE)
#' }

boxcarDetrend <- function(f, w, acqTime, nIntervals, plot = TRUE){
  if (w > (length(f)/2) | w <= 1){
    stop("'w' must be greater than 1 and smaller than half the length of 'f'")
  }
  w <- floor(w)
  
  if (missing(acqTime) | missing(nIntervals)) stop("'acqTime' and 'nIntervals' missing")

  x <- rep(0, length(f))
  f_zp <- c(f, rep(0, w - 1))
  for (i in 1:(length(f)-w+1)){
    x[i] <- mean(f_zp[i:i+w-1])
  }
  full.residuals <- f - x
  b <- binTimeSeries(x, acqTime, nIntervals, plot = F)
  y <- full.residuals + rpois(length(full.residuals), b$Counts[1])

  if (plot){
    a <- binTimeSeries(f, acqTime, nIntervals, plot = F)
    c <- binTimeSeries(y, acqTime, nIntervals, plot = F)
    plot(a$Counts~a$Time, xlab = "Time (s)", ylab = "Mean photon counts", type = "l", 
         ylim = c(min(b$Counts),max(c$Counts)), main = "Binned single-point data")
    lines(b$Counts~b$Time, col = "red")
    points(c$Counts~c$Time, col = "green", pch = 19)
    lines(c$Counts~c$Time, col = "blue")
    abline(h = b$Counts[1])
  }

  return(y)
}