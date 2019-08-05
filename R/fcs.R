#' @title Fluorescence Correlation Spectroscopy
#' 
#' @description Performs either the auto-correlation, or cross-correlation between vectors x and y, returning a correlation function.
#' @aliases fcs
#' @usage fcs(x , y = NULL, nPoints = 25000, pcf = FALSE)
#' @param x Numeric vector of length N.
#' @param y Numeric vector of length N.
#' @param nPoints The size of the sub-vectors in which the input vectors will be divided. This number must be less than N/2.
#' @param pcf A boolean parameter to determine if an alternate version of the correlation function is used for the calculation of de pCF and pComb functions.
#' @details Fluorescence correlation spectroscopy (FCS) is a technique with high spatial and temporal resolution used to analyze the kinetics of particles diffusing at low concentrations. The detected fluorescence intensity as a function of time is: F(t).
#' 
#' The correlation function is computed as the normalized autocorrelation function, G(tau) = <deltaF(t)deltaF(t+tau)>/(<F(t)>*<F(t)>), to the collected data set, where t refers to a time point of flourescence acquisition, and tau refers to the temporal delay between acquisitions and <...> indicates average. 
#' 
#' The correlation between deltaF(t) = F(t) - <F(t)> and deltaF(t+tau) = F(t+tau) - <F(t)> is calculated for a range of delay times.
#' For temporal acquisitions as FCS point, x takes the value of F(t) and y = NULL.
#' For cross-correlation experiments between two fluorescent signals x = F1(t) and y = F2(t), as channels, the correlation function is: G(tau) = <deltaF1(t) deltaF2(t+tau)> / (<F1(t)> <F2(t)>).
#' 
#' The function separate the original vector in sub-vectors of same length (n-points), then calculate an autocorrelation function form each sub-vector. The final result will be an average of all the autocorrelation functions.
#' @note The argument nPoints must be smaller than the total number of temporal observations N, it is recommended to set nPoints = 2^n, with n = {2, ..., infinity}.
#' @export
#' @return A numeric vector G containing either the autocorrelation for the input vector x, or the cross-correlation between x and y vectors, with a length of nPoints.
#' @references R.A. Migueles-Ramirez, A.G. Velasco-Felix, R. Pinto-Cámara, C.D. Wood, A. Guerrero. Fluorescence fluctuation spectroscopy in living cells.
#' Microscopy and imaging science: practical approaches to applied research and education, 138-151,2017.
#' @author Raul Pinto Camara, Adan O. Guerrero
#' 
#' @seealso \code{\link{gcf}}
#' @examples 
#' ## Choose between these two Data Sets: Cy5_1nM, Cy5_10nM.
#' ## These data sets are simulations of the molecular diffusion Cy5 in two distinct concentrations.
#' 
#' ##Here is the execution of the FCS method for the Cy5_1nM Data set
#' g <- fcs(x = Cy5_1nM$f)
#' len <- 1:length(g)
#' oldpar <- par(no.readonly = TRUE)
#' par(mfrow=c(1,1))
#' plot(y = g, x = Cy5_1nM$t[len], log = 'x', type = 'l',
#' xlab = expression(tau(mu~s)), ylab = expression(G(tau)), main = "Cy5 in 1nM")
#' par(oldpar)


fcs <- function(x , y = NULL, nPoints = 25000, pcf = FALSE){
  if(is.null(y)){
    y = x
  }
  nPointsint = nPoints*2
  if(!(is.vector(x)&&is.vector(y))){
    stop("x and y must be vectors")
  }
  if(!(length(x)==length(y))){
    stop("x and y must be vectors of the same length")
  }
  if(nPointsint > length(x)){
    stop("'nPoints' must be less than the length of x/2")
  }
  nVals <- length(x)
  nSect <- trunc(nVals / nPointsint)
  rest <- nVals%%nPointsint
  if(rest != 0){
    x <- x[-((nVals-rest+1):nVals)]; y <- y[-((nVals-rest+1):nVals)]
  }
  G <- array(data = 0, dim = c(nPointsint, nSect))
  for(i in 1:nSect){
    leftLim <- nPointsint * (i - 1) + 1
    rightLim <- nPointsint * i
    xmean <- mean(x[leftLim:rightLim], na.rm = T)
    ymean <- mean(y[leftLim:rightLim], na.rm = T)
    g <- array(data = 0, dim = nPointsint)
    if(xmean != 0 && ymean != 0){
      dx <- x[leftLim:rightLim]
      dy <- y[leftLim:rightLim]
      if(!pcf){
        dx <- dx - xmean
        dy <- dy - ymean
        g <- gcf(x = dx, y = dy, xmean = xmean, ymean = ymean)
      } else{
        g <- gcf(x = dx, y = dy, xmean = xmean, ymean = ymean)
        g <- g - 1
      }
    }
    G[, i] <- g
  }
  G <- apply(G, MARGIN = 1, mean, na.rm = T); G <- G[1:(nPoints+1)]
  return(G[-1])
}
