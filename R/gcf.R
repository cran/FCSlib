#' @title General Correlation Function

#' @aliases gcf
#' @description Performs either the auto-correlation, or cross-correlation between vectors x and y, returning a correlation function.
#' @usage gcf(x, y, xmean = 1, ymean = 1, c = 0)
#' @param x A numerical signal with dimensions M x N x Z.
#' @param y A numerical signal with dimensions M x N x Z.
#' @param xmean The mean value of the signal x.
#' @param ymean The mean value of the signal y.
#' @param c A numeric variable to restrict the correlation to positives values.
#' @details The number of emission events per unit time is determined and used to generate autocorrelation and cross-correlation curves from the intensity traces F(t) and the fluctuations deltaF(t) = F(t)-<F(t)>.
#' The auto-correlation function of the collected data set, is computed as the normalized auto-correlation function, when y=x. The general auto-correlation function is defined as:
#' G(tau) = (deltaF(t) deltaF(t+tau) )/(<F(t)> <F(t)>),
#' where t refers to a time point of fluorescence acquisition, and tau refers to the temporal delay between acquisitions. <...> is the temporal average of F(t); and deltaF(t) = F(t)-<F(t)>, deltaF(t+tau) = F(t+tau)-<F(t)>.
#' 
#' For temporal acquisitions such as point FCS, x and y are F(t).
#' The cross-correlation function between two channels of fluorescent signals, x = F1(t) and y = F2(t), the cross-correlation function is defined as:
#' G(tau) = (deltaF1(t) deltaF2(t+tau) )/(<F1(t)><F2(t)>),
#' where xmean = <F1(t)> and ymean = <F2(t)> are the mean values of the fluorescent signals.
#' 
#' @export
#' @importFrom stats convolve
#' @return G   A numerical signal with dimension N' x M' x Z'
#' @references Siegel, A. P., Hays, N. M., & Day, R. N. (2013). Unraveling transcription factor interactions with heterochromatin protein 1 using fluorescence lifetime imaging microscopy and fluorescence correlation spectroscopy. Journal of biomedical optics, 18(2), 025002.
#' @author Raul Pinto Camara.
#' 
#' @seealso \code{\link{fcs}}, \code{\link{convolve}}
#' @examples
#' oldpar <- par(no.readonly = TRUE)
#' g <- gcf(x = Cy5_1nM$f, y = Cy5_1nM$f, xmean = mean(Cy5_1nM$f), ymean = mean(Cy5_1nM$f))
#' len <- 1:length(g)
#' par(mfrow=c(1,1))
#' plot(y = g, x = Cy5_1nM$t[len], log = 'x', type = 'l',
#' xlab = expression(tau(mu~s)), ylab = expression(G(tau)),
#' main = "Cy5 in 1nM")
#' par(oldpar)

gcf <- function(x, y, xmean = 1, ymean = 1, c = 0){
  N <- length(x)
  G <- convolve(x, y) / (N * xmean * ymean)
  G <- G + c
  return(G)
}
