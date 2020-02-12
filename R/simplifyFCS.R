#' @title Simplify FCS
#' 
#' @description Reduces the amount of data in a data set without altering its overall structure
#' @aliases simplifyFCS
#' @usage simplifyFCS(g, tau)
#' @param g A vector containing the FCS data analysis
#' @param tau A vector that represents the time frame between data acquisitions
#' @details Allows to significantly reduce the points of the autocorrelation vector, maintaining its overall structure and allowing to further adjust physical models while obtaining consistent results.
#' @export 
#' @return A vector of the FCS data with reduced length
#' @author Adan O. Guerrero
#' 
#' @seealso \code{\link{gcf}}, \code{\link{var}, \link{mean}}
#' 
#' @examples
#' \donttest{
#' f <- Cy5_100nM$f
#' acqTime <- 2E-6
#' f <- as.vector(f)

#' time <- (1:length(f))*acqTime
#' cy5 <- data.frame(t = time, f)
#' 
#' len <- 1:length(g)
#' tau <-cy5$t[len]
#' G <- data.frame(tau,g)#' g <- fcs(x = cy5$f)
#' 
#' sfcs <- simplifyFCS(G$g, G$tau)
#' plot(df$g~df$tau, log = "x", type = "l",
#'      xlab = expression(tau(s)),
#'      ylab = expression(G(tau)), main = "Cy5")
#' }

simplifyFCS<-function(g, tau){
  df <- data.frame(g, tau)
  t1<-ceiling(log10(min(df$tau)))
  t2<-ceiling(log10(max(df$tau)))
  tR<-10^(t1:t2)
  ts<-NULL; for (i in tR){ ts<-c(ts, (2:10)*i)}
  s<-NULL
  for(i in 1:length(ts)){
    if(i==1) {
      idx<-which(df$tau <= ts[i])
      s<-df[idx,]
    } else{
      idx<-which(df$tau > ts[i-1] & df$tau <= ts[i])
      if(length(idx)) s<-rbind(s,apply(df[idx,],MARGIN = 2, mean) ) 
    }
  }
  return(s)
}