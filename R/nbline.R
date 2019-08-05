#' @title Number & Brightness (Single Image)
#' 
#' @description Performs the Number and Brightness Analysis (N&B) of a single image.
#' @aliases nbline
#' @usage nbline(img, plot = TRUE)
#' @param img The image to analyze.
#' @param plot by default True. This parameter allow to print the visual comparision of the Brightness and Number.
#' @details The N&B method provides molecular concentration and brightness for each pixel of a single image.
#' The average particle number and brightness are calculated directly from
#' the mean value <k> and variance (sigma^2) of the fluorescence intensity data (image) for a given pixel as follows: 
#' N = (<k>^2)/(sigma^2)
#' and
#' B = (sigma^2)/<k>
#' Leaving the plot parameter as True, the function will plot the two vectos in the same plot for a visual comparison, in addition to returning the two vectors.
#' 
#' @export
#' @importFrom stats var
#' @importFrom graphics axis mtext par
#' @return Brightness Number   A list containing two vectors, the Brightness and the Number of the image.
#' @author Raul Pinto Camara.
#' 
#' @seealso \code{\link{var}, \link{mean}}
#' @examples
#' 
#' dmv2 <- data.matrix(v2DataSet)
#' nb <- nbline(dmv2, plot = FALSE)
#' oldpar <- par(no.readonly = TRUE)
#' par(mar = c(5, 5, 3, 5))
#' plot(nb$Brightness, type = 'l', col = "blue", ylab = "Brightness (Blue)", xlab = "Pixel")
#' par(new=TRUE)
#' plot(nb$Number, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylab ="", xlab = "")
#' axis(4)
#' mtext("Number (Red)", side = 4, line = 2)
#' par(oldpar)


nbline <- function(img, plot = TRUE){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  di <- dim(img)
  if(length(di) != 2){
    stop("'img' must be two dimensional")
  }
  mn <- apply(img, MARGIN = 1, mean, na.rm = T)
  vr <- apply(img, MARGIN = 1, var, na.rm = T)
  pB <- vr/mn
  pV <- (mn^2)/vr
  if(plot){
    par(mar = c(5, 5, 3, 5))
    plot(pB, type = 'l', col = "blue", ylab = "Brightness (Blue)", xlab = "Pixel")
    par(new=T)
    plot(pV, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylab ="", xlab = "")
    axis(4)
    mtext("Number (Red)", side = 4, line = 2)
  }
  return(list("Brightness" = pB, "Number" = pV))
}
