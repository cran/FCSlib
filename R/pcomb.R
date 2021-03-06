#' @title Pair Correlation of Molecular Brightness
#' 
#' @description Performs the pair correlation of molecular brightness (pCOMB) analysis.
#' @aliases pcomb
#' @usage pcomb(img, nPoints = 25000, one.col = FALSE, dr = 1, w = 100)
#' @param img The image to analyze.
#' @param nPoints The size of the sub-vectors in which the input vectors will be divided. This number must be less than N/2.
#' @param one.col By default FALSE. If TRUE the correlation will be performed in the fixed colum mode, else the distance mode.
#' @param dr Is the distance between the two columns that will be correlated. For a value of deltar = 3, the columns are correlated as follows: (1,4), (2,5), ..., (n-3, n), with n as the last column.
#' @param w Range value that is used to calculate the brightness in the image.
#' @details With the Pair Correlation of Molecular Brightness (pCOMB) method, one can distinguish between different homo-oligomeric species of the same molecule coexisting in the same microenvironment, while separately and specifically tracking each species' moblity across the cellular compartments. This technique amplifies the signal from the brightest species present and filters the dynamics of the extracted oligomeric population based on arrival time between two locations. This method is suitable for mapping the impact oligomerization on transcription factor dynamics.
#' The resulting intensity fluctuations, pCF, are transformed into brightness fluctuations using B = (sigma^2)/mean, and the pair correlation analysis is then performed on the brightness fluctuations along the line scan , at a distance (delta(r)).
#' 
#' If the pcf is set as FALSE the pComb data will not be generated and will be NULL. In order to generate that data the pcf function must be used on the BCarpet data.
#' @export
#' @importFrom stats var
#' @return A list containing the Brightness Carpet and the Pair Correlation of that carpet
#' @author Raúl Pinto Cámara.
#' 
#' @seealso \code{\link{fcs}, \link{pcf}}
#' 
#' @examples
#' \donttest{
#' ### Load the FCSlib package
#' 
#' library(FCSlib)
#' 
#' # As an example, we will use a data set that corresponds to a population of Venus dimers
#' # diffusing in HEK-293 cells. Use the readFileTiff() function to extract the information
#' # from the '.tiff' files.
#' 
#' dmv2 <- data.matrix(V2)
#' pC <- pcomb(dmv2[1:32,1:2001], nPoints = 1000, type = 'd', dr = 10, w = 2, pcf = FALSE)
#' dmv2 <- data.matrix(v2DataSet)
#' pC <- pcomb(dmv2, nPoints = 5000, type = 'd', dr = 10, w = 100)
#' di <- dim(pC$pComb)
#' tau <- (1:(di[2]))
#' 
#' # Plot the result
#' library("fields")
#' image.plot( x = 1:di[1], y = log10(tau), z = pC$pComb, main = "pComb",
#' xlab = "Pixel", ylab = "Logarithmic tau",
#' cex.lab = 1.2, cex.main = 1.2, cex.axis = 1)
#' }

pcomb <- function(img, nPoints = 25000, one.col = FALSE, dr = 1, w = 100){
  di <- dim(img)
  if(length(di) != 2){
    stop("'img' must be two dimensional")
  }
  if((nPoints*2) > di[2]){
    stop("'nPoints' must be less than the length of the second dimension of image / 2")
  }
  if(one.col){
    if(dr <= 0 || dr > di[1]){
      stop(paste("'dr' must be dr > 0 and dr <= ", di[1], sep = ""))
    }
  } else{
    if(dr < 0 || dr > di[1]){
      stop(paste("'dr' must be dr >= 0 and dr <= ", di[1], sep = ""))
    }
  }
  if(w < 2){
    stop("w must be at least w = 2")
  }
  pB <- array(data = NA, dim = c(di[1], di[2]-w+1))
  for(i in 1:(di[2]-w+1)){
    pB[,i] <- apply(img[,i:(i+w-1)], MARGIN = 1, var)/apply(img[,i:(i+w-1)], MARGIN = 1, mean)
  }
  pC <- NULL
  pC <- pcf(pB, nPoints, one.col, dr)
  return(pC)
}