#' @title Pair Correlation Function
#' 
#' @description Calculates the correlation between the pixel i and pixel i + dr.
#' @aliases pcf
#' @usage pcf(img, nPoints = 1000, type = "d", dr = 1)
#' @param img The image to analyze 
#' @param nPoints The size of the sub-vectors in which the input vectors will be divided. This number must be less than N/2.
#' @param type Possible values are 'l' and 'c'. Set to 'l' by default
#' @param dr Distance between pixel at which the correlation is calculated. For a value of delta_r = 3, the columns are correlated as follows: (1,4), (2,5), ..., (n-3, n), with n being the last column.
#' @details The pair correlation function (pCF) analyzes data of a periodically scanned line, measuring the time it takes a particle to go from one pixel to another, i.e. calculates the spatial cross-correlation function
#' between pixels.
#' G(tau,deltar) = (<F(t,0) F(t + tau, deltar)>/<F(t,0)> <F(t,deltar)>)-1
#' 
#' @export
#' @return An image depicting the correlation between the pixel i and pixel i + dr.
#' @author Raul Pinto Camara.
#' 
#' @seealso \code{\link{fcs}, \link{pcomb}}
#' @examples 
#' \donttest{
#' ### Load the FCSlib package
#' 
#' library(FCSlib)
#' 
#' ### As an example, we will use a data set that corresponds to a population of Venus dimers
#' # diffusing in HEK-293 cells. Use the readFileTiff() function to extract the information
#' # from the '.tiff' files.
#' 
#' dmv2 <- data.matrix(V2)
#' pB <- pcf(dmv2, nPoints = 2500, dr = 10, type = 'd')
#' 
#' ### Plot the result
#' library(fields)
#' di <- dim(pB)
#' tau <- (1:(di[2]))
#' image.plot( x = 1:di[1], y = log10(tau), z = pB, main = "Column Distance 10",
#' xlab = "Pixel", ylab = "Logarithmic tau",
#' cex.lab = 1.2, cex.main = 1.2, cex.axis = 1)
#' 
#' ### Pair Correlation Function (Autocorrelation mode)
#' dmv2 <- data.matrix(v2DataSet)
#' pB <- pcf(dmv2, nPoints = 5000, type = 'a')
#' 
#' ### Plot the result
#' library("fields")
#' di <- dim(pB)
#' tau <- (1:(di[2]))
#' image.plot( x = 1:di[1], y = log10(tau), z = pB, main = "Autocorrelation",
#' xlab = "Pixel", ylab = "Logarithmic tau",
#' cex.lab = 1.2, cex.main = 1.2, cex.axis = 1)
#' 
#' ### Pair Correlation Function (Fixed column mode)
#' dmv2 <- data.matrix(v2DataSet)
#' pB <- pcf(dmv2, nPoints = 5000, dr = 10, type = 'c')
#' 
#' ### Plot the result
#' library("fields")
#' di <- dim(pB)
#' tau <- (1:(di[2]))
#' image.plot( x = 1:di[1], y = log10(tau), z = pB, main = "Fixed Colum 10",
#' xlab = "Pixel", ylab = "Logarithmic tau",
#' cex.lab = 1.2, cex.main = 1.2, cex.axis = 1)
#' }

pcf  <- function(img, nPoints = 1000, type = "d", dr = 1){
  di <- dim(img)
  if(length(di) != 2){
    stop("'img' must be two dimensional")
  }
  if((nPoints*2) > di[2]){
    stop("'nPoints' must be less than the length of the second dimension of image / 2")
  }
  if(tolower(type) != "d" && tolower(type) != "c" && tolower(type) != "a"){
    stop("'type' must be 'd' or 'c' or 'a'")
  }
  if(dr < 0 || dr > di[1]){
    stop(paste("'dr' must be dr >= 0 and dr <= ", di[1], sep = ""))
  }
  if(tolower(type) == "d"){
    pF <- array(data = NA, dim = c(di[1]-dr, nPoints))
    for( i in 1:(di[1]-dr)){
      pF[i,] <- fcs(x = img[i,], y = img[i+dr,], nPoints = nPoints, pcf = TRUE)
    }
  } else if(tolower(type) == "c"){
    pF <- array(data = NA, dim = c(di[1], nPoints))
    for( i in 1:di[1]){
      pF[i,] <- fcs(x = img[dr,], y = img[i,], nPoints = nPoints, pcf = TRUE)
    }
  } else if(tolower(type) == "a"){
    pF <- array(data = NA, dim = c(di[1], nPoints))
    for( i in 1:di[1]){
      t <- fcs(x = img[i,], nPoints = nPoints, pcf = TRUE)
      pF[i,] <- t
    }
  } else{
    stop("'type' must be 'd' or 'c' o 'a'")
  }
  return(pF)
}
