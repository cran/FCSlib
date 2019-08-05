#' @title Fitting FCS Data
#' 
#' @description Estimates the parameters based on a given equation, on the data generated with the fcs() function.
#' @aliases fitFCS
#' @usage fitFCS(data = parent.frame(), start, low = -Inf, up = Inf, type = "D3D", trace = TRUE)
#' @param data data frame in which to evaluate the variables in formula and weights.
#' @param start a named list or named numeric vector of starting estimates.
#' @param low,up a named list or named numeric vector of lower and upper bounds, replicated to be as long as start. If unspecified, all parameters are assumed to be -Inf and Inf.
#' @param type specification for the equation to model, is a character string. The default value is "D3D" equation for three-dimensional free diffusion.
#' Another possibles values are: "D2D" for two-dimensional free diffusion,  "D2DT" for two-dimensional free diffusion with triplet exited state, and "D3DT" for three-dimensional free diffusion with triplet exited state.
#' @param trace logical value that indicates whether the progress of the non-linear regression (nls) should be printed.
#' @details The autocorrelation function contains information about the molecular diffusion coefficient and the number of molecules occupying the observation volume. To interpret such data in molecular terms we need a model to describe the fluctuations. The model of free diffusion in three dimensions:
#' G(tau)=G(0)((1+(tau/tau_D))^(-1))(1+((s/u)^2)(tau/tau_D))^(-1/2), with tau_D=(s^2)/4D
#' where tau_D is the diffusion time, s is the radius and u = ks.
#' So, if "D3D" is selected, the variables that should be considered are G (resulting value of the fcs function), tau the vector of time in the correlation function, s the PSF (Point Spread Function) radius, k the scalar factor between the radius half-length of the PSF, D is the diffusion coefficent and the G(0) is inversely proportional to the molecular concentration.
#' 
#' @export
#' @importFrom stats nls formula
#' @return A nls object (from nls).
#' @author Raul Pinto Camara
#' 
#' @seealso \code{\link{nls}}, \code{\link{fcs}}
#' @examples 
#' g <- fcs(x = Cy5_1nM$f)
#' len <- 1:length(g)
#' 
#' df <- data.frame(g = g, tau = Cy5_1nM$t[len], s = 0.25, k = 3)
#' start <- list(D = 100, G0 = 0.1)
#' low <- list(D = 1E-1, G0 = 1E-2)
#' up  <- list(D = 500, G0 = 100)
#' 
#' modelFCS <- fitFCS(data = df, start = start, low = low, up = up, type = "D3D")
#' fit <- predict(modelFCS, list(tau = Cy5_1nM$t[len]))
#' 
#' plot(y = g, x = Cy5_1nM$t[len], log = 'x', type = 'l', xlab = expression(tau(mu~s)),
#' ylab = expression(G(tau)), main = "Symulation of Cy5 [1nM]")
#' lines(fit~Cy5_1nM$t[len], col = "blue")


fitFCS <- function(data = parent.frame(), start, low = -Inf, up = Inf, type = "D3D", trace = TRUE){
  if(!is.list(start)){
    stop("'start' must be a list")
  }
  if (!is.list(data) && !is.environment(data)){
    stop("'data' must be a list or an environment")
  }
  GTrip <- "(1+((B*exp(-(tau)/taub))/(1-B)))"
  D2D <- "((G0) * ((1 + ((4 * D * tau)/(s ^ 2))) ^ (-1))"
  D2DT <- paste(GTrip, "*", D2D)
  D3D <- "((G0) * ((1 + ((4 * D * tau)/(s ^ 2))) ^ (-1)) * ((1 + ((4 * D * tau)/((k^2) * (s^2)))) ^ (-1 / 2)))"
  D3DT <- paste(GTrip, "*", D3D)
  list_type <- list("D2D", "D2DT", "D3D", "D3DT"); list_model <- list(D2D, D2DT, D3D, D3DT)
  if(!type %in% list_type){
    stop(paste("The model", type, "doesn't exist"))
  } else{
    typePosList <- which(list_type == type); form <- paste("g ~ ", list_model[[typePosList]], sep = "")
  }
  express <- formula(form)
  modelFCS <- nls(formula = express, data = data, start = start, trace = trace, algorithm = "port", lower = low, upper = up)
  return(modelFCS)
}