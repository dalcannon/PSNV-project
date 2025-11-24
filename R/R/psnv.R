#' Piecewise Standard Normal Variate (PSNV)
#'
#' @description
#' Calculate piecewise standard normal variate correction (PSNV), i.e. SNV evaluated using a rolling average and rolling standard deviation.
#' For a window size w (an odd integer between 3 and 2:ncol(X)), the half-window size v is (w-1)/2. The moving window for channel j is equal to
#' max(1,j-v) to min(j+v,ncol(X)). As a result, at the edges of the window are asymmetric at the edges of the spectrum
#'
#' @param X vector or matrix of spectra to correct
#' @param w moving window size (must be an odd integer)
#' @param method Algorithm to calculate PSNV correction (see below for details)
#'
#' @return Matrix with PSNV-corrected spectra
#'
#' @details Two different implementations are included for PSNV: "fast" and "simple"
#' The fast algorithm uses C++ functions to evaluate the rolling average and rolling standard deviation, and is significantly more efficient
#' The simple algorithm calculates the mean and standard deviation of the local window corresponding to each column, and is far slower to calculate
#'
#' @author Cannon Giglio
#' @export
#' @references C. Giglio, J.-M. Roger, E. Andries. Piecewise Standard Normal Variate (PSNV) for scatter correction. J Chemometrics (2025) (Submitted)
#'
#' A. Rinnan, F. van den Berg, SB Engelsen, TrAC Trends in Analytical Chemistry, 2009, 28(10), 1201-1222. (Earliest paper to mention this method)
#' @examples
#' data("marzipan")
#' X.psnv <- psnv(marzipan$X,101) #PSNV with window size of 101 channels
#' matplot(as.vector(marzipan$ax),t(X.psnv),type="l",lty=1)
#'
#'
psnv <- function(X,w,method=c("fast","simple")){
  if(!is.matrix(X)){ X <- as.matrix(X)}
  if(is.null(nrow(X))||ncol(X)==1){ #x is a vector or a column matrix (1 x n)
    X <- matrix(X,nrow=1) #convert to a row matrix
  }
  n <- ncol(X)
  if(w%%2 != 1){stop("w must be an odd positive integer")}
  if(w<=1){stop("w too small")}
  if(w>(2*ncol(X))){stop("w too large")}
  if(missing(method)){method <- "fast"}
  if(method %in% c("fast","simple")==F){stop("Argument 'method' should be one of 'fast','simple'")}
  switch(match.arg(method), fast = {
    rollparams <- movmeanandstd(X,w) #movmeanstd.cpp
    X.psnv <- (X-rollparams$mean)/rollparams$sd
  }, simple = {
    v <- (w-1)/2
    X.psnv <- matrix(0,nrow=nrow(X),ncol=n)
    for(j in 1:n){
      local = max(1,j-v):min(j+v,n)
      xmean <- rowMeans(X[,local])
      xsd <- apply(X[,local],1,sd)
      X.psnv[,j] <- (X[,j]-xmean)/xsd
    }
  })
  return(X.psnv)
}
