#' Localized Standard Normal Variate (LSNV)
#'
#' @description
#' Perform localized standard normal variate (LSNV), i.e., SNV performed separately on k disjoint but continuous spectral segments.
#'
#' @param X row vector or matrix of spectra to correct
#' @param k number of segments
#'
#' @return Matrix of LSNV-corrected spectra
#'
#' @export
#' @examples
#' #data("marzipan")
#' #X.lsnv <- lsnv(marzipan$X,5) #LSNV with 5 segments
#' #matplot(as.vector(marzipan$ax),t(X.lsnv),type="l",lty=1)
#' @references C. Giglio, J.-M. Roger, E. Andries. Piecewise Standard Normal Variate (PSNV) for scatter correction. J Chemometrics (2025) (Submitted)
#'
#' Y. Bi, K. Yuan, W. Xiao, J. Wu, C. Shi, J. Xia, G. Chu, G. Zhang, G. Zhou, Analytica chimica acta 2016, 909, 30–40.
#'
lsnv <- function(X,k){
  if(!is.matrix(X)){ X <- as.matrix(X)}
  if(is.null(nrow(X))||ncol(X)==1){ #x is a vector or a column matrix (1 x n)
    X <- matrix(X,nrow=1) #convert to a row matrix
  }
  n <- ncol(X) #n = size(X,2);
  sidx <- c(0,which(diff(sort(((0:(n-1))%%k)+1))==1),n) # % Segment indices
  Xnew <- X
  for(i in 1:k){
    band <- (sidx[i]+1):sidx[i+1]
    S <- X[,band];
    A <- rowMeans(S);
    B <- apply(S,1,sd)
    Xnew[,band] <- (S-A)/B
  }
  return(Xnew)
}
