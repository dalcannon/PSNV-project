#' Fast Multiplicative Scatter Correction (MSC)
#'
#' @description
#' Fast code for evaluation of Multiplicative Scatter Correction (MSC)
#'
#' @param name description
#'
#' @return description
#'
#' @export
#' @examples
#' # example code
#' #Compare with the function msc from PLS package
#' #require(MASS)
#' #require(pls)
#' set.seed(1743)
#'
#' tmp <- matrix(rnorm(10000*1000),nrow=10000,ncol=1000)
#' ref <- colMeans(tmp)
#' system.time(fastmsc(tmp,colMeans(tmp)))
#'

fastmsc <- function(X,r){
  require(MASS)
  # X: (m,n) matrix of spectra---each spectrum is aligned row-wise
  # r: n-vector---reference spectrum
  if(is.data.frame(X)){ X <- as.matrix(X)}
  if(is.null(nrow(X))||ncol(X)==1){ #x is a vector or a column matrix (1 x n)
    X <- matrix(X,nrow=1) #convert to a row matrix
  }
  A = GetPinvM(r)#;           % Pseudoinverse of [ones r]
  #A <- ginv(cbind(rep(1,length(r)),r)) #A = (M'*M) \ M';
  Cmat <- X%*%t(A) ##C = X*A'; % Coefficients across all samples
  Xz <- sweep(X, MARGIN = 1, STATS = Cmat[1, ], FUN = "-", check.margin = FALSE)
  Xz <- sweep(Xz, MARGIN = 1, STATS = Cmat[2, ], FUN = "/", check.margin = FALSE)
}

GetPinvM <- function(r){
#n <- length(r);
#e <- rep(1,n)#ones(n,1);
#M <- cbind(e,r)# [e r]; nx2
M <- cbind(1,r)# [e r]; nx2
A <- ginv(M) #A = (M'*M) \ M';
}

# GetPinvM <- function(r){
#   n <- length(r);
#   e <- rep(1,length(r))#ones(n,1);
#   M <- cbind(rep(1,length(r)),r)# [e r]; nx2
#   A <- ginv(cbind(rep(1,length(r)),r)) #A = (M'*M) \ M';
# }

