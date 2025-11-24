mpmsc <- function( X, w, r ){
  if(!is.matrix(X)){ X <- as.matrix(X)}
  if(is.null(nrow(X))||ncol(X)==1){ #x is a vector or a column matrix (1 x n)
    X <- matrix(X,nrow=1) #convert to a row matrix
  }
  if(missing(r)){r <- colMeans(X)}
  n <- ncol(X)# = size(X,2);
  if(length(r)!=n){stop("ref must match number of columns of X")}
  if(w%%2 != 1){stop("w must be an odd positive integer")}
  if(w<=1){stop("w too small")}
  v <- (w-1)/2 #half-window size
  if(w>(2*ncol(X))){stop("w too large")}
  Xcorr <- X
  r <- matrix(r,nrow=1) #Ensure r is a ROW vector
  #Pre-compute intermediate moving results:
  Xmean <- movmean(X,w)
  rmean <- movmean(r,w)
  rsum <- movsum(r,w)
  rsum2 <- movsum(r^2,w)
  esum <- movsum(matrix(1,nrow=1,ncol=n),w)
  #Moving Numerator
  onecol <- matrix(1,nrow=nrow(X))
  Numer <- rsum2 - 2*rmean*rsum + rmean^2*esum
  #Moving Denominator
  Xr <- eachrow(X,r,"*") #X*onecol%*%r #xr <- Rfast::eachrow(X,r,"*")
  Xrsum <- movsum(Xr,w) #VERIFY THIS #  Xrsum = movsum(X.*r, w, 2);
  Xsum <- movsum(X,w)
  Denom <- Xrsum - eachrow(Xsum,rmean,"*") - eachrow(Xmean,rsum-rmean*esum,"*") #Xrsum - Xsum*onecol%*%rmean - Xmean*onecol%*%(rsum-rmean*esum)
  #Compute moving slope and correct spectra
  Slope <- 1/eachrow(Denom,Numer,oper="/")#Numer / Denom
  Slope[is.finite(Slope)==F] <- 0
  Xcorr <- eachrow((X - Xmean)*Slope,rmean,"+") #(X - Xmean)*Slope + onecol%*%rmean #scale(((X - Xmean)*Slope),center=-rmean,scale=F)
  return(Xcorr)
}
