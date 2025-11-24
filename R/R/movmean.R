#' Moving Mean
#'
#' Wrapper function for C++ code to evaluate a moving sample mean over each row for a matrix X. Equivalent to MATLAB's 'movstd' function.
#'
#' @param X matrix of data, dimension m x n. If X is a vector, convert to a
#' row matrix using matrix(X,nrow=1)
#' @param w window size, an odd integer. If we define the half-window size as v = (w-1)/2, and ncol(X)=n, then the moving window at channel j will span from
#' max(1,j-v) to min(n,j+v). In other words, the window is symmetric for the middle columns, and is asymmetric at the left and right edges.
#'
#' @return Matrix of data with moving-mean corrected data
#'
#' @examples
#' set.seed(1234)
#' #X <- matrix(sample(-10:10,3*10,replace=T),nrow=3,ncol=10) #generate synthetic data; equivalent to the line below:
#' X <- rbind(c(5,4,-5,-9,4,3,10,-8,-9,9),c(-6,-2,5,-4,3,-7,-3,-7,4,5),c(1,-6,-7,-5,9,-7,9,-6,-3,1))
#' #evaluate a 5-point moving window on each row of X:
#' X_mm <- movmean(X,5)
#'
#' If we load X into MATLAB, we can verify that it is identical to within machine tolerance:
#' x1 = [5,4,-5,-9,4,3,10,-8,-9,9];
#' x2 = [-6,-2,5,-4,3,-7,-3,-7,4,5];
#' x3 = [1,-6,-7,-5,9,-7,9,-6,-3,1];
#' X = [x1;x2;x3];
#' X_mm = movmean(X,5,2); %the third argument indicates to operate on the rows of X
#'
#' @export
movmean <- function(X, w) {
  .Call(`_PSNV_movmean`, X, w)
}
