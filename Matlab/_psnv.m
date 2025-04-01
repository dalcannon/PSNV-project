function X = psnv( X, w )
%---------------------------------------------------------------------
% PURPOSE: Perform piecewise SNV, i.e., SNV performed separately on
% a sliding band centered at the jth wave(length/number)
%---------------------------------------------------------------------
% INPUT
% [1] X: (m,n) matrix of spectra---each spectrum is aligned row-wise
% [2] w: width of sliding window (Must be an odd positive integer)
%---------------------------------------------------------------------
% OUTPUT
% [1] X: (m,n) matrix of transformed spectra
%---------------------------------------------------------------------
dim = 2;                 % Operate on rows (2nd dimension)
A = movmean(X,w,dim);   % Moving sample mean
B = movstd(X,w,0,dim);  % Moving sample standard deviation
X = ( X - A ) ./ B;      % Transform spectra
% Executable statments above could be replaced with the one line below
% X = (X-movmean(X,2*w+1,2))./movstd(X,2*W+1,0,2);
%---------------------------------------------------------------------
end
