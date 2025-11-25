function Z = mpsnv( X, w )
%---------------------------------------------------------------------
%
% MPSNV - Moving Piecewise SNV 
%         (SNV on sliding window without loops)
%
%---------------------------------------------------------------------
%
% PURPOSE: MATLAB implementation of moving SNV using vectorization
%
%---------------------------------------------------------------------
%
% USAGE:
% Z = mpsnv( X, w )
%
%---------------------------------------------------------------------
%
% INPUT
% [1] X: (n,p) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] w: Moving window size (must be an odd integer)
%
%---------------------------------------------------------------------
%
% OUTPUT
% [1] Z: (m,n) matrix of transformed spectra
%
%---------------------------------------------------------------------


%-----------------------------------------------------------
% Sanity checks on spectra 
%-----------------------------------------------------------
% Input validation and preprocessing
if ~ismatrix(X)
    error('X must be a matrix');
end
% Convert to matrix if vector
if isvector(X)
    X = X(:)'; 
end
[~, n] = size(X);

%-----------------------------------------------------------
% Sanity checks on window size 
%-----------------------------------------------------------
if mod(w, 2) ~= 1
    error('w must be an odd positive integer');
end
if w <= 1
    error('w too small');
end
if w > n
    error('w too large');
end


%-----------------------------------------------------------
% Transform spectra
%-----------------------------------------------------------
A = movmean(X,w,2);     % Moving sample mean
B = movstd(X,w,0,2);    % Moving sample standard deviation 
Z = ( X - A ) ./ B;    

%---------------------------------------------------------------------
end