function Z = mpsnv( X, v )
%---------------------------------------------------------------------
% PURPOSE: Moving Piecewise SNV, i.e., SNV on sliding window without 
% the need for loops
%---------------------------------------------------------------------
%
% INPUT
% [1] X: (n,p) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] v: integer bounded by 2 < v < p. 
%     Half-width of sliding window: max(1,j-v):min(j+v,p).
%
%---------------------------------------------------------------------
%
% OUTPUT
% Z: (m,n) matrix of transformed spectra
%
%---------------------------------------------------------------------
% Sanity check half-width
p = size(X,2);
if (v<3) || (v-floor(v)~=0) || (v>(p-1))
    error('v must be integer bounded by 3 <= v <= size(X,2).');
end   
% Moving stats: sample mean (A) and sample standard deviation (B)
w = 2*v+1;                                 % Maximal window width
A = movmean(X,w,2);     
B = movstd(X,w,0,2);    
% Transform spectra
Z = ( X - A ) ./ B;    
%---------------------------------------------------------------------
end