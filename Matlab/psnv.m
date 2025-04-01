function Z = psnv( X, v )
%---------------------------------------------------------------------
% PURPOSE: Piecewise SNV, i.e., SNV performed on sliding window.
%---------------------------------------------------------------------
%
% INPUT:
% [1] X: (n,p) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] v: integer bounded by 2 < v < p. 
%     Half-width of sliding window: max(1,j-v):min(j+v,p).
%
%---------------------------------------------------------------------
%
% OUTPUT:
% Z: Matrix of transformed spectra
%
%---------------------------------------------------------------------

%-------------------------------------------------
% SANITY CHECK HALF-WIDTH
%-------------------------------------------------
p = size(X,2);
if (v<3) || (v-floor(v)~=0) || (v>(p-1))
    error('v must be integer bounded by 3 <= v <= size(X,2).');
end    

%-------------------------------------------------
% SNV ON SLIDING WINDOW
%-------------------------------------------------
Z = zeros(size(X));
for j = 1:p
    window = max(1,j-v):min(j+v,p);        
    Zj = snv( X(:,window) );   
    Z(:,j) = Zj(:,window==j);     
end

%---------------------------------------------------------------------
end