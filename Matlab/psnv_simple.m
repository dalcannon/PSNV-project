function Z = psnv_simple( X, w )
%---------------------------------------------------------------------
%
% PSNV - Piecewise SNV 
%
%---------------------------------------------------------------------
%
% PURPOSE: Perform PSNV on a sliding band centered at the jth 
% channel (e.g., wavelength).
%
%---------------------------------------------------------------------
%
% USAGE:
% Z = psnv_simple( X, w );
% 
%---------------------------------------------------------------------
%
% INPUT:
% [1] X: (m,n) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] w: Moving window size (must be an odd integer)
%
%---------------------------------------------------------------------
%
% OUTPUT:
% [1] Z: Matrix of transformed spectra
%
%---------------------------------------------------------------------


%-----------------------------------------------------------
% Sanity checks on spectra
%-----------------------------------------------------------
% Input validation and preprocessing
if ~ismatrix(X)
    error('X must be a matrix');
end
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
% Calculate window half-width
v = (w - 1) / 2;


%-------------------------------------------------
% SNV on sliding window
%-------------------------------------------------
Z = zeros(size(X));
for j = 1:n
    band = max(1,j-v):min(j+v,n);        
    Zband = snv( X(:,band) );   
    Z(:,j) = Zband(:,band==j);     
end

%---------------------------------------------------------------------
end