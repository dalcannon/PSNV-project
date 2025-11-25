function Z = pmsc_simple( X, w, r, classic )
%---------------------------------------------------------------------
%
% PMSC - Piecewise Multiplicative Scatter Correction
%        (the old school slow way)
%
%---------------------------------------------------------------------
%
% PURPOSE: Perform PSMC on a sliding band centered at the jth 
% channel (e.g., wavelength).  Either the classic or Helland MSC 
% implementation can be selected.
%
% Since an external function msc.m is used, this implementation 
% computes regression coefficients based upon the entire window, 
% AND corrects all of the columns of X in the window (although 
% only the center column is actually extracted in the outer loop).  
%
%---------------------------------------------------------------------
%
% USAGE:
% Z = pmsc_simple( X, w, r, classic )
%
%---------------------------------------------------------------------
%
% INPUT:
% [1] X: Matrix of spectra to correct (m x n) where each row is 
%     a spectrum where m is the number of spectra and n is the 
%     number of channels.
%
% [2] w: Moving window size (must be an odd integer)
%
% [3] r: Reference spectrum (OPTIONAL)
%     Must be a vector of length n.  
%     Default value: vector of column means of X, i.e., r=mean(X,1).
%
% [4] classic:  logical. (OPTIONAL) 
%     If true, then do classical MSC, otherwise, do an
%     equivalent version formulated by Helland (see msc.m). 
%     Default value: true. 
%
%---------------------------------------------------------------------
%
% OUTPUT:
% [1] Z: Matrix with PMSC-corrected spectra (same size as X)
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

%-----------------------------------------------------------
% Sanity checks on reference mean
%-----------------------------------------------------------
% Default reference is column means 
if (nargin < 3) || isempty(r)
    r = mean(X, 1);
end
% Ensure r is a COLUMN vector
if ~iscolumn(r)
    r = r';
end
% Validate inputs
if length(r) ~= n
    error('length(r) must match size(X,2)');
end

%-----------------------------------------------------------
% Sanity checks on implementation logical 
%-----------------------------------------------------------
if (nargin < 4) || isempty(classic)
    classic = true;
end    
if ~islogical(classic)
    error('classic should be logical (true or false)')
end  

%-----------------------------------------------------------
% Perform full MSC on sliding window
%-----------------------------------------------------------
Z = zeros(size(X));
for j = 1:n
    band = max(1,j-v):min(j+v,n);        
    Zband = msc( X(:,band), r(band), classic );   
    Z(:,j) = Zband(:,band==j);     
end

%---------------------------------------------------------------------
end


