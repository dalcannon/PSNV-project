function Z = mpmsc( X, w, r )
%---------------------------------------------------------------------
%
% MPMSC - Moving Piecewise MSC 
%         (MSC on sliding window without loops)
%
%---------------------------------------------------------------------
%
% PURPOSE: Perform moving PSMC on a sliding band centered at the 
% jth channel (e.g., wavelength) using vectorization (no loops).  
% This script only implements the Helland version of MSC in 
% pmsc_ccc.m.  Moreover, this script uses implicit expansion 
% (used in new versions such as R2016+) as opposed to bsxfun 
% (used in older MATLAB versions such as R2013b).
%
%---------------------------------------------------------------------
%
% USAGE
% Z = mpmsc( X, w, r )
%
%---------------------------------------------------------------------
%
% INPUT
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
%---------------------------------------------------------------------
%
% OUTPUT
% [1] Z: Matrix with PMSC-corrected spectra (same size as X)
%
%---------------------------------------------------------------------

%-----------------------------------------------------------
% Sanity checks on spectra
%-----------------------------------------------------------
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

%-----------------------------------------------------------
% Sanity checks on reference mean 
%-----------------------------------------------------------
% Default reference
if (nargin < 3) || isempty(r)
    r = mean(X, 1)';
end
% Ensure r is a ROW vector
if iscolumn(r)
    r = r';
end
if length(r) ~= n
    error('length(r) must match size(X,2)');
end

%-----------------------------------------------------------
% Pre-compute intermediate moving results
%-----------------------------------------------------------
rmean = movmean(r, w);
rsum  = movsum(r, w);
esum  = movsum(ones(1, n), w);

%-------------------------------------------------
% Moving numerator
%-------------------------------------------------
rsum2 = movsum(r.*r, w);
Numer = rsum2 - 2*rmean.*rsum + rmean.^2.*esum;

%-------------------------------------------------
% Moving denominator
%-------------------------------------------------
Xmean = movmean(X, w, 2);
Xrsum = movsum(X.*r, w, 2);     
Xsum  = movsum(X, w, 2);       
Denom = Xrsum - Xsum.*rmean - Xmean.*(rsum - rmean.*esum);

%-------------------------------------------------
% Compute moving slope and correct spectra
%-------------------------------------------------
Slope = Numer ./ Denom;
Slope(~isfinite(Slope)) = 0;
Z = (X - Xmean).*Slope + rmean;

%---------------------------------------------------------------------
end
