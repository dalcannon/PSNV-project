function Z = pmsc_ccc_v2( X, w, r )
%---------------------------------------------------------------------
%
% PMSC - Piecewise Multiplicative Scatter Correction
%        (C++ MEX version of psmc_ccc)
%
%---------------------------------------------------------------------
%
% PURPOSE: Perform PSMC on a sliding band centered at the jth 
% channel (e.g., wavelength).  This script only implements the 
% Helland version of MSC in pmsc_ccc.m.
% 
% Outsources all computation to a C++ file.  Uses a windows-built 
% MEX file.  If you are using linux or MacOS, then the C++ 
% is provided such that you can build your own MEX file. 
% Compile the MEX file first: 
% >> mex pmsc_ccc_v2mex.cpp
% This implements only the Helland version of MSC in pmsc_ccc.m  
% (i.e., the 'classic' toggle is set to false).
%
%---------------------------------------------------------------------
%
% USAGE
% Z = pmsc_ccc_v2( X, w, r )
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
% [3] r: Reference spectrum (OPTIONAL), and a vector of length n.  
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
% Calculate window half-width
v = (w - 1) / 2;

%-----------------------------------------------------------
% Sanity checks on reference mean 
%-----------------------------------------------------------
% Default reference is column means 
if (nargin < 3) || isempty(r)
    r = mean(X,1)';
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
% Call fast Helland implementation
%-----------------------------------------------------------
Z = pmsc_ccc_v2mex(X, r, v);

%---------------------------------------------------------------------
end
