function Z = mpmsc_v2( X, w, r )
%---------------------------------------------------------------------
%
% MPMSC - Moving Piecewise MSC
%         (C++ MEX version of mpmsc.m)
%
%---------------------------------------------------------------------
%
% PURPOSE: Perform moving PSMC on a sliding band centered at the 
% jth channel (e.g., wavelength) using vectorization (no loops).  
% This scripts outsources all computations to a C++ file (see 
% lines 80-106 of mpmsc.m in C++).  
% 
% This uses a windows-built MEX file.  If you are using linux 
% or MacOS, then the original C++ file is provided such that 
% you can build your own MEX file. Compile the MEX file first: 
% >> mex mpmsc_v2mex.cpp
%
%---------------------------------------------------------------------
%
% USAGE
% Z = mpmsc_v2( X, w, r )
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
% Sanity checks on reference mean 
%-----------------------------------------------------------
% Default reference is column means 
if (nargin < 3) || isempty(r)
    r = mean(X, 1)';
end
% Ensure r is a ROW vector
if iscolumn(r)
    r = r';
end
% Validate inputs
if length(r) ~= n
    error('length(r) must match size(X,2)');
end

%-----------------------------------------------------------
% Call MEX function (implements lines 80-106 of mpmsc.m)
%-----------------------------------------------------------
Z = mpmsc_v2mex(X, r, w);

%---------------------------------------------------------------------
end
