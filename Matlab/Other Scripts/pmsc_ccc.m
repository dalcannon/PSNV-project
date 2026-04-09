function Z = pmsc_ccc( X, w, r, classic )
%---------------------------------------------------------------------
%
% PMSC - Piecewise Multiplicative Scatter Correction
%        Correct Centered Column
%
%---------------------------------------------------------------------
%
% PURPOSE: Perform PSMC on a sliding band centered at the jth 
% channel (e.g., wavelength).  Either the classic or Helland MSC 
% implementation can be selected.
% 
% Like pmsc_simple.m, this implementation computes regression 
% coefficients (slope and intercept) based upon the entire window, 
% but unlike pmsc_simple.m, the correction is only applied to 
% the centered vector X(:,j) as opposed to all of the columns 
% of X in the window done by pmsc_simple.m.
%
%---------------------------------------------------------------------
%
% USAGE
% Z = pmsc_ccc( X, w, r, classic )
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
% [4] classic:  logical. (OPTIONAL) 
%     If true, then do classical MSC, otherwise, do an
%     equivalent version formulated by Helland (see msc.m). 
%     Default value: true. 
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
% Window half-width
v = (w - 1) / 2;


%-----------------------------------------------------------
% Sanity checks on reference mean 
%-----------------------------------------------------------
% Default reference
if (nargin < 3) || isempty(r)
    r = mean(X,1)';
end
% Ensure r is a COLUMN vector
if ~iscolumn(r)
    r = r';
end
% Ensure reference has correct length
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
% Process each channel j
%-----------------------------------------------------------
if classic
    Z = classical_MSC(X,v,r,n);
else
    Z = helland_MSC(X,v,r,n);
end    


%---------------------------------------------------------------------
end


function Z = classical_MSC(X,v,r,n)
%---------------------------------------------------------------------
Z = X;
for j = 1:n
    % Define window boundaries and extract window data
    band = max(1,j-v) : min(n,j+v);
    Xband = X(:,band);     
    rband = r(band);    
    % Classic MSC on windowed data
    M = [ones(length(rband), 1), rband];                
    Mpinv = (M'*M) \ M';               
    c = Xband * Mpinv'; 
    a = c(:,1);
    b = c(:,2);
    % Apply correction only to column j
    Z(:,j) = ( X(:,j) - a ) ./ b;
end
%---------------------------------------------------------------------
end


function Z = helland_MSC(X,v,r,n)
%---------------------------------------------------------------------
Z = X;
for j = 1:n
    % Define window boundaries and extract window data
    band = max(1,j-v) : min(n,j+v);
    Xband = X(:,band);     
    rband = r(band);  
    % Perform Helland MSC on windowed data
    rmean = mean(rband);
    Xmean = mean(Xband,2); 
    Xc = Xband - Xmean;
    %Xc = bsxfun( @minus, Xband, Xmean );
    rc = rband - rmean;
    b  = (Xc*rc)/(rc'*rc);
    % Apply correction only to column j
    Z(:,j) = rmean  +  Xc(:,band==j) ./ b;
end
%---------------------------------------------------------------------
end
