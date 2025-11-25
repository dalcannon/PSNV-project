function Z = msc( X, r, classic )
%---------------------------------------------------------------------
%
% MSC - Multiplicative Scatter Correction
%
%---------------------------------------------------------------------
%
% PURPOSE: There are two equivalent MSC implementations that can 
% be performed: 
% a) The classic version where regression coefficients (slope and 
%    intercept) arfe compuited and the correction is applied to 
%    all columns of X.
% b) The implementation done by Helland et al. is decribed in [1].  
%    In short, no least squares procedure is needed to compute the
%    slope (b) and the correction on x is described by 
%        z = mean(r) + 1/b*(x-mean(x)).  
%
%---------------------------------------------------------------------
%
% USAGE 
% >> Z = msc( X, r, classic )
%
%---------------------------------------------------------------------
% INPUT 
% [1] X: (m,n) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] r: n-vector. (OPTIONAL) 
%     The reference spectrum. The default value for r is mean(X,1). 
%
% [3] classic:  logical. (OPTIONAL) 
%     If true, then do classical MSC, otherwise, do an
%     equivalent version formulated by Helland.
%     Default value: true. 
%
%---------------------------------------------------------------------
%
% OUTPUT
% [1] Z: (m,n) matrix of transformed spectra 
%
%---------------------------------------------------------------------
%
% REFERENCE:
% [1] "Related versions of the multiplicative scatter correction 
%     method for preprocessing spectroscopic data". Helland, Naes, 
%     Iasksson. Chemolab, 1995.
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
% Sanity checks on reference mean 
%-----------------------------------------------------------
% Default reference is column means 
if (nargin < 2) || isempty(r)
    r = mean(X,1);
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
if (nargin < 3) || isempty(classic)
    classic = true;
end    
if ~islogical(classic)
    error('classic should be logical (true or false)')
end  


%-------------------------------------------------
% CLASSICAL OR HELLAND-BASED MSC
%-------------------------------------------------
if classic
    Z = classic_MSC(X,r);
else    
    Z = helland_MSC(X,r);
end

%---------------------------------------------------------------------
end

function Z = classic_MSC(X,r)
%--------------------------------------------------------------------
M = [ones(length(r),1) r];
Mpinv = (M'*M) \ M';   
c = X*Mpinv';          
a = c(:,1);
b = c(:,2);    
Z = bsxfun( @rdivide, bsxfun(@minus,X,a), b );       
%---------------------------------------------------------------------
end


function Z = helland_MSC(X,r)
%---------------------------------------------------------------------
rmean = mean(r);
Xmean = mean(X,2);  
rc = r - rmean;
Xc = bsxfun( @minus, X, Xmean );
b  = (Xc*rc)/(rc'*rc);
Z  = rmean + bsxfun( @rdivide, Xc, b );
%---------------------------------------------------------------------
end
