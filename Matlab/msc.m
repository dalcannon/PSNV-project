function Z = msc(X,r,classic)
%---------------------------------------------------------------------
% PURPOSE: Multiplicatie Scatter Correction
%---------------------------------------------------------------------
%
% INPUT: 
% [1] X: (n,p) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] r: p-vector. (OPTIONAL) 
%     The reference spectrum. The default value for r is mean(X,1). 
%
% [3] classic:  (OPTIONAL) 
%     Logical. If true, then do classical MSC, otherwise, do an
%     equivalent version, i.e.,  
%         Z = mean(r) + 1/b*(x-mean(x)).  
%     described in [1]. The default value is true. 
%
%---------------------------------------------------------------------
% REFERENCE:
% [1] "Related versions of the multiplicative scatter correction 
%     method for preprocessing spectroscopic data". Helland, Naes, 
%     Iasksson. Chemolab, 1995.
%---------------------------------------------------------------------
%
% OUTPUT:
% Z: (n,p) matrix of transformed spectra 
%
%---------------------------------------------------------------------

%-------------------------------------------------
% DEFAULT INPUT ARGUMENTS
%-------------------------------------------------
p = size(X,2);
% Get reference spectrum 
if (nargin<2), r = []; end
if isempty(r), r = mean(X,1); end
r = r(:); % Force r to be column vector
p2 = length(r); 
if p2 ~= p
    error('length(r) does not equal size(X,2)'); 
end  
% Type of MSC implementation 
if (nargin<3), classic = []; end
if isempty(classic), classic = true; end
if ~islogical(classic)
    error('classic should be logical (true or false)')
end 
%-------------------------------------------------
% CLASSICAL OR HELLAND-BASED MSC
%-------------------------------------------------
switch classic
    case true,  Z = C_MSC(X,r,p);
    case false, Z = H_MSC(X,r);
end

%---------------------------------------------------------------------
end

function Z = C_MSC(X,r,p)
%--------------------------------------------------------------------
e = ones(p,1);
M = [r e];
Mpinv = (M'*M) \ M';    % size (2,p)
coef = X*Mpinv';        % size (n,2)
b = coef(:,1);
a = coef(:,2);    
Z = bsxfun( @rdivide, X-a*e' , b );       
%---------------------------------------------------------------------
end



function Z = H_MSC(X,r)
%---------------------------------------------------------------------
% Helland-based MSC in [1]
rmean = mean(r);
Xc = bsxfun( @minus, X, mean(X,2) );
rc = r - rmean;
b  = (Xc*rc)/(rc'*rc);
Z  = rmean + bsxfun( @rdivide, Xc, b );
%---------------------------------------------------------------------
end
