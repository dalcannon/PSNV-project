function Z = emsc( X, r, k, classic )
%---------------------------------------------------------------------
% PURPOSE: Extended Multiplicatie Scatter Correction
%---------------------------------------------------------------------
%
% INPUT: 
% [1] X: (n,p) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] r: p-vector. (OPTIONAL) 
%     The reference spectrum. The default value for r is mean(X,1). 
%
% [3] k: integer: 2 or 3.  (OPTIONAL) 
%     Type of EMSC expansion for a given spectrum x.
%     > If k=2, then x = b*r + a1*e + a2*lambda 
%     > If k=3, then x = b*r + a1*e + a2*lambda + a3*lambda.^2  
%     The default value for k is 3.
%
% [4] classic: logical. (OPTIONAL) 
%     > true: option produces the traditional EMSC transform.
%     > false: option produces the projection-based EMSC transform.
%     The default value is true.

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
% Get k 
if (nargin<3), k = []; end
if isempty(k), k = 3; end
if k~=2 && k~=3
    error('k should be 2 or 3'); 
end    
% Type of EMSC implementation 
if (nargin<4), classic = []; end
if isempty(classic), classic = true; end
if ~islogical(classic)
    error('classic should be logical (true or false)')
end   


%-------------------------------------------------
% EMSC BASIS VECTORS (excluding r)
%-------------------------------------------------
e = ones(p,1);
lambda = linspace(-1,1,p)'; 
switch k
    case 2, N = [ e  lambda  ];
    case 3, N = [ e  lambda  lambda.^2 ];
end

%-------------------------------------------------
% CLASSICAL OR PROJECTION-BASED EMSC
%-------------------------------------------------
switch classic
    case true,  Z = C_MSC(X,r,N);
    case false, Z = P_MSC(X,r,N);
end    

%---------------------------------------------------------------------
end



function Z = C_MSC(X,r,N)
%---------------------------------------------------------------------
M = [r N];
Mpinv = (M'*M) \ M';  % size (k,p)
coef = X*Mpinv';      % size (n,k)
b = coef(:,1);
a = coef(:,2:end);
Z = bsxfun( @rdivide, X-a*N', b );
%---------------------------------------------------------------------
end



function Z = P_MSC(X,r,N)
%---------------------------------------------------------------------
% Get Q, the orthonormal basis vectors of N
[Q,~] = qr(N,'econ'); 
% Project: P*r and P*X where P = I-Q*Q' projects onto nullspace of N
Pr = r  - Q*(Q'*r );
PX = X' - Q*(Q'*X');
% Compute least-squares estimates for slope coefficients
b = (PX'*Pr) / (Pr'*Pr);
% Transform
Z = bsxfun( @rdivide, PX', b );
%---------------------------------------------------------------------
end

