function Z = pemsc( X, v, r, k, classic )
%---------------------------------------------------------------------
% PURPOSE: Piecewise EMSC, i.e., EMSC performed on sliding window.
%---------------------------------------------------------------------
%
% INPUT:
% [1] X: (n,p) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] v: integer bounded by 2 < v < p. 
%     Half-width of sliding window: max(1,j-v):min(j+v,p).
%
% [3] r: p-vector. (OPTIONAL) 
%     The reference spectrum. The default value for r is mean(X,1). 
%
% [4] k: integer: 2 or 3.  (OPTIONAL) 
%     Type of EMSC expansion for a given spectrum x.
%     > If k=2, then x = b*r + a1*e + a2*lambda 
%     > If k=3, then x = b*r + a1*e + a2*lambda + a3*lambda.^2  
%     The default value for k is 3.
%
% [5] classic: logical. (OPTIONAL) 
%     > true: option produces the traditional EMSC transform.
%     > false: option produces the projection-based EMSC transform.
%     The default value is true.
%
%---------------------------------------------------------------------
%
% OUTPUT:
% Z: Matrix of transformed spectra
%
%---------------------------------------------------------------------

%-------------------------------------------------
% DEFAULT INPUT ARGUMENTS
%-------------------------------------------------
p = size(X,2);
% Sanity-check half width
if (v<3) | (v-floor(v)~=0) | (v>(p-1))
    error('v must be integer bounded by 3 <= v <= size(X,2).');
end  
% Get reference spectrum 
if (nargin<3), r = []; end
if isempty(r), r = mean(X,1); end
r = r(:); % Force r to be column vector
p2 = length(r); 
if p2 ~= p
    error('length(r) does not equal size(X,2)'); 
end
% Get k 
if (nargin<4), k = []; end
if isempty(k), k = 3; end
if k~=2 && k~=3
    error('k should be 2 or 3'); 
end    
% Type of EMSC implementation 
if (nargin<5), classic = []; end
if isempty(classic), classic = true; end
if ~islogical(classic)
    error('classic should be logical (true or false)')
end    


%-------------------------------------------------
% EMSC ON SLIDING WINDOW
%-------------------------------------------------
Z = zeros(size(X));
for j = 1:p
    window = max(1,j-v):min(j+v,p);     
    Zj = emsc( X(:,window), r(window), k, classic );    
    Z(:,j) = Zj(:, window==j );          
end

%---------------------------------------------------------------------
end
