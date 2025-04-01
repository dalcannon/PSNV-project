function Z = pmsc( X, v, r, classic )
%---------------------------------------------------------------------
% PURPOSE: Piecewise MSC, i.e., MSC performed on sliding window.
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
% [4] classic:  logical. (OPTIONAL) 
%     If true, then do classical MSC, otherwise, do an
%     equivalent version formulated by Inge Helland, i.e.,  
%         Z = mean(r) + 1/b*(x-mean(x)).  
%     described in [1] in the mfile msc.m.  The default value is true. 
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
% Sanity check half-width
if (v<3) || (v-floor(v)~=0) || (v>(p-1))
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
% Type of MSC implementation 
if (nargin<4), classic = []; end
if isempty(classic), classic = true; end
if ~islogical(classic)
    error('classic should be logical (true or false)')
end  

%-------------------------------------------------
% MSC ON SLIDING WINDOW
%-------------------------------------------------
Z = zeros(size(X));
for j = 1:p
    window = max(1,j-v):min(j+v,p);        
    Zj = msc( X(:,window), r(window), classic );   
    Z(:,j) = Zj(:,window==j);     
end

%---------------------------------------------------------------------
end