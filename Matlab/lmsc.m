function Z = lmsc(X,N,v,r,method)
%---------------------------------------------------------------------
% PURPOSE: Local MSC, i.e., MSC on sliding window
%---------------------------------------------------------------------
%
% INPUT
% [1] X: (n,p) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] N: positive integer.  The number of segments bounded 
%     between 2 <= k =< ceil(p/10)
%
% [3] v: integer bounded by 2 < v < p. 
%     Half-width of sliding window: max(1,j-v):min(j+v,p).
%     Only applicable for piecewise methods.
%
% [4] r: p-vector. (OPTIONAL) 
%     The reference spectrum. The default value for r is mean(X,1). 
% 
% [5] method: string (OPTIONAL) 
%     The type of MSC variant.  The valid methods are from
%         {'msc','pmsc','mpmsc','emsc','pemsc'}.
%     The piecewise methods are {'pmsc', 'pemsc', 'mpmsc'}.
%     The default value is 'msc'.
%
%---------------------------------------------------------------------
% 
% OUTPUT
% [1] Z: (m,n) matrix of transformed spectra.
% 
%---------------------------------------------------------------------


%-------------------------------------------------
% SANITY CHECK k
%-------------------------------------------------
p = size(X,2); 
% Sanity check for number of segments (N) 
if ( N < 2 ) || ( N-floor(N) ~= 0 ) || ( N > floor(p/2) )
    error('N must be integer bounded by 3 <= N <= size(X,2).');
end 
% Sanity check for window half-width (v). 
% Only meaningful for piecewise methods. 
if ( v < 3 ) || ( v-floor(v) ~= 0 ) || ( v > (p-1) )
    error('v must be integer bounded by 3 <= v <= size(X,2).');
end  
% Get reference spectrum (r)
if (nargin<4), r = []; end
if isempty(r), r = mean(X,1); end
r = r(:)'; % Force r to be column vector
p2 = length(r); 
if p2 ~= p
    error('length(r) does not equal size(X,2)'); 
end 
% Make sure method is valid string
if (nargin<5), method = []; end
if isempty(method), method = 'msc'; end
method = lower(method);
if ~ismember(method,{'msc','pmsc','mpmsc','emsc','pemsc'})
    error('MSC variant you selected is not valid.');
end  


%-------------------------------------------------
% SEGMENT START/STOP INDICES
%-------------------------------------------------
% Make vector of consecutive segment indices: 
% [1...1,2...2,...,N...,N]      
rindex = sort(rem(0:(p-1),N)+1);
% Get indices where change occurs in rindex
tindex = find(diff(rindex));
% Segment indices
sindex = [0 tindex p]; 


%-------------------------------------------------
% CYCLE THROUGH DISJOINT WINDOWS
%-------------------------------------------------
k = 3;          % EMSC expansion type (see emsc.m)
Z = zeros(size(X));
for i = 1:N
    window = (sindex(i)+1) : sindex(i+1);   
    XW = X(:,window);
    rW = r(window);
    switch method
        case 'msc',   XWt = msc(   XW,    rW    );
        case 'pmsc',  XWt = psmc(  XW, v, rW    );
        case 'mpmsc', XWt = mpmsc( XW, v, rW    );  
        case 'emsc',  XWt = emsc(  XW,    rW, k ); 
        case 'pemsc', XWt = pemsc( XW, v, rW, k );
    end  
    Z(:,window) = XWt; 
end

%---------------------------------------------------------------------
end


