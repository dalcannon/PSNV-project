function Z = lsnv(X,N)
%---------------------------------------------------------------------
% PURPOSE: Perform local SNV, i.e., SNV on sliding window
%---------------------------------------------------------------------
%
% INPUT
% [1] X: (n,p) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] N: positive integer.  The number of segments bounded 
%     between 2 <= N =< ceil(p/10)
%
%---------------------------------------------------------------------
% 
% OUTPUT
% [1] Z: (m,n) matrix of transformed spectra
% 
%---------------------------------------------------------------------


%-------------------------------------------------
% SANITY CHECK k
%-------------------------------------------------
p = size(X,2); 
if ( N < 2 ) || ( N-floor(N) ~= 0 ) || ( N > floor(p/2) )
    error('v must be integer bounded by 3 <= v <= size(X,2).');
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
% SNV ON DISJOINT WINDOWS
%-------------------------------------------------
Z = zeros(size(X));
for i = 1:N
    window = (sindex(i)+1) : sindex(i+1);   
    Z(:,window) = snv( X(:,window) );    
end

%---------------------------------------------------------------------
end

