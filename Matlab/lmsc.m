function X = lmsc(X,k,r)
%---------------------------------------------------------------------
% PURPOSE: Perform local MSC, i.e., MSC performed separately on k 
% disjoint but continguous spectral segments.
%---------------------------------------------------------------------
% INPUT
% [1] X: (m,n) matrix of spectra---each spectrum is aligned row-wise
% [2] k: number of segments
% [3] r: (n,1) vector: The reference spectrum
%---------------------------------------------------------------------
% OUTPUT
% [1] X: (m,n) matrix of transformed spectra
%---------------------------------------------------------------------
n = size(X,2);                                    
sidx = [0 find(diff( sort(rem(0:(n-1),k)+1))) n];  % Segment indices
for i = 1:k
    band = (sidx(i)+1):sidx(i+1);
    X(:,band) = MSC2( X(:,band), r(band) );  
end
%---------------------------------------------------------------------
end


