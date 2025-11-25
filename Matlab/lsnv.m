function Z = lsnv( X, k )
%---------------------------------------------------------------------
%
% LMSC - Local SNV
%
%---------------------------------------------------------------------
%
% PURPOSE: Chop up the n wavelengths into k disjoint wavelength 
% contiguous bands or segments such that SNV is performed 
% separately on each segment.  
% 
%---------------------------------------------------------------------
%
% USAGE
% >> Z = lsnv( X, k )
%
%---------------------------------------------------------------------
%
% INPUT
% [1] X: (m,n) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] k: The number of disjoint segments ( 2 <= k =< ceil(n/10) )
%
%---------------------------------------------------------------------
% 
% OUTPUT
% [1] Z: (m,n) matrix of transformed spectra.
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
% Sanity checks on number of segments 
%-----------------------------------------------------------
if ( round(k) ~= k ) && ( k > 0 )
    error('w must be an odd positive integer');
end
if k <= 2
    error('w too small');
end
if k > ceil(n/10)
    error('w too large');
end


%-----------------------------------------------------------
% Segment start/stop indices
%-----------------------------------------------------------
% Make vector of consecutive segment indices: 
% [1...1,2...2,...,k...,k]      
consec_index = sort(rem(0:(n-1),k)+1);
% Get indices where change occurs in cindex
change_index = find(diff(consec_index));
% Segment indices
segment_index = [0 change_index n]; 


%-------------------------------------------------
% SNV ON DISJOINT WINDOWS
%-------------------------------------------------
Z = zeros(size(X));
for i = 1:k
    band = (segment_index(i)+1) : segment_index(i+1);
    Z(:,band) = snv( X(:,band) );    
end

%---------------------------------------------------------------------
end

