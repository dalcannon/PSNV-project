function Z = lmsc( X, k, r, classic )
%---------------------------------------------------------------------
%
% LMSC - Local MSC 
%
%---------------------------------------------------------------------
%
% PURPOSE: Chop up the n wavelengths into k disjoint wavelength 
% contiguous bands or segments such that MSC is performed 
% separately on each segment.  Either the classic or Helland MSC 
% implementation can be selected.
% 
%---------------------------------------------------------------------
%
% USAGE
% >> Z = lmsc( X, k, r, classic )
%
%---------------------------------------------------------------------
%
% INPUT
% [1] X: (m,n) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] k: The number of disjoint segments ( 2 <= k =< ceil(n/10) )
%
% [3] r: p-vector. (OPTIONAL) 
%     The reference spectrum. The default value for r is mean(X,1). 
%
% [4] classic:  logical. (OPTIONAL) 
%     If true, then do classical MSC, otherwise, do an
%     equivalent version formulated by Helland (see msc.m). 
%     Default value: true. 
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
% Segment start/stop indices
%-----------------------------------------------------------
% Make vector of consecutive segment indices: 
% [1...1,2...2,...,k...,k]      
consec_index = sort(rem(0:(n-1),k)+1);
% Get indices where change occurs in cindex
change_index = find(diff(consec_index));
% Segment indices
segment_index = [0 change_index n]; 


%-----------------------------------------------------------
% Perform MSC on each of the k disjoint segments 
%-----------------------------------------------------------
Z = zeros(size(X));
for i = 1:k
    band = (segment_index(i)+1) : segment_index(i+1); 
    Z(:,band) = msc( X(:,band), r(band), classic ); 
end

%---------------------------------------------------------------------
end


