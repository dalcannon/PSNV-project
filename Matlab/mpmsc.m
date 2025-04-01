function Z = mpmsc(X,v,r)
%---------------------------------------------------------------------
% PURPOSE: Moving Piecewise MSC, i.e., piecewise MSC via moving 
% sums and means
%---------------------------------------------------------------------
%
% INPUT
% [1] X: (n,p) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] v: integer bounded by 2 < v < p. 
%     Half-width of sliding window: max(1,j-v):min(j+v,p).
%
% [3] r: p-vector. (OPTIONAL) 
%     The reference spectrum. The default value for r is mean(X,1). 
%
%---------------------------------------------------------------------
%
% OUTPUT
% Z: (m,n) matrix of transformed spectra
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
r = r(:)'; % Force r to be row vector
p2 = length(r); 
if p2 ~= p
    error('length(r) does not equal size(X,2)'); 
end 

%-------------------------------------------------
% PRE-COMPUTE MOVING SUMS AND MEANS
%-------------------------------------------------
w = 2*v+1;     % Maximal width of sliding window 
Xmean = movmean( X, w, 2);
rmean = movmean( r, w ); 
rsum  = movsum( r, w );
esum  = movsum( ones(1,p), w );

%-------------------------------------------------
% COMPUTE THE MOVING SLOPE
%-------------------------------------------------
Denom = movsum( r.*r, w ) - 2*rmean.*rsum + rmean.*(rmean.*esum);
Numer = movsum( bsxfun(@times,X,r), w, 2 ) - ...
        bsxfun( @times, movsum( X, w, 2 ), rmean) - ...
        bsxfun( @times, Xmean, rsum) + ...
        bsxfun( @times, Xmean, rmean.*esum );
Slope = bsxfun(@rdivide,Denom,Numer);

%-------------------------------------------------
% CORRECTED SPECTRA 
%-------------------------------------------------
Z = bsxfun( @plus, (X-Xmean).*Slope, rmean );

%---------------------------------------------------------------------
end


