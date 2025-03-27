function Xt = pmsc2( X, w, r )
%---------------------------------------------------------------------
% PURPOSE: Perform piecewise MSC, i.e., MSC performed separately on
% a sliding band centered at the jth wave(length/number)
%---------------------------------------------------------------------
% [1] X: (m,n) matrix of spectra---each spectrum is aligned row-wise
% [2] w: half-width of sliding band of wavelengths/wavenumbers
% [3] r: n-vector---the reference spectrum
%---------------------------------------------------------------------
% [1] Xt: (m,n) matrix of transformed spectra 
%---------------------------------------------------------------------
n = size(X,2);                                                    
Xt = zeros( size(X) );                   
for j = 1:n
    band = max(1,j-w):min(j+w,n);        % wave(length/number) band
    Xmsc = MSC2( X(:,band), r(band) );    % MSC on windowed spectra
    Xt(:,j) = Xmsc(:,band==j);           % Slice out jth wavelength 
end
%---------------------------------------------------------------------
end