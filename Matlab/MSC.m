function X = MSC(X,r)
%---------------------------------------------------------------------
% X: (m,n) matrix of spectra---each spectrum is aligned row-wise
% r: n-vector---reference spectrum 
%---------------------------------------------------------------------
% X: (m,n) matrix of transformed spectra 
%---------------------------------------------------------------------
r = r(:);
A = GetPseudoinverseM(r);           % Pseudoinverse of [ones r]
C = X*A';                           % Coefficients across all samples
X = bsxfun(@minus,  X, C(:,1));     % Mean-center 
X = bsxfun(@rdivide,X, C(:,2));     % Re-scale by slope
%---------------------------------------------------------------------
end

function A = GetPseudoinverseM(r)
%---------------------------------------------------------------------
n = length(r);
e = ones(n,1);
M = [e r];
A = (M'*M) \ M';
%---------------------------------------------------------------------
end
