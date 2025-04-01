function Z = snv(X)
%---------------------------------------------------------------------
% IN:  X: (m,n) matrix of spectra---each spectrum is aligned row-wise
%---------------------------------------------------------------------
% OUT: X: (m,n) matrix of transformed spectra 
%---------------------------------------------------------------------
A = mean(X,2);                % Sample mean of each row of X
B = std(X,[],2);              % Sample std dev of each row of X
Z = bsxfun( @rdivide, bsxfun(@minus,X,A), B );
%---------------------------------------------------------------------
end