function Z = snv( X )
%---------------------------------------------------------------------
%
% SNV - Standard Normal Variate
%
%---------------------------------------------------------------------
%
% PURPOSE: Perform SNV transformation.  This script uses implicit 
% expansion (used in new versions such as R2016+) as opposed to 
% bsxfun (used in older MATLAB versions such as R2013b).
%
%---------------------------------------------------------------------
%  
% USAGE
% Z = snv( X )
%
%---------------------------------------------------------------------
%
% INPUT
% [1] X: (m,n) matrix of spectra---each spectrum is aligned row-wise
%
%---------------------------------------------------------------------
%
% OUTPUT
% [1] Z: (m,n) matrix of transformed spectra 
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


%-----------------------------------------------------------
% Classic SNV
%-----------------------------------------------------------
A = mean(X,2);                % Sample mean of each row of X
B = std(X,[],2);              % Sample std dev of each row of X
Z = (X - A) ./ B;

%---------------------------------------------------------------------
end