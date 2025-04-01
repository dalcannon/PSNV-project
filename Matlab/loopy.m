function [Z,RRSSQ] = loopy( X, v, N, method )
%---------------------------------------------------------------------
% PURPOSE: Perform Loopy XX-MSC
%---------------------------------------------------------------------
%
% INPUT:
% [1] X: (n,p) matrix of spectra.  Each spectrum in the matrix X 
%     is aligned row-wise.
%
% [2] v: integer bounded by 2 < v < p. 
%     Half-width of sliding window: max(1,j-v):min(j+v,p).
%    
%
% [3] N: integer bounded by 2 < v < p
%     Number of loopy iterations bounded by 1 < N <= 10
%     Only valid for piecewise methods:  
%
% [4] method: string (OPTIONAL) 
%     The type of MSC variant.  The valid methods are from
%         {'msc','pmsc','mpmsc','emsc','pemsc'}.
%     The piecewise methods are {'pmsc', 'pemsc', 'mpmsc'}.
%     The default value is 'msc'
%
%---------------------------------------------------------------------
% 
% INPUT:
% [1] Z: (N+1,1) cell array. Z{1} = Z, Z{2} is transformed spectra 
%     of Z{1}, Z{3} is transformed spectra of Z{2}, and so on.
% [2] RRSSQ: (N,1) vector.  RRSSQ(i) is the relative root sum of 
%     squares: sum(sum( (Z{i+1}-Z{i}).^2 )) / sum(sum( Z{i}.^2 ));
%      
%---------------------------------------------------------------------
%
% REFERENCE: 
% [1] Loopy MSC: A Simple Way to Improve Multiplicative Scatter
%     Correction. WILLEM WINDIG, JEREMY SHAVER, and RASMUS BRO.
%     APPLIED SPECTROSCOPY, 2008
% 
%---------------------------------------------------------------------

               
%-------------------------------------------------
% DEFAULT INPUT ARGUMENTS
%-------------------------------------------------
p = size(X,2);
% Sanity check for window half-width (v). 
% Only meaningful for piecewise methods. 
if (v<3) || (v-floor(v)~=0) || (v>(p-1))
    error('v must be integer bounded by 3 <= v <= size(X,2).');
end  
% Sanity check number of loopy iterations (N)
if (N<2) || (N-floor(N)~=0) || (N>10)
    error('N must be integer bounded by 2 <= N <= 10.');
end  
% Make sure method is valid string
if (nargin<4), method = []; end
if isempty(method), method = 'msc'; end
method = lower(method);
if ~ismember(method,{'msc','pmsc','mpmsc','emsc','pemsc'})
    error('MSC variant you selected is not valid.');
end    


%-------------------------------------------------
% INITIALIZE z
%-------------------------------------------------
Z = cell(N+1,1);
Z{1} = X;
for i = 2:(N+1), Z{i} = zeros(size(X)); end

%-------------------------------------------------
% CYCLE THROUGH LOPY ITERATIONS
%-------------------------------------------------
k = 3;              % EMSC expansion type (see emsc.m)
RRSSQ = zeros(N,1); % Initialize RRSSQ
for i = 1:N
    Xold = Z{i};
    r = mean(Xold,1);
    switch method
        case 'msc',   Xnew = msc(   Xold,    r    );
        case 'pmsc',  Xnew = psmc(  Xold, v, r    );
        case 'mpmsc', Xnew = mpmsc( Xold, v, r    );  
        case 'emsc',  Xnew = emsc(  Xold,    r, k ); 
        case 'pemsc', Xnew = pemsc( Xold, v, r, k );
    end    
    Z{i+1} = Xnew;
    RRSSQ(i) = sum(sum( (Xnew-Xold).^2 )) / sum(sum( Xold.^2 ));
end

%---------------------------------------------------------------------
end


