function X = lsnv2(X,k)
%---------------------------------------------------------------------
% PURPOSE: Perform local SNV, i.e., SNV performed separately on k 
% disjoint but contingouos spectral segments.
%---------------------------------------------------------------------
% INPUT
% [1] X: (m,n) matrix of spectra---each spectrum is aligned row-wise
% [2] k: number of segments
%---------------------------------------------------------------------
% OUTPUT
% [1] X: (m,n) matrix of transformed spectra
%---------------------------------------------------------------------
n = size(X,2);                  
dim = 2;                            % Operate on rows (2nd dimension)
sidx = [0 find(diff(sort(rem(0:(n-1),k)+1))) n]; % Segment indices
for i = 1:k
    band = ( sidx(i)+1 ):sidx(i+1); % band of sliding window
    S = X(:,band);                  % Spectral block on X
    A = mean(S,dim);                % Mean of each row of S
    B = std(S,[],dim);              % SD. of each row of S
    X(:,band) = bsxfun(@rdivide, bsxfun(@minus,S,A), B);    
end
%---------------------------------------------------------------------
end

