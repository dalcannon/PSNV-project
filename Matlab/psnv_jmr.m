function [xsnv]=psnv_jmr(x,w)
%Function to calculate piecewise standard normal variate correction
%x= row vector or matrix of spectra
%w=moving window size (an odd positive integer)
%Output: xsnv = matrix of PSNV-corrected spectra
%Note: This code is significantly slower than the main version, but 
% is compatible with Matlab versions prior to R2016A 
[m,n]=size(x);
v = (w-1)/2;
xsnv = zeros(m,n);
for j=1:n
    local = max(1,j-v):min(j+v,n);
    xsnv(:,j) = (x(:,j)-mean(x(:,local)))./std(x(:,local));
end
end