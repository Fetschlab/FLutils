function [y,mu] = demean(x,dim)

% first non-singleton dimension
if nargin < 2, dim = find(size(x)~=1,1); end 

mu = nanmean(x,dim);
y  = bsxfun(@minus, x, mu);