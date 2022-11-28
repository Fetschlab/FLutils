function [longest,pos] = findones(x)
% [longest,pos] = findones(x)
% return length and starting index of longest series of 'ones' (true
% values) in a logical vector

x = x(:);
[longest,pos] = max(accumarray(nonzeros((cumsum(~x)+1).*x),1));