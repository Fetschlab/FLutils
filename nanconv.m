function C = nanconv(A,B,shape)
%
%   area under convolution kernel B should be 1 to get correct results.
%
%   C = NANCONV(A,B,shape) performs the one dimensional convolution of
%   matrices A and B.
%   The advantage of this function is that unlike the usual convolution
%   where the borders of the convolution result are corrupted by the zero
%   padding used in conv or convn, there will be no sharp changes at the
%   borders. To achieve this goal, instead of calculating the sum over
%   the full length of kernel B at each point along A, only that part of B
%   that falls on A at each point will be used. The function assumes that
%   the area under convolution kernel B is equal to 1 at all times.
%

if (nargin < 3)
  shape = 'full';
end

A = A(:);
B = B(:);

N = convn(A,B,shape);
D = convn(ones(size(A)),B,shape);
C = N./D;

