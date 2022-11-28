
function newd = local_sum ( oldd , sumlen )

%
% function  newd = local_sum(oldd, sumlen)
%
%locally summates the vector oldd
%for instance B=local_sum(A,L) breaks the vectors A into chuncks of length
%L and claulates the sum of each elements in each chunck. The length of B
%will be 1/L of the length of A. L should be an integer.
%
% RK, 2006
%
% modified to work on matrices, treating each column as a vector to be localsummed -CF


newd_len = floor(size(oldd,1)/sumlen);
oldd_len = newd_len * sumlen;
for m = 1:size(oldd,2)    
    % newd = sum(reshape(oldd(1:oldd_len),[sumlen newd_len])',2);
    newd(:,m) = sum(reshape(oldd(1:oldd_len,m),[sumlen newd_len]))';
end