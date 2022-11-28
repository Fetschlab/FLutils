function Xm = mod2(X,m)


% SJ 08-2022

% adapted version of matlab in-built mod function, where modulo value is in place
% e.g. mod2([1 2 3 4 5 6],3) will return [1 2 3 1 2 3];

% usually find this useful for things like plotting, to determine figure or
% subplot number based on an index in X

% just two lines, but found myself wanting it often enough to make it a
% function...

Xm = mod(X,m); Xm(Xm==0) = m;



