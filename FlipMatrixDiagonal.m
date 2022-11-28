function [ Xflip ] = FlipMatrixDiagonal(X)
% [ Xflip ] = FlipMatrixDiagonal(X)
% Averages a square matrix across diagonal

assert(ismatrix(X) && size(X,1)==size(X,2),'Input matrix must be 2D and square');

% ignore main diagonal, this will be maintained in final matrix
BottomHalf = tril(X,0);
TopHalf    = triu(X,0); 

BottomHalfArranged = ( rot90( fliplr( BottomHalf )));           % flip and rotate tril half to overlay triu half
NewTopHalf         = ( BottomHalfArranged + TopHalf) / 2;       % average

Xflip = NewTopHalf + tril( rot90( fliplr( NewTopHalf )), -1 ); % flip back

