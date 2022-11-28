function Xred = PCsmooth(X,numPCs)
% 
% Churchland et al 2012 'PCsmooth' method
% Smooth individual neural responses by reconstructing each condition's response from main PCs
%
% X      - t x c matrix of responses for individual unit (t = time, c = conditions)
% numPCs - number of PCs to use (default = 6)
%
% TODO: does this require demeaned responses??

if nargin == 1
    numPCs = min(6,size(X,1));
end

[V,score] = pca(X);

Xred = score(:,1:numPCs) * V(:,1:numPCs)';


