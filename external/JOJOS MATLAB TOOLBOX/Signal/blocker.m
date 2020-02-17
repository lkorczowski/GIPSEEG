function [Y] = blocker(X, N)
% function [Y] = blocker(X, N)
% This function creates a zero-blocked version of input matrix X, expanding
% it to the size N. Value of possible remaining samples at the end of Y 
% is taken from the last sample of X.
%
% INPUTS
% - X       --> data matrix of size [s x d]
% - N       --> size of output vector Y
%
% OUTPUTS
% - Y       --> Output zero-blocked version of X of size [s x N]
%
% EXAMPLE
% x=randn(5,8), y=blocker(x,12),
%
% HISTORY 
% First version: 13/07/2011
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

blockLength = floor(N/size(X,2));
Y = zeros(size(X,1),N);

for l_ix = 1: size(X,1)
    y = reshape(ones(blockLength,1) * X(l_ix,:),1,[]); 
    Y(l_ix,:) = [y X(l_ix,end)*ones(1,N-length(y))];
end


