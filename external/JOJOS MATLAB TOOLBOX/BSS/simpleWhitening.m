function [H, Q, M, W] = simpleWhitening(COV,PERCENT_VAR)
% [H, Q] = simpleWhitening(COV,PERCENT_VAR)
% *************************************************************************
% This function computes whitening matrix (H) from symmetric matrix COV
% so that H*COV*H'=I
%
% Inputs:
% - COV         --> symmetric matrix to whiten (N by N)
% - PERCENT_VAR --> % of the variance (range: [0 1]) to keep out in the
%                   subspace reduction. default: 0 (no subspace reduction)
%
% Outputs:
% - H           --> whitening transformation matrix (M by N)
% - Q       	--> inverse whitening transformation matrix (N by M)
% - M           --> number of output components after data reduction
% - W           --> components eigenvalues (1 by N)
%
%  *** History *** 
% Jonas Chatel-Goldman @ GIPSA-Lab, 24/04/2013
%
% *** Test ***
% Ns=1000; Nc=10; X=rand(Nc,Nc)*rand(Nc,Ns) + .1*randn(Nc,Ns); 
% C=cov(X'); figure,imagescMatrixSelection(C); 
% H=simpleWhitening(C,.01); Cw=H*C*H'; figure,imagescMatrixSelection(Cw);


% check if input matrix is symmetric 
epsilon = 1e-9;
is_sym = max(max(abs(COV - COV'))) <= epsilon;    
if(~is_sym)
    error('simpleWhitening.m --> input matrix is not symmetric.');
end
if nargin<2
   PERCENT_VAR = 0; 
end

% EVD (eigenvalues decomposition)
[U, W] = eig(COV); 
[W, sort_ix] = sort(diag(W),'descend'); % sort eigenvalues
U = U(:,sort_ix);                       % sort eigenvectors accordingly
H = diag(1./W.^(1/2)) * U'; 
Q = U * diag(W.^(1/2));

% data reduction
M = find(cumsum(W)/sum(W)>=(1-PERCENT_VAR),1,'first');
H = H(1:M,:);
Q = Q(:,1:M);

% H = cov(X)^(-1/2); % this is the simplest way to perform whitening on the 
% data covariance matrix, see
% http://en.wikipedia.org/wiki/Whitening_transformation


