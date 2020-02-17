function [Wtilde,sw]=gp_whitening(s,epsilon)
% [Wtilde,sw]=gp_whitening(s,epsilon)
%
% Compute the whitening matrix with data reduction (PCA + sphering).
%
% INPUTS
%          s : signal [samples x electrodes]
%    epsilon : rejected variance (default : 1e-9). Set 0 for no data
%              rejection.
% OUTPUTS
%     Wtilde : "hard"-whitening matrix (with data reduction)
%         sw : whitened data (optionnal, speedup the function if skipped)
%
% Author : Louis Korczowski (louis.korczowski@gmail.com)
if nargin<2,epsilon=1e-9;end
SIGMA=cov(s);
[V D]=eig(SIGMA);
[D,ind] = sort(diag(D),'descend');
V = V(:,ind);%sort spatial eigenvector

if any(D<=0); warning('Rank Deficient Covariance Matrix'); end

% if max(eigWhite1)>indWhite1tmp
%     eigWhite1=1:indWhite1tmp;
% end
D(D<0)=0;
sumD=cumsum(D);%ignore negative eigenvalues in the cumsum
indW=find(sumD>max(sumD)*(1-epsilon),1);
if isempty(indW),indW=size(s,2);end
Dtilde=D(1:indW);
Vtilde=V(:,1:indW); % PCA matrix
Wtilde=diag(1./sqrt(Dtilde + epsilon)) * Vtilde'; %hard whitening matrice
if nargout>1, sw=Wtilde * s';end %output the whitened data
