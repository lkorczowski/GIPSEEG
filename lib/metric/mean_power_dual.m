function [W,e,niter] = mean_power_dual(COV,i_cfg)

% Function to estimate the power mean (Lawson and Lim, 2013)
% and its dual version, using the algorithms described by Marco
% in his pdf tutorials.
% The mandatory parameters are:
% - COV: the p x p x N array of N SPD matrices
% The optional parameters (need to be a struct) to be set are:
% - tol: the tolerance for the convergence criterion
% - N_itermax: maximum number of iterations
% - W: the initial guess for mean
% - t: the power present in the definition of power mean
% - phi: the power used in the algorithm with the matrix H
%
% Author : Paolo Zanini @Gipsa-lab 2016, modified by Louis Korczowski. Based on the work of Marco Congedo

N = size(COV,3);
p = size(COV,1);

Wdef = 0;
% init
COVinv = COV;
for i=1:N
    COVinv(:,:,i)=inv(COV(:,:,i));
end
for i=1:N
    Wdef = Wdef + COVinv(:,:,i)^(0.5);
end
Wdef = (Wdef/N)^-2;

% Wdef = 0;
% for i=1:N
%     Wdef = Wdef + COV(:,:,i)^(0.5);
% end
% Wdef = (Wdef/N)^2;

cfg = struct('tol',10^(-3),'N_itermax',200,'W',Wdef,'t',0.5,'phi',0.8);

% read the acceptable names
optionNames = fieldnames(cfg);


if  nargin<3 || isempty(i_cfg)
%     warning('default cfg')
else
    nameArgs=fieldnames(i_cfg);
    nArgs = length(nameArgs);
    
    for indO = 1: nArgs
        inpName = nameArgs{indO}; %case sensitive
        if any(strcmp(inpName,optionNames))
            cfg.(inpName) = i_cfg.(inpName);
        else
            warning('%s is not a recognized parameter name',inpName)
        end
    end
end

niter=0;


P=cfg.W;

while (niter<cfg.N_itermax)
    niter = niter + 1;
    
    H = 0;
    %         COVinv = COV;
    %         for i=1:N
    %             COVinv(:,:,i)=inv(COV(:,:,i));
    %         end
    for i=1:N
        H = H + (P*COVinv(:,:,i)*transpose(P))^(cfg.t);
    end
    H = H/N;
    Ptmp = H^(-cfg.phi)*P;
    
    conv = distance(Ptmp/p,P/p, 'euclid');
    P = Ptmp;
    
    e(niter) = conv;
    if conv < cfg.tol % break if the improvement is below the tolerance
        break;
    end
end

W = transpose(P)*P;

if niter >= cfg.N_itermax
    warning(['Can not converge, stopped at iteration ' num2str(cfg.N_itermax)]);
end
end