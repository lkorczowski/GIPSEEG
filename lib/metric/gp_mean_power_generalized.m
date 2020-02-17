function [P,conv,niter] = mean_power_generalized(COV,p,i_cfg)
%
% Generalized Power Mean of Symetrical Positive Definite Matrices (SPD) using
% fixed-point algorithm proposed by [1]
% ------ INPUTS: 
%     COV: the NxNxK array of K SPD matrices
%       p: the power present in the definition of power mean
%
% i_cfg is a structure containing the optional parameters
%          tol: the tolerance for the convergence criterion
%    N_itermax: maximum number of iterations
%           P0: the initial guess for mean
%
% ------ OUTPUTS: 
%       P: the power mean
%    conv: the stopping criterion
%   niter: number of iteration required
%
% *** History: 07-March-2016
% *** Author: Paolo Zanini and Louis Korczowski @Gipsa-lab 2016
% *** Reference: 
% [1] Marco Congedo et. al, 2016 
%   "A Fixed-Point Algorithm for Estimating Power Means 
%   of Positive Definite Matrices"


%load parameters

K = size(COV,3);
N = size(COV,1);

cfg = struct('tol',10^(-10),'N_itermax',100,'P0',P0,'w',ones(K,1));

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

% normalize weights (L1 norm)
cfg.w=abs(cfg.w);cfg.w=cfg.w/sum(cfg.w);

%inititalization P
P0 = 0;
for indN=1:K
    P0 = P0 + w(indN)* gp_powm(COV(:,:,indN),p); % eq. (13)
end
P0 = gp_powm(P0,1/p); % eq. (13)

niter=0;

if sign(p)>0
X = gp_powm(P0,-1/2);%eye(N);
else
    X = gp_powm(P0,1/2);%eye(N);
end

while (niter<cfg.N_itermax)
    niter = niter + 1;
    
    H = 0;
    for indN=1:K
        H = H + w(indN)*gp_powm(X*gp_powm(COV(:,:,indN),sign(p))*(X'),abs(p));
    end
%     H = H/N;
    phi= 0.375/abs(p);
    
    Ptmp = H^(-phi)*X;
    
    %     conv = distance(Ptmp/p,X/p, 'euclid');
    conv(niter)=norm(H-eye(N),'fro')/sqrt(N);
    X = Ptmp;
    
    if conv(niter)<cfg.tol % break if the improvement is below the tolerance
        break;
    end
end

if sign(p)<0
    Y =X;
else
    Y=pinv(X)';
end

P=Y'*Y;

if niter >= cfg.N_itermax
    warning(['Can not converge, stopped at iteration ' num2str(cfg.N_itermax)]);
end


    function B=gp_powm(A,alpha)
        % B=gp_powm(A,alpha)
        % Fast matrix power using svd
        %
        % *** History: 04-March-2016
        % *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2016
        %
        % see also: powm
        
        if abs(alpha)<1
            [U,S,V]=svd(A);
            B=U*diag((diag(S)).^(alpha))*V';
        elseif alpha==1
            B=A;
        else
            B=A^(alpha);
        end
        
    end
end