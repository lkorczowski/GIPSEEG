function [W,e,niter] = mean_wasserstein(COV,i_cfg)

% Function to estimate the Wasserstein mean for a set of N SPD matrices
% of dimension p x p gathered in the pxpxN array COV.
% The optional parameters (need to be a struct) to be set are:
% - method: the algorithm used to estimate the mean.
% - tol: the tolerance for the convergence criterion
% - N_itermax: maximum number of iterations
% - W: the initial guess for mean
% The two implemented method are:
% fp_old: point fixed algorithm based on (W^(1/2) C_i W^(1/2))^(1/2) (by Alexandre)
% fp_new: point fixed algorithm based on (W C_i)^(1/2) (by Paolo)

N = size(COV,3);
p = size(COV,1);

Wdef = 0;
for i=1:N
    Wdef = Wdef + COV(:,:,i)^(0.5);
end
Wdef = (Wdef/N)^2;

cfg = struct('method','fp_new','tol',10^(-6),'N_itermax',200,'W',Wdef);

% read the acceptable names
optionNames = fieldnames(cfg);

if  nargin<2 || isempty(i_cfg)
        warning('default cfg')
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

niter = 0;

if strcmp(cfg.method,'fp_old')
    
    K = cfg.W^(0.5);

    while (niter<cfg.N_itermax)
        niter = niter+1;
   
        Ktmp = 0;
        for i=1:N
            Ktmp = Ktmp + (K*COV(:,:,i)*K)^(0.5); 
        end
        Ktmp = sqrt(1/N)*(Ktmp)^(0.5);
    
        conv = distance(Ktmp/p,K/p,'euclid');
        K = Ktmp;
    
        e(niter) = conv;
        if conv<cfg.tol % break if the improvement is below the tolerance
            break; 
        end
    end

W = K^2;

elseif strcmp(cfg.method,'fp_new')
    
    W=cfg.W;
    
    while (niter<cfg.N_itermax)
        niter = niter + 1;
    
        Wtmp = 0;
        for i=1:N
            Wtmp = Wtmp + (W*COV(:,:,i))^(0.5);
        end
        Wtmp = 0.5*(Wtmp/N + transpose(Wtmp)/N);
        
        conv =  distance(Wtmp/p,W/p,'euclid');
        W=Wtmp;
        
        e(niter) = conv;
        if conv<cfg.tol % break if the improvement is below the tolerance
            break; 
        end
    end
    
    
else 
    error('You have to choose an admissible method');
end

if niter==cfg.N_itermax
    disp('Warning : Nombre d''iterations maximum atteint');
end

end
    
    