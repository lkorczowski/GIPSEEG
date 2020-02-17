function [B,P,CRIT] = JBSS_li_procr(C, S, eps, max_it, B_init)
% function [B,P,CRIT] = JBSS_li_procr(C, S, eps, max_it, B_init)
%
% ************************************************************************
% This function performs Joint Blind Source Separation (JBSS) using an
% implementation of Li's procruste algorithm for Generalized Joint
% Approximate Diagonalization
% 
% It produces a mixing matrix B (M components * N channels) 
% such as:   C_out(k) = B.C(k).B' has a "stripe-diagonal" form
%
%
% Input:
% ------
%   C   :   3D array containing the different matrices, may be given in
%           unfolded version (channels*channels*cov_struct)
%   S   :   the size of the square submatrices D_ii, i.e. S is a vector
%   eps : (optional)criterion threshold
%   max_it:(optional) maximum iteration number
%   B_init: (optional) initialization for demixing matrix 
%
% Output:
% -------
%   B   : demixing matrix whose product B*C*B' is stripe-diagonal
%   P   : the stripes diagonalised version of input P
%   CRIT: stopping criterion, based on off-diag of P matrices
%               CRIT(1,:) absolute change criterion (convergence)
%               CRIT(2,:) sign (-1 = converge, +1 = diverge)
%
% History:
% --------
% *** 2011-04-21
% Created by J. Chatel-Goldman @ GIPSA Lab, based on work of:
% X.-L. Li, T. Adal?, and M. Anderson, "Joint blind source separation by 
% generalized joint diagonalization of cumulant matrices,"
% Signal Processing, vol. 91, no. 10, pp. 2314–2322, 2011.

% assertion on inputs
if sum(S)~=size(C,1) error('[JBSS_li_proc] Incoherent submatrices size description'); end
if find(S-S(1)) error('[JBSS_li_proc] Input matrices should have the same size'); end
if length(S)<=1 error('[JBSS_li_proc] Only one set is given as input!'); end
S = assert_col(S); % verify if S is a column vector as it should be
if nargin < 4   max_it  = 5000;     end
if nargin < 3   eps = 1e-12;        end
iter    = 1;

% initialization
I   = length(S);            % number of subjects
N   = sum(S);               % total number of channels
Cix = [cumsum(S)-S+1, cumsum(S)]; % index of Cij matrices 
S   = S(1);                 % same size fo each subject
K   = size(C,3);            % number of covariance structures
Lcross = S*(I-1);           % length of cross structures
% Fix = [kron([0:I*I*K-1]',S)+1 , kron([1:I*I*K]',S)]; % index of E matrices
PI_diag     = zeros(1,max_it);  % init stripe-diags Performance Index
CRIT  = zeros(3,max_it);
% CRIT(1,:) --> convergence absolute value
% CRIT(2,:) --> convergence absolute change value
% CRIT(3,:) --> convergence absolute change sign
if (nargin < 5) || isempty(B_init)   
    B = eye(N); 
else
    B = B_init;         
end

% memory allocation   
D       = zeros(N,N,K);
P       = zeros(N,N,K);
F       = zeros(2*K*S*(I-1),S,I);
E       = zeros(2*K*S*(I-1),S,I);

% stack C(i,j,k)-->E. Each E(:,:,p) is constructed for one specific set 'p'
for p = 1:I 
    for k = 1:K 
        off = (k-1)*2*Lcross ;
        E(off+1:off+Lcross,:,p) = [C(Cix(p,1):Cix(p,2),1:Cix(p,1)-1,k) C(Cix(p,1):Cix(p,2),Cix(p,2)+1:end,k)]' ;
        E(off+Lcross+1:off+2*Lcross,:,p) = [C(1:Cix(p,1)-1,Cix(p,1):Cix(p,2),k) ; C(Cix(p,2)+1:end,Cix(p,1):Cix(p,2),k)] ;
    end
end

% Performance Indexes calculation 
PI_diag(iter)   = DIAG_crit_proc(); 


% ******************************* MAIN LOOP *******************************
while iter<=1 || ((CRIT(2,iter)>eps && iter<max_it))
    % UI display
    fprintf('%s','.'); if(mod(iter,50) == 0) fprintf('\n'); end
    %    disp(['Processing Li-JBBS, iteration #' int2str(iter)]);
    iter = iter+1; 
    
    % process D (only keep diag part of each matrix set)
    for k = 1:K
        for i = 1:I 
            for j = 1:I 
                D(Cix(i,1):Cix(i,2),Cix(j,1):Cix(j,2),k) = ... 
                    diag(diag(  B(Cix(i,1):Cix(i,2),Cix(i,1):Cix(i,2))  * ...
                                C(Cix(i,1):Cix(i,2),Cix(j,1):Cix(j,2),k) * ...
                                B(Cix(j,1):Cix(j,2),Cix(j,1):Cix(j,2))' )) ;
            end
        end
    end

    for p = 1:I 
        % stack D(,,k)*B(p) --> F(,,p) 
        for k = 1:K 
            off = (k-1)*2*Lcross ;
            index = 1;
            for i = setdiff(1:I,p)%[(1:p-1) (p+1:I)] 
                F(off+(index-1)*S+1:off+index*S,:,p) =  ...
                    B(Cix(i,1):Cix(i,2),Cix(i,1):Cix(i,2))' * D(Cix(p,1):Cix(p,2),Cix(i,1):Cix(i,2),k); 
                F(off+Lcross+(index-1)*S+1:off+Lcross+index*S,:,p) = ...
                    B(Cix(i,1):Cix(i,2),Cix(i,1):Cix(i,2))' * D(Cix(i,1):Cix(i,2),Cix(p,1):Cix(p,2),k);
                index = index+1;
            end
        end
        
        % perform SVD on (F'*E)
        [U,ddd,V] = svd(F(:,:,p)'*E(:,:,p));
        B(Cix(p,1):Cix(p,2),Cix(p,1):Cix(p,2)) = U*V';     
    end

    % Performance Indexes calculation 
    PI_diag(iter)   = DIAG_crit_proc(); 
    % convergence is attained when off-diag does not change significantly
    CRIT(1,iter) = PI_diag(iter);
    CRIT(2,iter) = abs(PI_diag(iter)-PI_diag(iter-1));
    CRIT(3,iter) = sign(PI_diag(iter)-PI_diag(iter-1));
end
% *************************************************************************

% Convergence attained?...
fprintf('\n');
if iter==max_it
    disp(['JBSS (Li) algorithm has not converged! (eps:' num2str(eps) ')'])
else
    disp(['JBSS (Li) algorithm has converged. (eps:' num2str(eps) ')'])
end


% *** OffDiag/Diag PI calculation ***
function crit = DIAG_crit_proc()
    crit = 0;
    for kk = 1:K
        P(:,:,kk) = B*C(:,:,kk)*B';
        norm_all = trace(P(:,:,kk)'*P(:,:,kk)); % total Frobenius norm
        norm_diag = 0;
        for jj = 1 : I
            for ii = 1 : jj-1 % norm of diags only on the upper triangular part
                norm_diag = norm_diag + 2*sum(diag(P(Cix(jj,1):Cix(jj,2),Cix(ii,1):Cix(ii,2),kk)).^2);
            end
        end
        norm_diag = norm_diag + sum(diag(P(:,:,kk)).^2); % norm of main diag
        crit = crit + (norm_all-norm_diag); % /norm_diag;
    end
    crit = crit / (K*N*(N-1)) ; 
end

end


function x = assert_col(x)
    if size(x,2) > 1    x = x';     end        
end