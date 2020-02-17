function [B,P,CRIT] = JBSS_li_grad(C, S, eps, max_it, B_init)
% function [B,P,CRIT] = JBSS_li_grad(C, S, eps, max_it, B_init)
%
% ************************************************************************
% This function performs Joint Blind Source Separation (JBSS) using an
% implementation of Li's gradient descent algorithm for Generalized Joint
% Approximate Diagonalization.
%
% /!\ here we use a home-made gradient descend with optimal step (therefore 
% results are not the same as provided by algorithm in Li's paper)
% 
% It produces a mixing matrix B (M components * N channels) 
% such as:   C_out(k) = B.C(k).B' has a "stripe-diagonal" form
%
% Input:
% ------
%   C   :   3D array containing the different matrices, may be given in
%           unfolded version (channels*channels*cov_struct)
%   S   :   the size of the square submatrices D_ii, i.e. S is a vector
%   eps :   criterion threshold
%   max_it: iteration number
%   B_init: (optional) initial demixing matrix
%
% Output:
% -------
%   B   : demixing matrix whose product B*C*B' is stripe-diagonal
%   P   : the stripes diagonalised version of C
%   CRIT: stopping criterion, based on off-diag of P matrices
%               CRIT(1,:) absolute change criterion (convergence)
%               CRIT(2,:) sign (-1 = converge, +1 = diverge)
%
% History:
% --------
% *** 2011-05-16
% Created by J. Chatel-Goldman @ GIPSA Lab, based on work of:
% X.-L. Li, T. Adal?, and M. Anderson, "Joint blind source separation by 
% generalized joint diagonalization of cumulant matrices,"
% Signal Processing, vol. 91, no. 10, pp. 2314–2322, 2011.

S = assert_col(S); % verify if S is a column vector as it should be
% minStep = 1e-16;                % minimum step reached at convergence
I       = length(S);            % number of subjects
N       = sum(S);               % total number of channels 
K       = size(C,3);            % number of covariance structures
cix     = [cumsum(S)-S+1, cumsum(S)]; % index of Cij matrices for any fixed i
iter    = 1;

% Memory allocation
PI_diag     = zeros(1,max_it);  % init stripe-diags Performance Index
CRIT  = zeros(3,max_it);
% CRIT(1,:) --> convergence absolute value
% CRIT(2,:) --> convergence absolute change value
% CRIT(3,:) --> convergence absolute change sign
P       = zeros(N,N,K);

% init Bii to unit values
B       = zeros(N);
if (nargin < 6) || isempty(B_init)
    for i=1:I
        B(cix(i,1):cix(i,2),cix(i,1):cix(i,2)) = eye(S(i)); % initialization to identity
    end
else
	B = B_init;    
end

% (first) Performance Index calculation 
PI_diag(iter)   = DIAG_crit_proc();


% ******************************* MAIN LOOP *******************************
while iter<=1 || ((CRIT(2,iter)>eps && iter<max_it))
    % UI display
    fprintf('%s','.'); if(mod(iter,50) == 0) fprintf('\n'); end
    iter = iter+1; 

    grad    = zeros(N,N);    
    for i= 1:I % p
        
        offB = zeros(S(i));
        offA = zeros(S(i));
        
        % Gradient calculation
        for k = 1:K % l
            for j = 1:I % q
                tempOffA =  B(cix(i,1):cix(i,2),cix(i,1):cix(i,2))' * ...
                        C(cix(i,1):cix(i,2),cix(j,1):cix(j,2),k) * ...
                        B(cix(j,1):cix(j,2),cix(j,1):cix(j,2));
                tempOffA = tempOffA-diag(diag(tempOffA));
                grad(cix(i,1):cix(i,2),cix(i,1):cix(i,2)) = grad(cix(i,1):cix(i,2),cix(i,1):cix(i,2)) + ...
                    tempOffA * B(cix(j,1):cix(j,2),cix(j,1):cix(j,2))' * C(cix(i,1):cix(i,2),cix(j,1):cix(j,2),k)';

                off =   B(cix(i,1):cix(i,2),cix(i,1):cix(i,2))' * ...
                        C(cix(j,1):cix(j,2),cix(i,1):cix(i,2),k)' * ...
                        B(cix(j,1):cix(j,2),cix(j,1):cix(j,2));
                off = off-diag(diag(off));
                grad(cix(i,1):cix(i,2),cix(i,1):cix(i,2)) = grad(cix(i,1):cix(i,2),cix(i,1):cix(i,2)) + ...
                    off * B(cix(j,1):cix(j,2),cix(j,1):cix(j,2))' * C(cix(j,1):cix(j,2),cix(i,1):cix(i,2),k);
                
                % offB is used next in optimal step calculation
                tempOffB =  grad(cix(i,1):cix(i,2),cix(i,1):cix(i,2))* ...
                            C(cix(i,1):cix(i,2),cix(j,1):cix(j,2),k) * ...
                            B(cix(j,1):cix(j,2),cix(j,1):cix(j,2));
                tempOffB = tempOffB-diag(diag(tempOffB));
                offB = offB + tempOffB;
                offA = offA + tempOffA;
            end
            
        end
        
        % Calculation of optimal step 
        stepOpt = trace(offA*offB') / trace(offB*offB');
        NewBi = B(cix(i,1):cix(i,2),cix(i,1):cix(i,2)) - stepOpt*grad(cix(i,1):cix(i,2),cix(i,1):cix(i,2))';
        B(cix(i,1):cix(i,2),cix(i,1):cix(i,2)) = NewBi*(NewBi'*NewBi)^-.5;     % orthogonalization
        
        % Gradient descent
%         step = 10;
%         Bgrad = B;
%         currentCost = wholeCost(B);
%         while(step > minStep) 
%             NewBi = B(cix(i,1):cix(i,2),cix(i,1):cix(i,2)) - step*grad(cix(i,1):cix(i,2),cix(i,1):cix(i,2))';
%             NewBi = NewBi*(NewBi'*NewBi)^-.5;     % orthogonalization
%             Bgrad(cix(i,1):cix(i,2),cix(i,1):cix(i,2)) = NewBi;  % zeros everywhere for Bi
%             if( wholeCost(Bgrad) < currentCost )
%                 B = Bgrad; 
%                 step = 0;
%             else
%                 step = step * .5;
%             end        
%         end
    end
    
    % Performance Index calculation (offdiagonality)
    PI_diag(iter)   = DIAG_crit_proc(); 
    
    % convergence is attained when off-diag does not change significantly
    CRIT(1,iter) = PI_diag(iter);
    CRIT(2,iter) = abs(PI_diag(iter)-PI_diag(iter-1));
    CRIT(3,iter) = sign(PI_diag(iter)-PI_diag(iter-1));
end

%  ********* FINALIZATION *********
% VERY IMPORTANT: B transposition for respect of convention C_out = B*C*B';
B = B';
% Last thing: process P for every k slice
for k = 1:K 
    P(:,:,k) = B*C(:,:,k)*B';
end
% Convergence attained?...
fprintf('\n');
if iter==max_it
    disp(['JBSS (Li''s gradient descent) algorithm has not converged! (eps:' num2str(eps) ')'])
else
    disp(['JBSS (Li''s gradient descent) algorithm has converged. (eps:' num2str(eps) ')'])
end
% *********************************


% *** OffDiag/Diag PI calculation ***
function crit = DIAG_crit_proc()
    crit = 0;
    for kk = 1:K
        P(:,:,kk) = B'*C(:,:,kk)*B; 
        norm_all = trace(P(:,:,kk)'*P(:,:,kk)); % total Frobenius norm
        norm_diag = 0;
        for jj = 1 : I
            for ii = 1 : jj-1 % norm of diags only on the upper triangular part
                norm_diag = norm_diag + 2*sum(diag(P(cix(jj,1):cix(jj,2),cix(ii,1):cix(ii,2),kk)).^2);
            end
        end
        norm_diag = norm_diag + sum(diag(P(:,:,kk)).^2); % norm of main diag
        crit = crit + (norm_all-norm_diag); % /norm_diag;       % WHY /norm_diag ?
    end
    crit = crit / (K*N*(N-1)) ; 
end

% % *** Whole cost calculation *** 
% function cost = wholeCost(X)
%     cost = 0;
%     for kk = 1:K
%         P(:,:,kk) = X'*C(:,:,kk)*X;
%         norm_all = trace(P(:,:,kk)'*P(:,:,kk)); % total Frobenius norm
%         norm_diag = 0;
%         for jj = 1 : I
%             for ii = 1 : jj-1 % norm of diags only on the upper triangular part
%                 norm_diag = norm_diag + 2*sum(diag(P(cix(jj,1):cix(jj,2),cix(ii,1):cix(ii,2),kk)).^2);
%             end
%         end
%         norm_diag = norm_diag + sum(diag(P(:,:,kk)).^2); % norm of main diag
%         cost = cost + (norm_all-norm_diag);
%     end 
% end

end


function x = assert_col(x)
    if size(x,2) > 1    x = x';     end        
end




