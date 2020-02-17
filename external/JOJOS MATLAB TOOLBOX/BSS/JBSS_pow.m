function [B,C_out,CRIT] = JBSS_pow(C, S, eps, max_it, ortho, B_init)
% function [B,C_out,CRIT] = JBSS_pow(C, S, eps, max_it, ortho, B_init)
%
% ************************************************************************
% This function performs Joint Blind Source Separation (JBSS) using an
% implementation of OJoB (orthogonal) and NOJoB (non-orthogonal) algorithms
% 
% It produces a mixing matrix B (M components * N channels) 
% such as:   C_out(k) = B.C(k).B' has a "stripe-diagonal" form
%
% Input:
% ------
%   C   :   3D array containing the different matrices, may be given in
%           unfolded version channels*channels*(freqs/lags/condition/...)
%   S   :   the size of the square submatrices D_ii, i.e. S is a vector
%   eps :   criterion threshold
%   max_it: iteration number
%   ortho : (optional) specifies which version of power iteration is used
%           (orthogonal=1 or not=0, default is orthogonal)
%   B_init: (optional) initialization of the demixing matrix
%
% Output:
% -------
%   B   : demixing matrix whose product B*C*B' is stripe-diagonal
%   C_out   : the stripes diagonalised version of C
%   CRIT: stopping criterion, based on off-diag of C_out matrices
%               CRIT(1,:) absolute change criterion (convergence)
%               CRIT(2,:) sign (-1 = converge, +1 = diverge)
%
% History:
% --------
% *** 2011-10-26
% Created by J. Chatel-Goldman @ GIPSA Lab, 
% based on work of M. Congedo and R. Phlypo


% assertion on inputs
if sum(S)~=size(C,1) error('JBBS_pow2: Incoherent submatrices size description'); end
S = assert_col(S); % verify if S is a column vector as it should be

% initialization
iter = 1;                   % iteration number
I   = length(S);            % number of subjects
N   = sum(S);               % total number of channels    
K   = size(C,3);            % number of covariance layers
cix = [cumsum(S)-S+1, cumsum(S)]; % index of Cij matrices 
C_out       = zeros(N,N,K);
PI_diag = zeros(1,max_it);  % init stripe-diags Performance Index
CRIT  = zeros(3,max_it);
% CRIT(1,:) --> convergence absolute value
% CRIT(2,:) --> convergence absolute change value
% CRIT(3,:) --> convergence absolute change sign
if nargin < 5   ortho = 1;  end     % default is orthogonal version

% init Bii to unit values
B       = zeros(N);
if (nargin < 6) || isempty(B_init)
    for i=1:I
        B(cix(i,1):cix(i,2),cix(i,1):cix(i,2)) = eye(S(i)); % initialization to identity
    end
else
	B = B_init';
end

% Performance Indexe calculation 
PI_diag(1)   = DIAG_crit_proc(); 


% ******************************* MAIN LOOP *******************************
while iter<=1 || ((CRIT(2,iter)>eps && iter<max_it))
    
    % UI display
%    disp(['Processing AJSVD with power iterations, iteration #' int2str(iter)]);
    fprintf('%s','.'); if(mod(iter,50) == 0) fprintf('\n'); end
    iter = iter+1; 
        
    % for each subject
    for i = 1:I
        Bii_new = B(cix(i,1):cix(i,2),cix(i,1):cix(i,2)); % select Bii matrix
        if(ortho == 0)      
            Mtot    = zeros(S(i));              % sum of all Mi(n)
            Mi_all  = zeros(S(i),S(i)*S(i));    % all the Mi(n), stored for consequent normalization
        end   
        
        % for each column of Bjj
        for n = 1:S(i)  
            M = zeros(S(i));    % for each column, initialize Mi(n) at zero
            % for each slice k
            for k = 1:K 
                for j = 1:I     
                    if n<=S(j)  % columns n outside Bij matrix are discarded
                        Ctemp = C(cix(i,1):cix(i,2), cix(j,1):cix(j,2), k); % select C(k,ij)
                        v = Ctemp * B(cix(j,1):cix(j,2), cix(j,1)+n-1);     % v = C(k,ij)*b(jn)
                        M = M + v*v';        % Mi(n) = Mi(n) + ( C(k,ij)*b(jn) )^2
                    end
                end
            end
           
            if(ortho == 1)  % Orthogonal case  
                Bii_new(:,n) = M * Bii_new(:,n); % power iteration on Bii,n  
            
            else            % Non-orthogonal case
                Mtot = Mtot + M;                        % process sum(Mi(n))
                Mi_all(:,(((n-1)*S(i))+1):n*S(i)) = M;  % store Mi(n)   
            end
        end     
        
        if(ortho == 1)      
            % orthogonalize (also normalize) new Bii and replace old Bii
            B(cix(i,1):cix(i,2),cix(i,1):cix(i,2)) = orthogonalize(Bii_new);
        else
            % Linear system solving using Cholesky decomposition 
            %Bii_new = linsolve(Mtot,Bii_new); % Mtot * Bii_new = Bii_new (-->solution Ronald a priori bad)
            L = chol(Mtot)';
            % Normalization: each column of Bjj with its corresponding matrix Mi(n)
            % B(cix(i,1):cix(i,2),cix(i,1):cix(i,2)) = Bii_new*(Bii_new'*Bii_new)^(-1/2);    % <-- A REFAIRE!!!!!
            for n = 1:S(i)  
                x = linsolve( L, Mi_all(:,(((n-1)*S(i))+1):n*S(i)) * B(cix(i,1):cix(i,2), cix(i,1)+n-1) );  % solve L*x  = Mi(n)*Bi(n)
                y = linsolve( L', x);                                                                       % solve L'*x = y
                Bii_new(:,n) = y * ((y'*Mi_all(:,(((n-1)*S(i))+1):n*S(i))*y)^(-1/2));                       % bi(n)= y / (y'*Mi(n)*y)^(1/2)
            end
            B(cix(i,1):cix(i,2),cix(i,1):cix(i,2)) = Bii_new;
        end
    end  
    
    % Performance Index calculation (offdiagonality)
    PI_diag(iter)   = DIAG_crit_proc(); 
              
    % convergence is attained when off-diag does not change significantly
    CRIT(1,iter) = PI_diag(iter);
    CRIT(2,iter) = abs(PI_diag(iter)-PI_diag(iter-1));
    CRIT(3,iter) = sign(PI_diag(iter)-PI_diag(iter-1));
end
% *************************************************************************


%  ********* FINALIZATION *********
% B transposition for respect of convention C_out = B*C*B';
B = B';

% Last thing: process C_out for every k slice
for k = 1:K 
    C_out(:,:,k) = B*C(:,:,k)*B';
end

% Convergence attained?...
fprintf('\n');
if ortho == 1
    algo_str = 'OJoB';
elseif ortho == 0
    algo_str = 'NOJoB'; 
end
if iter==max_it
    disp(['JBSS (' algo_str ' algorithm) has not converged! (eps:' num2str(eps) ')'])
else
    disp(['JBSS (' algo_str ' algorithm) has converged. (eps:' num2str(eps) ')'])
end
% **********************************



% *** orthogonalization *** 
function B = orthogonalize(B)
    [U,Dsvd,V] = svd(B);
    B = U*V';
end
    
% *** OffDiag/Diag PI calculation ***       
function crit = DIAG_crit_proc()
    crit = 0;
    for k = 1:K
        C_out(:,:,k) = B'*C(:,:,k)*B;  
        norm_all = trace(C_out(:,:,k)'*C_out(:,:,k)); % total Frobenius norm
        norm_diag = 0;
        for j = 1 : I
            for i = 1 : j-1 % norm of diags only on the upper triangular part
                norm_diag = norm_diag + 2*sum(diag(C_out(cix(j,1):cix(j,2),cix(i,1):cix(i,2),k)).^2);
            end
        end
        norm_diag = norm_diag + sum(diag(C_out(:,:,k)).^2); % norm of main diag
        crit = crit + (norm_all-norm_diag); % /norm_diag;       % WHY /norm_diag ?
    end
    crit = crit / (K*N*(N-1)) ; 
end

end

function x = assert_col(x)
    if size(x,2) > 1    x = x';     end        
end
