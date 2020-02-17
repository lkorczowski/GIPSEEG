% A = riemann_mean(B,epsilon,tol)
%
% Calcul du barycentre des matrice de covariances.
% A : baricentre des K matrices NxN
%
% B : Matrice NxNxK
% epsilon : Pas de la descente de gradient
% tol : arret de la descente si le crit�re < tol


function [A, critere, niter] = riemann_mean(B,tol,N_itermax,A,iterDiv)
if (nargin<5)||(isempty(iterDiv))
    
    iterDiv=10;
end
if (nargin<4)||(isempty(A))
    
    A = mean(B,3);
end
if (nargin<3)||(isempty(N_itermax))
    N_itermax = 100;
end
if (nargin<2)||(isempty(tol))
    tol = 10^-3;
end



niter = 0;
n=size(A,1);
TA =zeros(n*(n+1)/2,1);
while (niter<N_itermax)
    niter = niter+1;
    % Tangent space mapping
    T = Tangent_space(B,A);
    
    % sum of the squared distance (frobenius norm)
    % improvement
    
    % arithmetic mean in tangent space
    
    TAn = 1*mean(T,2);
    N=size(TAn,1);
    conv(niter) = norm((TAn-TA)/sqrt(N),'fro');%/norm(Pnew,'fro') sum of squared distance
    TA=TAn;
    if niter>iterDiv & (conv(niter)-conv(niter-1)>tol) % break if the improvement is below the tolerance
        warning(['solution degenated at iter ' num2str(niter) ])
        break;
    end
    % back to the manifold
    A = UnTangent_space(TA,A);
    if conv(niter)<tol % break if the improvement is below the tolerance
        break;
    end
    
end

% figure
% plot(conv)
if niter==N_itermax
    disp('WARNING: Maximum iterations reached, riemann_mean not converging');
end

critere = conv;
