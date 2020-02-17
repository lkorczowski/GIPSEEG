% A = riemann_mean(B,epsilon,tol)
%
% Calcul du barycentre des matrice de covariances.
% A : baricentre des K matrices NxN
%
% B : Matrice NxNxK
% epsilon : Pas de la descente de gradient
% tol : arret de la descente si le crit�re < tol


function [A, critere, niter] = riemann_mean_stepwise(B,args)

N_itermax = 100;
if (nargin<2)||(isempty(args))
    tol = 10^-9;
    A=logdet_mean(B);
else
    tol = args{1};
    A = args{2};
end

niter = 0;
fc = 0;
sigma=1;
conv=realmax;
while (niter<N_itermax)
    niter = niter+1;
    % Tangent space mapping
    % update M=1/K*sum(ln(M^(-1/2)*Ck*M^(-1/2))
    T = sigma*Tangent_space(B,A); 
    % sum of the squared distance (frobenius norm) 
    fcn = sqrt(sum(sum(T.^2)));
    % improvement
    if length(conv)>1
        fc=conv(niter-1);
    else
        fc=conv;
    end
    %if fcn<fc % if there is a progression
        conv(niter)=fcn;
        sigma=0.95*sigma;
        
        % arithmetic mean in tangent space
        TA = sigma*mean(T,2);
        % back to the manifold
        A = UnTangent_space(TA,A); %
    %else % if no progression
    %    sigma=0.5*sigma;
    %    conv(niter)=fc;
    %end
    
        if conv(niter)<tol||sigma<tol % break if the improvement is below the tolerance
           break; 
        end
    
end

figure
plot(conv)
if niter==N_itermax
    disp('WARNING: Maximum iterations reached, riemann_mean_stepwise not converging');
end

critere = conv(end);
