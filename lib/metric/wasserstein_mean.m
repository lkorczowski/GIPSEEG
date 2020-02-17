function [A, critere, niter]=wasserstein_mean(B,args)
% A = wasserstein_mean(B,epsilon,tol)
%
% Compute geometric mean according to wasserstein metric of non-singular symetric
% matrices
%
% INPUTS
% B : Matrices [NxNxK]
% tol : convergence criterion (defaut : 1e-3)
% A : geometric mean init
% w : weights [K x 1]

% A : baricenter of B [NxN]
N_itermax = 100;
if (nargin<2)||(isempty(args))
    w=ones(size(B,3),1)/size(B,3);
    tol = 10^-9;
    tmp=zeros(size(B));
    for indB=1:size(B,3)
    tmp(:,:,indB) = w(indB)*gp_powm(B(:,:,indB),0.5);
    end
    A=sum(tmp,3)^2;
%     A=
else
    tol = args{1};
    A = args{2};
    w=args{3};
end

niter = 0;
fc = 0;
while (niter<N_itermax)
    niter = niter+1;
    Asqrt=gp_powm(A,0.5);
    % Tangent space mapping
    tmp=zeros(size(B));
    for indB=1:size(B,3)
        tmp(:,:,indB) = w(indB)*gp_powm(Asqrt*B(:,:,indB)*Asqrt,0.5);
    end
    
    Anew=sum(tmp,3);
%     compt(niter)=distance_ws(Anew,A);
     compt(niter) = norm(Anew-A,'fro');%/norm(Anew,'fro');
    A=Anew;
    
    if compt(niter)<tol % break if the improvement is below the tolerance
        break;
    end    
end

% figure
% plot(conv)
if niter==N_itermax
    warning('Maximum iterations reached, wasserstein_mean not converging');
end

critere = compt;
