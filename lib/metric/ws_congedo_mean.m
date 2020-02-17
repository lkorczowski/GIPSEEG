function [A, critere, niter]=ws_congedo_mean(B,args)
% A = wasserstein_mean(B,epsilon,tol)
%
% Compute geometric mean according to wasserstein metric of non-singular symetric
% matrices
%
% ws_congedo_mean(B,tol*,A*,w*,alpha*)
% *optional
%
% INPUTS
% B : Matrices [NxNxK]
% tol : convergence criterion (defaut : 1e-3)
% A : geometric mean init
% w : weights [K x 1]
% alpha : power law coefficient (default -1/2)
%
% OUTPUTS
% A : baricenter of B [NxN]
% critere : converge criteria for each iteration
% niter : number of iterations required
N_itermax = 100;
if (nargin<2)||(isempty(args))
    w=ones(size(B,3),1)/size(B,3);
    tol = 10^-3;
    
    % init
    tmp=zeros(size(B));
    for indB=1:size(B,3)
        tmp(:,:,indB) = w(indB)*gp_powm(B(:,:,indB),1/2);
    end
    A=sum(tmp,3);
    alpha=-1/2;
else
    tol = args{1};
    A = args{2};
    w=args{3};
    alpha=args{4};
end

niter = 0;
fc = 0;
P=pinv(A);
N=size(P,1);
while (niter<N_itermax)
    
    niter = niter+1;
    tmp=zeros(size(A));
    for indB=1:size(B,3)
        tmp =tmp+ w(indB)*gp_powm(P*B(:,:,indB)*P',1/2);
    end
    
    Pnew=gp_powm(tmp,alpha)*P;
    
    %     compt(niter)=distance_ws(Pnew,P);
    conv(niter) = norm((Pnew-P)/N,'fro');%/norm(Pnew,'fro');
    if niter>10 & (conv(niter)-conv(niter-1)>tol) % break if the improvement is below the tolerance
        warning(['solution degenated at iter ' num2str(niter) ])
        break;
    end
    P=Pnew;
    alpha=-((-alpha)^0.75);
    if conv(niter)<tol % break if the improvement is below the tolerance
        break;
    end
    
end

A=pinv(P'*P);

% figure
% plot(conv)
if niter==N_itermax
    warning('Maximum iterations reached, wasserstein_mean not converging');
end

critere = conv;
