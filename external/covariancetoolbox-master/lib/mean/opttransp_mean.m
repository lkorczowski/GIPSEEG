% A = riemann_mean(B,epsilon,tol)
%
% Calcul du barycentre des matrice de covariances.
% A : baricentre des K matrices NxN
%
% B : Matrice NxNxK
% epsilon : Pas de la descente de gradient
% tol : arret de la descente si le crit�re < tol


function [A critere niter] = opttransp_mean(B,args)
Imat = size(B,3);
N_itermax = 100;
if (nargin<2)||(isempty(args))
    tol = 10^-3;
    %A = mean(B,3);
    A = eye(size(B,1));
else
    tol = args{1};
    A = args{2};
end

niter = 0;
% fc = 0;
K = A^(0.5);
N=size(A,1);
while (niter<N_itermax)
    niter = niter+1;
   
    Ktmp = 0;
    for i=1:Imat
       Ktmp = Ktmp + gp_powm(K*B(:,:,i)*K,0.5); 
    end
    Ktmp = gp_powm((Ktmp),(0.5));
    
%     fcn = norm(Ktmp-K,'fro');
%         conv(niter) = abs((fcn-fc)/fc);;
     conv(niter) = distance(Ktmp/N,K/N,'fro');%/norm(Pnew,'fro');

%         conv(niter)=distance_ws(Ktmp,K);
    if niter>5 & (conv(niter)-conv(niter-1)>tol) % break if the improvement is below the tolerance
        warning(['solution degenated at iter ' num2str(niter) ])
       break; 
    end
    K = Ktmp;
    % improvement
    if conv(niter)<tol % break if the improvement is below the tolerance
       break; 
    end

%     fc = fcn;
end

A = ((1/Imat)*K)^2;

if niter==N_itermax
    disp('Warning : Nombre d''itérations maximum atteint');
end

critere = conv;
