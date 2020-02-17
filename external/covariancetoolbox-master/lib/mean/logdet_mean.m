% from Means of Hermitian positive-definite matrices based on the log-determinant α-divergence function
% Zeineb Chebbi , Maher Moakher.


function [A conv niter] = logdet_mean(B,tol,MaxIter,A)

K = size(B,3); % Nombre de matrices

if (nargin<4)||(isempty(A))
    
    A = mean(B,3);
end
if (nargin<3)||(isempty(MaxIter))
    MaxIter = 100;
end
if (nargin<2)||(isempty(tol))
    tol = 10^-3;
end


conv = tol + 1;
niter=1;
N=size(A,1);
while conv(niter)>tol && niter<MaxIter
    niter=niter+1;
    fc = zeros(size(B,1));
    
    for i=1:K
        fc = fc + inv(0.5*B(:,:,i) + 0.5*A);
    end
    
    Anew = inv(fc/K);
%     conv(compt) = distance_ld(Anew,A);
     conv(niter) = norm((Anew-A)/N,'fro');%/norm(Pnew,'fro');
    if niter>15 & (conv(niter)-conv(niter-1)>tol) % break if the improvement is below the tolerance
        warning(['solution degenated at iter ' num2str(niter) ])
       break; 
    end
    A = Anew;
end
if niter>=MaxIter
    disp('Maximum iteration reached, Logdet_mean cannot converge')
else
    %disp(['Mean cov computed in ' int2str(compt) ' iterations'])
end


function B = mypseudoinv(A)
[U S] = svd(A);
k=1;
K = size(A,1);
B = zeros(K,K);
while (k<=K) && (S(k,k)>10^-15)
    B = B+U(:,k)*(1/S(k,k))*U(:,k)';
    k=k+1;
end



