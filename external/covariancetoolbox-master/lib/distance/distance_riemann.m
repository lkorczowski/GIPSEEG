function [a EIG] = distance_riemann(A,B,method_dist_eig)
EIG=eig(A,B);
validEIG=find(EIG>0);if (length(validEIG)~=length(EIG)), warning('singular or negative eigenvalue(s) and discarded in riemannian distance'); end
a = sqrt(sum(log(EIG(validEIG)).^2));