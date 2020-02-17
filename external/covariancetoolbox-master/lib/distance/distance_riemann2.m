function a = distance_riemann2(A,B)
EIG=eig(A,B);EIG(end);
L=length(EIG);
a = sqrt(sum(log(EIG(end)).^2));