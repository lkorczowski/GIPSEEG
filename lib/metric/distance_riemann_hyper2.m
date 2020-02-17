function a = distance_riemann_hyper2(A,B,NbE)
%distance riemann

size1=1:NbE;
size2=NbE+1:NbE*2;

Cp1=A(size1,size1);
Cp2=A(size2,size2);
C1p1=B(size1,size1);
C1p2=B(size2,size2);
EIG1=eig(C1p1,Cp1);L1=length(EIG1);
EIG2=eig(C1p2,Cp2);L2=length(EIG2);
a = sqrt(sum(log([EIG1(1:floor(1)); EIG2(1:floor(1))]).^2));