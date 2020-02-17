function [a EIG] = distance_riemann_hyper(A,B,NbE,method_dist_eig)
%distance riemann
if nargin==4 % intra stats only
    size1=1:NbE;
    size2=NbE+1:NbE*2;
    
    %first block
    Cp1=A(size1,size1);
    Cp2=A(size2,size2);
    
    %second block
    C1p1=B(size1,size1);
    C1p2=B(size2,size2);
    
    %compute each eigenvalues
    EIG1=eigsort(C1p1,Cp1);
    EIG2=eigsort(C1p2,Cp2);
    
    if method_dist_eig==0 %all eigenvalues
        EIG=[EIG1;EIG2];
        a = sqrt(sum(log(EIG).^2));
    else
        EIG=[EIG1(method_dist_eig); EIG2(method_dist_eig)];
        a = sqrt(sum(log(EIG).^2));
    end
    
elseif nargin==3 % intra + inter stats
    method_dist_eig=NbE;
    EIG=eigsort(A,B);
    if method_dist_eig==0%all eigevalues
        a = sqrt(sum(log(EIG).^2));
    else
        EIG=EIG(method_dist_eig);
        a = sqrt(sum(log(EIG).^2));
        
    end
    
else % normal riemannian distance
    EIG=eigsort(A,B);
    a = sqrt(sum(log(EIG).^2));
    
end