function B=gp_powm(A,alpha)
% B=gp_powm(A,alpha)
% Fast matrix power using svd
%
% *** History: 04-March-2016
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2016
%
% see also: powm

if abs(alpha)<1
    [U,S,V]=svd(A);
    B=U*diag((diag(S)).^(alpha))*V';
else
    B=A^(alpha);
end

end
