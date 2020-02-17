function [c, crit, a, logdet, decr] = jadiag(c, a, logdet)
% syntaxe       [c, crit, a, logdet, decr] = jadiag(c, a, logdet)

% Performs approximate joint diagonalization of several matrices.
% The matrices to be diagonalised are given in concatenated form by a
% m x n matrix c, n being a multiple of m. They are tranformed to
% more diagonal with an iteration sweep. The function returns the
% transformed matrices in c. The transformation is also applied (through
% pre-multiplication only) to the matrix a (default to the identity matrix
% if not provided) and the result is returned. Further the variable logdet
% is added by twice the number of matrices to be diagonalized times the
% logarithm of the determinant of the transformation (logdet default to 0
% if not given) and is returned. Finally crit contains the logarithm of
% the product of the diagonal elements of all the transformed matrices
% minus the new value of logdet.

[m, n] = size(c);
nmat = fix(n/m);
if (n > nmat*m)
  error('argument must be the concatenation of square matrices')
end

if (nargin < 3)
  logdet = 1;
  if (nargin < 2); a = eye(m); endsd
end
one = 1 + 10e-12;
det = 1;
decr = 0;
for i = 2:m
  for j=1:i-1
     c1 = c(i,i:m:n);
     
    c2 = c(j,j:m:n);
    g12 = mean(c(i,j:m:n)./c1);
    g21 = mean(c(i,j:m:n)./c2);
    omega21 = mean(c1./c2);
    omega12 = mean(c2./c1);
    omega = sqrt(omega12*omega21);
    tmp = sqrt(omega21/omega12);
    tmp1 = (tmp*g12 + g21)/(omega + 1);
    omega = max(omega, one);
    tmp2 = (tmp*g12 - g21)/(omega - 1);
    h12 = (tmp1 + tmp2);
    h21 = (tmp1 - tmp2)/tmp;
    decr = decr + m*(g12*h12 + g21*h21)/2;
    tmp = 1 + (h12*h21 - conj(h12*h21))/4;   
    T = eye(2) - [0 h12; h21 0]/(tmp + sqrt(tmp^2 - h12*h21));    
    c([i j],:) = T*c([i j],:); 
    for k=0:m:n-m
      c(:,[i+k j+k]) = c(:,[i+k j+k])*T';
    end
    a([i j],:) = T*a([i j],:);
    det = det*abs(1 - T(1,2)*T(2,1));
  end
end
logdet = logdet + 2*nmat*log(det);
crit = 1;
for k=1:m:n
   crit = crit*prod(diag(c(:,[k:k-1+m])));
   end
crit = log(crit) - logdet;
