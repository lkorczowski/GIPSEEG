function [x,err] = genearma(a, b, nm, varinn)

% Syntaxe	  x = genearma(a, b, nm, varinn)
% Generate stationary Gaussian ARMA signals defined by
%    x(i) = a(1)*x(i-1) + ... + a(na)*x(i-na) + e(t) +
%           b(1)e(t-1) + .... +  b(nb)*e(t-nb)
% where e(1:n) is a Gaussian white noise of variance varinn (default to 1)
% The output x is a column vector of length nm if nm is scalar,
% otherwise it is a nm(1) by nm(2) matrix containing nm(2) signals
% of length nm(1). The inputs a and b must be column vectors

% Note: not working on octave (except when nm is a scalar) because the
% filter function there does not accept matrix argument
err=0;
x=[];
if nargin < 3
  uiwait(errordlg('genearma requires at least 3 arguments'));
  err=1;
  return;
end

[na,m] = size(a);
if m > 1
  uiwait(errordlg('AR coefficients must be a column vector'));
  err=1;
  return;
end
[nb,m] = size(b);
if m > 1
  uiwait(errordlg('MA coefficients must be a vector'));
  err=1;
  return;
end

if length(nm) > 1
  n = nm(1); m = nm(2);
else
  n = nm; m = 1;
end
if n < na+nb
  uiwait(errordlg('Sorry: length of the series too short'));
  err=1;
  return;
end

if nargin < 4; v = 1; else v = varinn; end

if na > 0
  for k = na-1:-1:1			% Backward Levinson Durbin algorithm
    r = a(k+1);
    a(1:k) = (a(1:k) + r*a(k:-1:1))/(1 - r^2);
  end
  if abs(a) < 1				% good, |correlation| < 1
    v = v/prod(1 - a.^2);		% variance of the y process
    y = randn(na+nb,m);
    y(1,:) = sqrt(v)*y(1,:);		% 1st y sample
    v = v*(1 - a(1)^2);
    for k = 1:na-1
      y(k+1,:) = a(1:k)'*y(k:-1:1,:) + sqrt(v)*y(k+1,:);
      r = a(k+1);
      a(1:k) = a(1:k) - r*a(k:-1:1);
      v = v*(1 - r.^2);
    end
    for k = na+1:na+nb
      y(k,:) = a(1:na)'*y(k-1:-1:k-na,:) + sqrt(v)*y(k,:);
    end
  else
    uiwait(errordlg('AR polynomial unstable'));
    err=1;
    return;
  end
  x = toeplitz([1; zeros(na-1,1)],[1; -b; zeros(na-1,1)])*y;
  z = hankel(a)*x;
  if nb > 0
    if nb > na; z = [z; zeros(nb-na,m)]; end
    z(1:nb,:) = z(1:nb,:) + hankel(b)* ...
                toeplitz([1; zeros(nb-1,1)],[1; -a; zeros(nb-1,1)])*y;
  end
  x = [x
       filter([1; b], [1; -a]/sqrt(v), randn(n-na,m), z)];
else					% no AR part
  if nb == 0				% simple white noise
    x = sqrt(v)*randn(n,m);
    return
  end
  z = hankel(b)*sqrt(v)*randn(nb,m);
  x = filter([1; b], 1/sqrt(v), randn(n,m), z);
end
