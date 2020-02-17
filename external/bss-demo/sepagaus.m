function [sep, seps] = sepagaus(data, nfreq, bloclen, sep, eps, maxiter)

% Syntaxe  [sep, seps] = sepagaus(data, nfreq, bloclen, sep, eps, maxiter)
%
% Separation of source through the Gaussian mutual information criterion.
% * data is a n by K matrix, K being the number of sources
% * nfreq is the number of positive frequency channels to use
% * bloclen is the "stationarity length" of the data, it must be
%   greater than nfreq*K and preferably a divisor of n, otherwise the
%   last incomplete block of data is discard.
% * sep is the initial separation matrix.
% * seps, if present, is a (iter+1)*K by K matrix containing the sequence
%   of iterated separation matrices (iter = the number of iterations)
% * eps defines the stoping criterion, which is the square norm
%   (with respect to a certain metric) of the relative "gradient" is
%   less it. This square norm also equals approximatively the decrease
%   of the criterion for this step
% * maxiter is maximum number of iteration of the algorithm
%
% nfreq defaults to 2, bloclen to fix(n/2), sep to the identity matrix, 
% eps to K*(K-1)*1e-8 and maxiter to 20 

[n,K] = size(data);
if nargin < 6; maxiter = 20; end
if nargin < 5; eps = K*(K-1)*1e-8; end
if nargin < 3
    bloclen = fix(n/2);
elseif n < bloclen
    error('data length should not be less than block length')
end
if nargin < 2; nfreq = 2; end
if bloclen <= nfreq*K
    error(['block length too low, '...
            'too many frequency chanels or data length too short'])
end

% Compute the matrices to be diagonalized

filt = [0.5:nfreq-.5, nfreq-.5:-1:.5];		% length 2*freq
scale = sqrt((4*nfreq^2 - 1)*nfreq/6);

spec = [];
if nfreq > 1
  for nf = 1:nfreq			% loop through frequencies
        
    % Demodulation
        
    %    tmp = exp((1:n)*(i*(nf-.5)*pi/nfreq))';
    z = exp(i*(nf-.5)*pi/nfreq);
    %    tmp = cumprod(z(ones(n,1)));
    %    x = tmp(:,ones(1,K)).*data;
    tmp = cumprod(z(ones(4*nfreq,1)));
    x =      tmp(rem((0:n-1),4*nfreq)+1,ones(1,K)).*data;
    x = filter(filt,scale,[x; zeros(nfreq,K)]);	% add nfreq elements
    x(1:nfreq,:) = [];				% discard 1st nfreq elements
        
    for j = bloclen:bloclen:n			% loop through blocks
      xr = real(x(j-bloclen+1:j,:));
      xi = imag(x(j-bloclen+1:j,:));
      spec = [spec, xr'*xr + xi'*xi];
    end
  end
else						% use time diversity only
  for j = bloclen:bloclen:n			% loop through blocks
    xr = data(j-bloclen+1:j,:);
    spec = [spec, xr'*xr];
  end
end
if nargin < 4
  sep = eye(K);
elseif any(any(sep - eye(K)))			% do inital transformation
  spec = sep*spec;
  for j = K:K:length(spec)
    spec(:,j-K+1:j) = spec(:,j-K+1:j)*sep';
  end
end

if nargout > 1, seps = sep; end
for it = 1:maxiter
  [spec, sep, decr] = jadiag1(spec, sep);
  if nargout > 1, seps = [seps; sep]; end
  if decr < eps; break; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c, a, decr] = jadiag1(c, a)
% ***************** Internal function ***********************
% same as the stand-alone version, but no output crit, logdet,
% _no defaults, no check_ and assume REAL matrices

% Performs approximate joint diagonalization of several matrices.
% The matrices to be diagonalised are given in concatenated form by a
% m x n matrix c, n being a multiple of m. They are tranformed to
% more diagonal with _one_ iteration sweep. The function returns the
% transformed matrices in c. The transformation is also applied
% (through pre-multiplication only) to the matrix a (default to the
% identity matrix if not provided) and the result is returned. Further
% decr provides an estimate decrease of the critertion for this step
% which is also the squared norm off the gradient with respect to a
% certain metric, to be use a a stopping criterion

[m, n] = size(c);
nmat = fix(n/m);

one = 1 + 10e-12;			% considered as equal to 1

decr = 0;
for i = 2:m
  for j=1:i-1
    c1 = c(i,i:m:n);
    c2 = c(j,j:m:n);
    g12 = mean(c(i,j:m:n)./c1);	% this is g_{ij}
    g21 = mean(c(i,j:m:n)./c2);	% this is the conjugate of g_{ji}
    omega21 = mean(c1./c2);
    omega12 = mean(c2./c1);
    omega = sqrt(omega12*omega21);
    tmp = sqrt(omega21/omega12);
    tmp1 = (tmp*g12 + g21)/(omega + 1);
    omega = max(omega, one);
    tmp2 = (tmp*g12 - g21)/(omega - 1);
    h12 = tmp1 + tmp2;			% this is twice h_{ij}
    %      h21 = conj((tmp1 - tmp2)/tmp);	% this is twice h_{ji}
    %      decr = decr + nmat*(g12*conj(h12) + g21*h21)/2;
    %      tmp = 1 + 0.5i*imag(h12*h21);% = 1 + (h12*h21 - conj(h12*h21))/4
    %      T = eye(2) - [0 h12; h21 0]/(tmp + sqrt(tmp^2 - h12*h21));
    %
    % this section replaces the above since the matrices are assumed real
    %
    h21 = (tmp1 - tmp2)/tmp;                  % this is twice h_{ji}
    decr = decr + nmat*(g12*h12 + g21*h21)/2; % above since c is real
    T = eye(2) - [0 h12; h21 0]/(1 + sqrt(1 - h12*h21));
    
    c([i j],:) = T*c([i j],:);
    for k=0:m:n-m
      c(:,[i+k j+k]) = c(:,[i+k j+k])*T';
    end
    a([i j],:) = T*a([i j],:);
  end
end
