function [y, z, B, Bs] = icainfpnl(x, np, maxiter, bin, mu_B, mu_z);

% Parameter check and defaults

[n, K] = size(x);

%if nargin < 2, np = 8; end
%if nargin < 3, maxiter = 40; end;
%if nargin < 4, bin = 1; end
%%if nargin < 5, mu_B = .5; end;
%if nargin < 5, mu_B = .6; end;
%%if nargin < 6, mu_z = .2; end;
%if nargin < 6, mu_z = .3; end;

filt=1;
bin = bin*(352*sqrt(pi)/(15*n))^.2;	% 352 = 2^5*11
% more generally filt and bin can be set as below, for some integer nb
% filt = cumsum(1:nb-1);
% filt = [filt, (nb*(nb-1)/2+(1:nb).*(nb:-1:1)), filt(nb-1:-1:1)]/nb^3;
% bin = bin*(352*sqrt(pi)/(15*n))^.2/nb;

ALF = 1.0e-4;
TOLX2 = 10e-4/n;
EULER = .5772156649;            % the Euler constant

% the spline basis of 1st degree

basis = zeros(n,np);

i1 = fix((n+1-n*cos(pi/(np+1)))/2);
u = acos(1+1/n-(2/n)*(1:i1)')*((np+1)/pi);		% u <= 1

basis(1:i1,1) = u;
for k=2:np
  i0 = i1;
  i1 = fix((n+1 - n*cos(k*pi/(np+1)))/2);
  u = acos(1+1/n-(2/n)*(i0+1:i1)')*((np+1)/pi) - (k-1);	% u <= k - (k-1)
  basis(i0+1:i1,k-1) = (1 - u);
  basis(i0+1:i1,k) = u;
end
u = acos(1+1/n-(2/n)*(i1+1:n)')*((np+1)/pi) - np;	% u <= k - (np-1)
basis(i1+1:n,np) = (1-u);

% Initialisation

[rx, p] = sort(x);			% p temporally stores the indexes
for k=1:K
  rx(p(:,k),k) = (1:n)';
end

p = exp(-erfinv((1-n:2:n-1)'/n).^2)/sqrt(2*pi);	% normal distrib.
%p = ones(n,1);
p0 = (basis'*basis)\(basis'*p);
p = basis*p0;

%compute the quantile z from p, the auxilary array lambda and the entropy
z = diff(p);				% save diff(p) temporally
lamda = z./p(1:n-1);
non0 = abs(lamda) > 1e-5;
z(non0) = log(1+lamda(non0))./(z(non0)*n);
z(~non0) = (1/n)./p(~non0);
entz = sum(log(z))/(n-1);
lamda(non0) = 1./log(lamda(non0)+1) - 1./lamda(non0);
lamda(~non0) = .5;
z = [0; cumsum(z)];
z = z - sum(z)/n;

%make K identical copies
p0 = p0(:,ones(1,K));
p = p(:,ones(1,K));
z = z(:,ones(1,K));
lamda = lamda(:,ones(1,K));
for k=1:K
  y(:,k) = z(rx(:,k),k);
end

%if nargin < 7
%  B = eye(K);			% start with B = identity matrix
%else
%  y = y*B';			% start with B provided by the user
%end
B = eye(K);

% Inital estimation of the covariance matrix of Y, of its
% entropies and score functions and of the criterion

cv = y'*y/n;
for k = 1:K
  [enty(k), psiy(:,k)] = kentropy(y(:,k),bin*sqrt(cv(k,k)),filt);
end

crit0 = sum(enty) - K*(entz+log(n-1)+EULER) - log(abs(det(B)));
%fprintf(' 0 Criterion %.5g\n', crit0);

% Iteration loop

if nargout > 3; Bs = B; end
for iter = 1:maxiter

  % Linear part of the gradient

  G = psiy'*y/n;

  m = sum(psiy)/n;
  psiy = psiy - m(ones(n,1),:);			% centrer psiy
  m = ((1 - diag(G))./diag(cv))';		% for the "correction" of psiy
  psiy = psiy + y.*m(ones(n,1),:);		% psiy corrected

  G = G + m(ones(K,1),:)'.*cv - eye(K);		% account for the correction
%  g = G./((sum(psiy.^2)/n)'*diag(cv)');
  omega = (sum(psiy.^2)/n)'*diag(cv)';
  g = (omega'.*G - G')./(omega.*omega' - 1);

  % Non linear part of the gradient

  psiy = psiy*B;			% psiy no longer needed
  for k = 1:K
    psiy(rx(:,k),k) = psiy(:,k);
  end
%  dz = -cumsum(psiy(1:n-1,:)).*diff(z);
  dz = -cumsum(psiy(1:n-1,:));

%  figure(11);
%  for k = 1:K
%    subplot(2,K,k)
%    stairs(z(:,k), [dz(:,k); dz(n-1,k)]/n); grid on
%    hold on; plot(z(:,k), p(:,k),'r'); hold off
%  end
  dz = dz.*diff(z);

  m = sum(dz)/(n-1);
  dz = dz - m(ones(n-1,1),:);		% centering
  for k = 1:K
    T = basis./p(:,k(ones(1,np)));
    T = lamda(:,k(ones(1,np))).*T(1:n-1,:) + ...
	(1-lamda(:,k(ones(1,np)))).*T(2:n,:);
    dp0(:,k) = (T'*T)\(T'*dz(:,k));
    decr(k) = (dz(:,k)'*T)*dp0(:,k)/n;
  end

%  figure(11);
%  for k = 1:K
%    subplot(2,K,k+K);
%    stairs(z(:,k), [dz(:,k); dz(n-1,k)]); grid on
%  end

  slope = -(mu_B*sum(sum(G.*g)) + mu_z*sum(decr));

  % Transformation (trial until the criterion decreases)
  % adapted from "Line search and backtrack" of "numer. recipe."

  alamin = sqrt(TOLX2*K*(K+np-2)/(-slope));

  alam = 1.0;				% full step first
  while 1				% begin line search loop
    if (alam < alamin);			% convergent or step becomes too small
      for k = 1:K
        z(:,k) = z(rx(:,k),k);
      end
      return;
    else
    % new sources and the criterion
    
      Bn = B - (alam*mu_B)*g*B;			% linear part
%      p0n = p0.*(1 + tanh(alam*mu_z*dp0./p0));
      p0n = p0.*exp(alam*mu_z*dp0./p0);

      p = basis*(p0n);

      %compute the quantile z from p, the auxilary array lambda and the entropy
      z = diff(p);				% save diff(p) temporally
      lamda = z./p(1:n-1,:);
      non0 = abs(lamda) > 1e-5;
      z(non0) = log(1+lamda(non0))./(z(non0)*n);
      z(~non0) = (1/n)./p(~non0);
      m = sum(log(z))/(n-1);
      z = z.*exp(entz - m(ones(n-1,1),:));	% rescale
      lamda(non0) = 1./log(lamda(non0)+1) - 1./lamda(non0);
      lamda(~non0) = .5;
      z = [zeros(1,K); cumsum(z)];
      m = sum(z)/n;
      z = z - m(ones(n,1),:);			% center

      for k = 1:K				% reconstructed sources
        y(:,k) = z(rx(:,k),k);
      end
      y = y*Bn';

      cv = y'*y/n;
      for k = 1:K
        [enty(k), psiy(:,k)] = kentropy(y(:,k), bin*sqrt(cv(k,k)),filt);
      end

      crit = sum(enty) - K*(entz+log(n-1)+EULER) - log(abs(det(Bn)));

      if (crit <= crit0 + ALF*alam*slope);	% sufficient decrease
        crit0 = crit;
        B = Bn;
	if nargout > 3; Bs = [Bs; B]; end
        p0 = p0n;

	%fprintf('%2d Criterion %.5g\n', iter, crit);
        break;					% JUMP OUT OF LINE SEARCH here
      else
        if (alam == 1.0)
          tmplam = - slope/(2*(crit-crit0-slope));	% first time
        else
          rhs1 = crit - crit0 - alam*slope;
          rhs2 = crit2 - crit0 - alam2*slope;
          a = (rhs1/(alam*alam) - rhs2/(alam2*alam2))/(alam - alam2);
          b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
          if (a == 0.0)
            tmplam = -slope/(2.0*b);
          else
            disc = b*b-3.0*a*slope;
            if (disc < 0)
              error('Roundoff problem in line search');	% exit automatically
            else
              tmplam=(-b+sqrt(disc))/(3.0*a);
            end		% if (disc < 0)
            if (tmplam > 0.5*alam), tmplam = 0.5*alam; end
          end			% if (a == 0.0)
        end                   % if (alam == 1.0)
        alam2 = alam;
        crit2 = crit;
        alam = max([tmplam 0.1*alam]);
	%fprintf('crit = %.5g => alam = %.4g\n', crit, alam);
      end                     % if (crit <= crit0 ...
    end                       % if (alam < alamin)       
  end                         % while (end of line search loop)

% this part is to follow the evolution of the algorithm (2 sources only)

%  figure(12);
%  plot(z(rx(:,1),1), z(rx(:,2),2), '+'); grid on
%  mini = min(y); maxi = max(y);
%  A = inv(B);
%  line([ [mini(1) 0; maxi(1) 0]*A(1,:)' ...
%	[0 mini(2); 0 maxi(2)]*A(1,:)'], ...
%       [ [mini(1) 0; maxi(1) 0]*A(2,:)' ...
%	[0 mini(2); 0 maxi(2)]*A(2,:)'])

%  pause;

end                           % for iter = ... (end of iteration loop)

for k = 1:K
   z(:,k) = z(rx(:,k),k);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Internal function %%%%%%%%%%%%%

function [entropy, psi] = kentropy(data, bin, filt)

% Internal function: ALL INPUTS MUST BE SPECIFIED CORECTLY, NO CHECK IS DONE.
%
% Estimate the entropy of a random variable and its score function
% based on the partial derivatives of the former with respect to the
% data values. The entropy is estimated through the kernel density
% estimator with the kernel being the density of the sum of 3
% independent uniform random variables in [-1/2,1/2]. The kernel bandwidth
% is bin and the density is evaluated at bin apart. If filt is not equal
% to 1, it specifies the post smoothing.
%
% data must be a column vector (assumed centered), the output psi is
% also a column vector (not centered, not normalised)

n = length(data);

% Grouping the data into cells, index gives the index of the cell
% containing a datum, r gives its relative distance to the leftmost
% border of the cell

r = data/bin;				% column vector here
index = floor(r);
r = r - index;
m = min(index);
lfilt = length(filt);
lenf = max(index) - m + 2+lfilt;   % anticipated length of the density vector
index = index - m + 2;			% 2 <= index <= lenf - lfilt

% Compute the probability at grid cell, a data at distance r to i
% (0 <= r <= 1) contribute
%	(1-r)^2/2		to the cell [i-1,i]
%	1/2 + r(1-r)		to the cell [i,i+1]
%	r^2/2			to the cell [i+1,i+2]

f = full(sparse([index-1; index; index+1], 1, ...
                [(1-r).^2/2; .5 + r.*(1-r); r.^2/2], lenf, 1));
f = filter(filt, n, f);			% filter and division by n

% compute the entropy

logf = log(max(f,realmin));		% avoid zero
entropy = -sum(f.*logf) + log(bin);

% Compute the score:

logf = filter(filt, bin, logf);		% filter and division by bin
logf(1:lfilt-1) = [];

psi = (1-r).*logf(index-1) + (2*r-1).*logf(index) - r.*logf(index+1);
