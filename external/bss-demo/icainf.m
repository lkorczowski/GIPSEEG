function [sep, source, seps] = icainf(data, sep, maxiter, bdwidth, nb)

% Syntax:  [sep, source, seps] = icainf(data, sep, maxiter, bdwidth, nb)
%
% Separation of K linear instantaneous mixture of K sources based on
% the minimisation of the mutual information. The entropy is estimated
% via a kernel density estimator with the kernel being the density of
% the sum of 3 independent uniform random variables in [-1/2,1/2]. The
% kernel bandwidth is set to bdwidth times the optimal bandwith for
% Gausian density (= 2*(11*sqrt(pi)/(12*n))^.2 = 2.107683/n^.2) and the
% entropy is computed by discretizing at a step bandwidth/nb. 
%
% data must be a n by K matrix and sep a K by K matrix, sep, maxiter,
% bdwidth, nb are optional:
% sep at input is the the initial separation matrix (defaults to the
% identity matrix) and at output is the obtained separation matrix
% maxiter defaults to 15,
% bdwidth default 1
% nb defaults to 1
% seps, if present, is a iter*K by K matrix containing the sequence
% of iterated separation matrices, including the initial one.

% Check conformity and set default

if nargin < 1 error('At least one argument is required'); end

[n, K] = size(data);

if nargin < 2, sep = eye(K); end
if nargin < 3, maxiter = 15; end
if nargin < 4, bdwidth = 1; end
bin = bdwidth*2.107683/n^.2;
if nargin < 5
  filt = 1;
else
  filt = cumsum(1:nb-1);
  filt = [filt, (nb*(nb-1)/2+(1:nb).*(nb:-1:1)), filt(nb-1:-1:1)]/nb^3;
  bin = bin/nb;
end

ALF = 1.0e-4;
TOLX2 = 10e-4/n;

% Compute the initial transformation and the criterion

source = (data - ones(n,1)*sum(data)/n);
source = source*sep';
cova = source'*source/n;				% for use later

entropy = zeros(1,K);
psi = zeros(n,K);
for k=1:K
  [entropy(k),psi(:,k)] = kentropy(source(:,k),bin*sqrt(cova(k,k)),filt);
end
Ent = sum(entropy);

%fprintf(' 0 Entropy %.5g, criterion %.5g\n', Ent, Ent-log(abs(det(sep)))); 

% Iteration loop

if nargout > 2, seps = sep; end

for iter = 1:maxiter

  % compute gradient and the Newton step

  grad = psi'*source/n;

  omega = (sum(psi.^2)'/n - (sum(psi)'/n).^2 + ...
           (1-diag(grad).^2)./diag(cova))*diag(cova)';

  grad = grad + diag((ones(K,1) - diag(grad))./diag(cova))*cova - eye(K);

  newton = -(omega'.*grad - grad')./(omega.*omega' - 1);

  slope = sum(sum(grad.*newton));

% Line search and backtrack

%  alamin = sqrt(TOLX2*2/ ...
%      (newton(1)^2*cova(2,2)/cova(1,1) + newton(2)^2*cova(1,1)/cova(2,2)));
  alamin = sqrt(TOLX2*K*(K-1)/(-slope));

  %fprintf('Grad'); %fprintf(' %.4g', grad);
  %fprintf('\nstep'); %fprintf(' %.4g', newton);
  %fprintf('\nslope %.4g alamin %.4g\n', slope, alamin);

  alam = 1.0;			% full Newton step first
  search = 1;
  while search
    if (alam < alamin);		% convergent or step becomes too small
% comment this out if do not want plot
%     if K == 2
%	plot(data(:,1)-sum(data(:,1))/n, data(:,2)-sum(data(:,2))/n, '+');
%        A = inv(sep);
%        mini = min(source); maxi = max(source);
%        line([[mini(1) 0; maxi(1) 0]*A(1,:)' ...
%              [0 mini(2); 0 maxi(2)]*A(1,:)'], ...
%             [[mini(1) 0; maxi(1) 0]*A(2,:)' ...
%              [0 mini(2); 0 maxi(2)]*A(2,:)'])
%        grid on
%      end
      return
    else
    % new source and its entropy
      u = eye(K) + alam*newton;
      s = source*u';
      cv = u*cova*u';
      for k=1:K
        [entropy(k),psi(:,k)] = kentropy(s(:,k),bin*sqrt(cv(k,k)),filt);
      end
      EntN = sum(entropy);

      crit = EntN - log(abs(det(u)));
      if (crit <= Ent + ALF*alam*slope);	% sufficient decrease
        Ent = EntN;
        sep = u*sep;
	if nargout > 2, seps = [seps; sep]; end
        source = s;
        cova = cv;

	%fprintf('%2d Entropy %.5g, criterion %.5g\n', ...
	%	iter, Ent, Ent -log(abs(det(sep))));

        search = 0;				% get out of line search
      else
        if (alam == 1.0)
          tmplam = - slope/(2*(crit-Ent-slope));	% first time
        else
          rhs1 = crit - Ent - alam*slope;
          rhs2 = crit2 - Ent - alam2*slope;
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
            end			% if (disc < 0)
            if (tmplam > 0.5*alam), tmplam = 0.5*alam; end
          end			% if (a == 0.0)
        end			% if (alam == 1.0)
        alam2 = alam;
        crit2 = crit;
        alam = max([tmplam 0.1*alam]);
	%fprintf('crit = %.4g => alam = %.4g\n', crit, alam);
      end                     % if (crit <= Ent ...
    end                       % if (alam < alamin)       
  end                         % while search
% comment this out if do not want plot during the iteration
%  if K == 2
%    plot(data(:,1)-sum(data(:,1))/n, data(:,2)-sum(data(:,2))/n, '+');
%    A = inv(sep);
%    mini = min(source); maxi = max(source);
%    line([[mini(1) 0; maxi(1) 0]*A(1,:)' ...
%         [0 mini(2); 0 maxi(2)]*A(1,:)'], ...
%         [[mini(1) 0; maxi(1) 0]*A(2,:)' ...
%          [0 mini(2); 0 maxi(2)]*A(2,:)'])
%    grid on
%  end
end                           % for iter = ...

% should not go there ! Comment this out if do not want plot
%if K == 2
%  plot(data(:,1)-sum(data(:,1))/n, data(:,2)-sum(data(:,2))/n, '+');
%  A = inv(sep);
%  mini = min(source); maxi = max(source);
%  line([[mini(1) 0; maxi(1) 0]*A(1,:)' ...
%        [0 mini(2); 0 maxi(2)]*A(1,:)'], ...
%       [[mini(1) 0; maxi(1) 0]*A(2,:)' ...
%        [0 mini(2); 0 maxi(2)]*A(2,:)'])
%  grid on
%end

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
% is bin and the density is evaluated at bin apart. filt, if present and
% of langth > 1, specifies the post smoothing.
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
