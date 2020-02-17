function  [sep, source, seps] = icansng(data, sep, nbloc, maxiter, eps)

% Syntaxe	[sep, source, seps] = icansng(data, sep, nbloc, maxiter)
% Separation of non stationary non gaussian sources, batch method
% Use linear combinations of s, s^3 as separating function
%
%  - data and source are m x n matrices, containing m mixtures and
%    m reconstructed sources.
%  - sep in input is the initial separating matrix
%    (default = identity matrix).
%  - sep in output is the estimated separating matrix.
%  - nbloc is the number of blocks in each of which the sources are
%    considered to be stationary (default = 3)
%  - maxiter = number of iteration, default = 20
%  - eps: (squared) stoping threshold, default 1e-2*K(K-1)/n
%  - seps if present contains the sequences of separating matrices.

[m, n] = size(data);

if nargin < 2, sep = eye(m); end
if nargin < 3, nbloc = 3; end
if nargin < 4, maxiter = 20; end
if nargin < 5, eps = 1e-4*m*(m-1)/n; end

if n < 5*nbloc, error('data too short'), end

% Initialisation: Compute the boundary of blocks and the initial sources

bloc = round(linspace(0, n, nbloc+1));		% boundary of block
source = sep*data;

% Iteration loop

%coef = zeros(2*m,nbloc);	%for the coefficients of the score functions.

if nargout > 2; seps = sep; end

for iter = 1:maxiter

  omega = zeros(m);
  G = zeros(m);

  for nb = 1:nbloc

    s = source(:,bloc(nb)+1:bloc(nb+1))';

    % Computation of moments of a block
    m0 = bloc(nb+1) - bloc(nb);
    m2 = sum(s.^2);
    m4 = sum(s.^4);
    m6 = sum(s.^6);
  
    % Computation of coefficients of the score function of a block        
    det = m6.*m2 - m4.^2;
    c = [(m0*m6 - 3*m4.*m2)./det
	 (3*m2.^2 - m0*m4)./det];
    psi = s.*c(ones(m0,1),:) + s.^3.*c(2*ones(m0,1),:);

    %coef(:,nb) = c(:);
  
    % Compute the gradient

    G = G + psi'*s;

    % Compute the hessian matrix

    psip = sum([ones(1,m); 3*m2/m0].*c);
    omega = omega + psip'*m2;

  end

  % Compute the transformed gradient

  omega = omega/n;
  G = G/n - eye(m);
  H = (omega'.*G - G')./(omega.*omega' - 1);

  % New separating matrix
  
  sep = sep - H*sep;
  if nargout > 2; seps = [seps; sep]; end
  source = source - H*source;

  if sum(sum(G.*H)) < eps; return; end

end
