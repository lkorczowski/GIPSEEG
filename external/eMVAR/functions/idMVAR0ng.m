%% IDENTIFICATION OF EXTENDED MVAR MODEL: Y(n)=B(0)Y(n)+B(1)Y(n-1)+...+B(p)Y(n-p)+W(n)
% MAKES USE OF ICA (no need of prior knowledge!)

%%% input:
% Y, M*N matrix of time series (each time series is in a row)
% p, model order
% idMode, determines strcitly causal estimation algorithm (0:builtin least squares, else other methods [see mvar.m from biosig package])

%%% output:
% Am=[A(1)...A(p)], M*pM matrix of the estimated coeffs of strictly causal MVAR model
% B0, M*M estimated matrix of instantaneous effects
% Bm=[B(1)...B(p)], M*pM matrix of the estimated coeffs of extended MVAR model
% Su, Sw, estimated M*M input covariance matrices
% Up, Wp, estimated innovations
% ki, estimated causal order for instantaneous effects

function [Bm,B0,Sw,Am,Su,Up,Wp,percup,ki]=idMVAR0ng(Y,p,idMode);

[M,N]=size(Y);

error(nargchk(1,4,nargin));
if nargin < 4, idMode=0; end % default least squares identification
if nargin < 3, ki=[1:M]'; end % default series order: unaltered
if nargin < 2, p=10; end % default model order

%% 1) strictly causal MVAR model describing the data
[Am,Su,Yp,Up]=idMVAR(Y,p,idMode);
% Up=Y-Yp; Up=Up(:,p+1:N); % residuals of strictly causal model


%% 2) LiNGAM analysis on the residuals U
%%% a) ICA on residuals
[Sica, Mm, Q] = fastica( Up, 'approach', 'symm', 'g', 'tanh', 'epsilon', 1e-14, 'displayMode', 'off');  

%%% b) row-permutation of Q to avoid zeros on the main diagonal
if M<=8
    [Qbar,rigaperm] = permnozeribrutal ( Q ); % Try all row permutations, find best solution
else
    [Qbar,rigaperm] = permnozerihungarian( Q ); % Find best row permutation by hungarian algorithm
end

%%%c) Divide each row of Qbar by the diagonal element and find B0
Qtilda = Qbar./(diag(Qbar)*ones(1,M));
B0=eye(M)-Qtilda;

% This is to give an idea of lower triangularity of permuted B0, i.e. of acyclicity of B0
if M<=8
    [Bopt,ki,bestvalore] = permslowertriagbrutal( B0 );% Identically permute the rows and columns of Qtilda so as to get an approximately lower triangular matrix
else
    [Bopt,ki] = sltprune( B0 );
end
D = tril(Bopt,-1)-Bopt; scoreuptri = sum(sum(D.^2)); % strictly lower triangular score
percup = scoreuptri/sum(sum(Bopt.^2));


%% 3) estimate lagged effects and residual covariance
Ak=NaN*ones(M,M,p); Bk=Ak;
Bm=[];
for k=1:p
    Ak(:,:,k)=Am(:,(k-1)*M+1:k*M);
    Bk(:,:,k)=(eye(M)-B0)*Ak(:,:,k);
    Bm=[Bm Bk(:,:,k)];
end

Sw=(eye(M)-B0)*Su*(eye(M)-B0)';
Sw=diag(diag(Sw)); %prune off-diagonal elements of Sw (to have exactly zero)

Wp=(eye(M)-B0)*Up; % residuals of extended model

