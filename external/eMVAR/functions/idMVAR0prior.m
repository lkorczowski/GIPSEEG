%% IDENTIFICATION OF EXTENDED MVAR MODEL: Y(n)=B(0)Y(n)+B(1)Y(n-1)+...+B(p)Y(n-p)+W(n)
% MAKES USE OF PRIOR KNOWLEDGE ABOUT DIRECTION OF INSTANTANEOUS EFFECTS

%%% input:
% Y, M*N matrix of time series (each time series is in a row)
% p, model order
% ki, row vector of causal ordering
% idMode, determines strcitly causal estimation algorithm (0:builtin least squares, else other methods [see mvar.m from biosig package])

%%% output:
% Am=[A(1)...A(p)], M*pM matrix of the estimated coeffs of strictly causal MVAR model
% B0, M*M estimated matrix of instantaneous effects
% Bm=[B(1)...B(p)], M*pM matrix of the estimated coeffs of extended MVAR model
% Su, Sw, estimated M*M input covariance matrices
% Up, Wp estimated innovations

function [Bm,B0,Sw,Am,Su,Up,Wp]=idMVAR0prior(Y,p,ki,idMode);

[M,N]=size(Y);

error(nargchk(1,4,nargin));
if nargin < 4, idMode=0; end % default least squares identification
if nargin < 3, ki=[1:M]'; end % default series order: unaltered
if nargin < 2, p=10; end % default model order

%% strictly causal MVAR model describing the data
[Am,Su,Yp,Up]=idMVAR(Y,p,idMode);
% Up=Y-Yp; Up=Up(:,p+1:N); % residuals of strictly causal model

%% permutations according to causal ordering
P=zeros(M,M); %permutation matrix according to causal order ki
for i=1:M
    P(i,ki(i))=1;
end
% Upt=P*Up; % permute according to causal order ki
Sut=P*Su*P'; %Su tilda (reordered)

%% Estimate B0 through Cholesky decomposition of input covariance
[Lt,Swt]=choldiag(Sut);
% Lt1=inv(Lt);

%% return to original causal ordering and get outputs of instantaneous model
L=P'*Lt*P;
B0=eye(M)-inv(L); %instantaneous effects
Wp=(eye(M)-B0)*Up; % residuals of extended model


%% estimate lagged effects and residual covariance
Ak=NaN*ones(M,M,p); Bk=Ak;
Bm=[];
for k=1:p
    Ak(:,:,k)=Am(:,(k-1)*M+1:k*M);
    Bk(:,:,k)=(eye(M)-B0)*Ak(:,:,k);
    Bm=[Bm Bk(:,:,k)];
end

Sw=(eye(M)-B0)*Su*(eye(M)-B0)';
Sw=diag(diag(Sw)); %prune off-diagonal elements of Sw (to have exactly zero)
 


