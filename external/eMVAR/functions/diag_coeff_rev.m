%% ESTIMATES STRICTLY CAUSAL MVAR MODEL GIVEN EXTENDED MVAR MODEL

% output: Am=[A(1)...A(p)], M*pM matrix of strictly causal MVAR coeffs
% output: Su, input covariance matrix of strictly causal MVAR model

% input: Bm=[B(1)...B(p)], M*pM matrix of extended MVAR coeffs
% input: B0, M*M  matrix of instantaneous effects (must have unitary diagonal)
% input: Sw, input covariance matrix of extended MVAR model (must be diagonal)

function [Am, Su]=diag_coeff_rev(Bm,B0,Sw);

M=size(Bm,1);
p=size(Bm,2)/M; %MVAR model order

L=inv(eye(M)-B0); % Eq. 20
Su=L*Sw*L'; % Eq. 21

Ak=NaN*ones(M,M,p); Bk=Ak;
for k=1:p
    Bk(:,:,k)=Bm(:,(k-1)*M+1:k*M);
    Ak(:,:,k)=L*Bk(:,:,k); % Eq. 20
end

Am=[];
for k=1:p
    Am=[Am Ak(:,:,k)];
end

