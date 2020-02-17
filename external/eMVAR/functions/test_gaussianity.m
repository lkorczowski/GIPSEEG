%% JARQUE-BERA TEST FOR GAUSSIANITY [from Luetkepohl 2005, New introduction to multiple time series analysis, pag 175]
% null hypothesis: joint gaussianity
% if ptot<alpha, reject -> non gaussian

function [ptot,ps,pk,lambdas,lambdak,crittresh,stringout,stringflag]=test_gaussianity(U,alpha);

% %%% input:
% u: matrice dei residui (in riga!!!)
% alpha=0.05; %significance
% u=u'; %residui in riga


%% test interno
% clear;close all;clc;
% alpha=0.05;
% T= 500; K=3; % size voluto
% % Nonlinearity exponent, selected to lie in [0.5, 0.8] or [1.2, 2.0]. (<1 gives subgaussian, >1 gives supergaussian)
% q = rand(K,1)*1.1+0.5;    
% ind = find(q>0.8);           
% q(ind) = q(ind)+0.4;     
% % This generates the disturbance variables, which are mutually independent, and non-gaussian
% U = randn(K,T); % gaussiani
% % U = sign(U).*(abs(U).^(q*ones(1,T))); % non gaussiani, commenta se li voglio gaussiani


%% vectors for Jarque-Bera analysis
T=size(U,2); K=size(U,1); % notazione di Lutkepohl

Um=mean(U')'; % media di U

Su=zeros(K,K); % covarianza di U
for t=1:T
    Su=Su+(U(:,t)-Um)*(U(:,t)-Um)';
end
Su=Su./(T-1);

Ps=chol(Su)';

for t=1:T
    V(:,t)=inv(Ps)*(U(:,t)-Um);
end

b1=zeros(K,1); b2=b1;
for k=1:K
    b1(k)=sum(V(k,:).^3) / T;
    b2(k)=sum(V(k,:).^4) / T;
end


%% statistica di Jarque-Bera (estensione a multivariate series)
lambdas=T*b1'*b1/6;
lambdak=T*(b2-3*ones(K,1))'*(b2-3*ones(K,1))/24;
dg=K;

% critlo=chi2inv(alpha,dg);
% crithi=chi2inv(1-alpha/2,dg);
crittresh=chi2inv(1-alpha,dg);

% ipotesi nulla: è gaussiana (third moment=0, fourth moment=3)

ps=1-chi2cdf(lambdas,dg); %se alpha/2<pp<1-alpha/2, don't reject -> OK is gaussian!
pk=1-chi2cdf(lambdak,dg);
ptot=1-chi2cdf(lambdas+lambdak,2*dg); %statistica totale (joint test)

if ptot<0.05
    stringout='rejection: signals are NOT GAUSSIAN';
    stringflag=1;
else
    stringout='non-rejection: signals are GAUSSIAN';
    stringflag=0;
end


% ps
% pk
% ptot
