%% PORTMANTEAU TEST FOR AUTOCORRELATION [from Luetkepohl 2005, New introduction to multiple time series analysis]
% null hypothesis: joint whiteness
% if pp<alpha, reject -> non white

function [pp,Qh,critlo,crithi,crit1tail,C0,stringout,stringflag]=test_whiteness(u,p,h,alpha);

% %%% input: u, p, alpha, h
% u: matrice dei residui (in riga!!!)
% h=10; %number of lags to test
% alpha=0.05; %significance
% p=7; % model order

% u=u'; %residui in riga

T=size(u,2); K=size(u,1); % notazione di Lutkepohl


C=zeros(K,K,h+1); %il primo è C0, l'(h+1)esimo è Ch

%% stima matrici di correlazione dei residui
% LAG 0
C0=zeros(K,K);
for t=1:T
    C0=C0+u(:,t)*u(:,t)';
end
C0=C0./T;

%LAG 1,...,h
C=zeros(K,K,h);
for i=1:h
    for t=i+1:T
        C(:,:,i)=C(:,:,i)+u(:,t)*u(:,t-i)';
    end
    C(:,:,i)=C(:,:,i)./T;
%     C(:,:,i)=C(:,:,i)./(T-i);
end

%% portmanteau statistic (extension of Ljung Box test to multivariate time series)
Qh=0;
for j=1:h
    Qh=Qh+trace(C(:,:,j)'*inv(C0)*C(:,:,j)*inv(C0)) / (T-j);
end
Qh=Qh*T*(T+2);

dg=K*K*(h-p);

critlo=chi2inv(alpha/2,dg);
crithi=chi2inv(1-alpha/2,dg);
crit1tail=chi2inv(1-alpha,dg);

pp=1-chi2cdf(Qh,dg); %se alpha/2<pp<1-alpha/2, don't reject -> OK is white!

if pp<0.05
    stringout='rejection: signals are NOT WHITE';
    stringflag=1;
else
    stringout='non-rejection: signals are WHITE';
    stringflag=0;
end
