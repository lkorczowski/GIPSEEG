% test SYNC Riem
clear all
Directory= 'D:\data\Hyperscanning\BI-multiplayers\'
load([Directory 'Training\ALLp1.mat'])
[P1 r m]=zscore(P1,2);

close all
figure
hold all
allusers=1:size(P1,3);
for user=37;
for j=1:size(P1,1)
for i=allusers
SYNC(i)=SYNC_Riem(P1(:,:,user),P1(:,:,i));
[tmp LAGS]=xcorr(P1(j,:,user),P1(j,:,i));
COR(i,j)=tmp(LAGS==0).^2;
%cov(P1(:,:,1)',P1(:,:,i)')
end
end
plot(SYNC,mean(COR.^2,2),'*')
diffSYNC(user)=SYNC(allusers==user)-max(SYNC(~(allusers==user)));
hold off
end
diffSYNC=diffSYNC'

%% compute Riemmanien distance between P1
clear all
%close all
Directory='D:\data\Hyperscanning\BI-multiplayers\Training\';
load([Directory 'ALLp1.mat'])
load([Directory 'ROC71.mat'])
AUsers=Generate_Users_Numbers(1:71);
nusers=2;
allcombi=nchoosek(AUsers,nusers);
%[P1 r m]=zscore(P1,2);
for Col=1:size(allcombi,1)
    try
        Users=allcombi(Col,:);
        IND=find(ismember(AUsers,allcombi(Col,:)));
        RND=[];
        P1_riem_dist(Col)=distance(cov(P1(:,:,IND(1))'),cov(P1(:,:,IND(2))'),'riemann');
       	Cref=blkdiag(cov(P1(:,:,IND(1))'),cov(P1(:,:,IND(2))'));
        C=cov([P1(:,:,IND(1));P1(:,:,IND(2))]');
        P1toP2_riem_dist(Col)=distance(Cref,C,'riemann'); % SYNC_riem
        for p=1:size(P1,1)
        for j=1:size(P1,1)
        [tmp LAGS]=xcorr(P1(j,:,IND(1)),P1(p,:,IND(2)));
        COR(j,p)=abs(tmp(LAGS==0));
        end
        end
        P1toP2_CorrINTRA(Col)=mean(diag(COR));
        P1toP2_CorrALL(Col)=mean(mean(COR));
        MeanCovINTRA(Col)=mean(diag(abs(P1(:,:,IND(1))*P1(:,:,IND(2))')));
        MeanCovALL(Col)=mean(mean(abs(P1(:,:,IND(1))*P1(:,:,IND(2))')));
        MeanCorrINTRA(Col)=mean(diag(abs(supressMean(P1(:,:,IND(1)))*supressMean(P1(:,:,IND(2)))')));
        MeanCorrALL(Col)=mean(mean(abs(supressMean(P1(:,:,IND(1)))*supressMean(P1(:,:,IND(2)))')));
        MeanNormCorrINTRA(Col)=mean(mean(abs(zscore(P1(:,:,IND(1)))*zscore(P1(:,:,IND(2)))')));
        MeanNormCorrALL(Col)=mean(mean(abs(zscore(P1(:,:,IND(1)))*zscore(P1(:,:,IND(2)))')));

    catch e
        disp(['error occurs at ' num2str(Col)])
        e
    end
end
%%
disp(P1toP2_CorrINTRA)
disp(P1toP2_CorrALL)
disp(MeanCovINTRA)
disp(MeanCovALL)
disp(MeanCorrALL)
disp(MeanCorrALL)
plot(LAGS,tmp)
RangeP=mean(P1toP2_riem_dist)+3*std(P1toP2_riem_dist);
RangeM=mean(P1toP2_riem_dist)-3*std(P1toP2_riem_dist);
Range=(RangeM<P1toP2_riem_dist)&(P1toP2_riem_dist<RangeP);
hist(P1toP2_riem_dist(Range))

%%
close all
figure;plot(P1_riem_dist,P1toP2_riem_dist,'*')
figure;plot(P1toP2_riem_dist.^2,P1toP2_CorrALL,'*');xlabel('SYNC RIEM^2');ylabel('Average Average Corr(0)')
figure;plot(P1toP2_riem_dist.^2,MeanCovALL,'*');xlabel('SYNC RIEM');ylabel('Average Corr(0) electrod per electod')
figure;plot(P1toP2_riem_dist.^2,MeanCorrALL,'*');xlabel('SYNC RIEM');ylabel('Average Corr(0) electrod per electod')
figure;plot(P1toP2_riem_dist.^2,MeanNormCorrINTRA,'*');xlabel('SYNC RIEM');ylabel('Average Corr(0) electrod per electod')
figure;plot(P1toP2_riem_dist.^2,MeanNormCorrALL,'*');xlabel('SYNC RIEM');ylabel('Average Corr(0) electrod per electod')

%%
clear all
close all
max=100;
for p=1:1
for i=0:99
t=0:0.1:2;
Siz=length(t);
maxt=max-i;
x=sin(t*2*pi);
%x=zscore(x);
x1=[x zeros(1,max-Siz)];
x2=[zeros(1,max/2) x zeros(1,max/2)];
x3=[zeros(1,max) x];

NoiseV=0;
At=x1;Bt=x3;A=At+NoiseV*randn(1,length(At));B=Bt+NoiseV*randn(1,length(Bt));
[R LAGS]=xcorr(A,B);
Correlation(i+1,p)=R(LAGS==0);
SYNC(i+1,p)=SYNC_Riem(A,B);
end
end

%
figure;subplot(211);plot(At);title('signals');subplot(212);plot(Bt);
figure;subplot(211);plot(A);title('signals+noise');subplot(212);plot(B);
figure;subplot(211);plot(abs(mean(Correlation,2)));title('Coor CS SYNC');subplot(212);plot(mean(SYNC,2));

%%
clear all
close all
ms=10;
mMax=2;
mt=0:mMax/ms:mMax-1/ms;
m=sin(mt*2*pi);
subplot(211);plot(m)

allMax=100;
p0=dirac((0:allMax-1));p0(p0==inf)=1;
x0=zscore(conv(p0,m));
for ite=1:50
for tau=0:allMax-1
p=dirac((0:allMax-1)-tau);
p(p==inf)=1;
x=zscore(conv(p,m));
NoiseV=2;
At=x0;Bt=x;A=At+NoiseV*randn(1,length(At));B=Bt+NoiseV*randn(1,length(Bt));
[R LAGS]=xcorr(A,B);
Correlation(tau+1,ite)=R(LAGS==0);
SYNC(tau+1,ite)=SYNC_Riem(A,B);
end
subplot(212);plot(x)
end
figure;subplot(211);plot(At);title('signals');subplot(212);plot(Bt);
figure;subplot(211);plot(A);title('signals+noise');subplot(212);plot(B);
figure;subplot(211);plot(abs(mean(Correlation,2)));title('Corr VS SYNC');subplot(212);plot(mean(SYNC,2));
