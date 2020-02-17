clear all
%close all
Directory='D:\data\Hyperscanning\BI-multiplayers\Training\';
load([Directory 'ALLp1.mat'])
%%
tic
IND=[1 2];
for p=1:size(P1,1)
        for j=1:size(P1,1)
        [tmp LAGS]=xcorr(P1(j,:,IND(1)),P1(p,:,IND(2)));
        COR(j,p)=abs(tmp(LAGS==0));
        end
end
mean(mean(COR))
        toc
        
        %figure;plot(tmp)
        %%
        tic
                SYNC_Riem(P1(:,:,IND(1)),P1(:,:,IND(2)))

        toc
