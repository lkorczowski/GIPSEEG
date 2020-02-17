clear all
Directory='D:\data\Hyperscanning\BI-multiplayers\Training\';
load([Directory 'ALLsolo.mat'])
%% check the data
for i=1:length(ALLdata.Xte)
NbSamples(i)=size(ALLdata.Xte{i},2); % check the epochs' sizes
SizeDifference(i)=-size(ALLdata.Xte{i},3)+size(ALLdata.Yte{i},1); % check the number of epochs
P1(:,:,i)=mean(ALLdata.Xte{i}(:,:,ALLdata.Yte{i}(1:size(ALLdata.Xte{i},3))==1),3);
NaNs(i)=length(find(isnan(ALLdata.Xte{i})));
end


Nb_errors=length(find(NbSamples~=128))+length(find(SizeDifference))
indNaN=find(NaNs)
%%
close all
eegplot(P1(:,:,16),'winlength',10,'srate',128)