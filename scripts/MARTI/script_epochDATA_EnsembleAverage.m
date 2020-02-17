%%% (0) PREPARE Load the data player, filter, preprocessing
clear all
Directory='D:\data\Hyperscanning\MARTI\Groups\' %office
% Directory='E:\data\MARTI\' %home

load(['MARTI_GroupsName.mat'])
LABELS=getfield(load(['MARTI_ElectrodesName.mat']),'MARTImapping');

for indGroup=1:length(GroupsName)
    for indPlayer=1:2
        try
            clear X Xf Xcont Trigger TriggerC
            load([Directory 'mat\epoch\' GroupsName{indGroup} '_' num2str(indPlayer) '.mat'],'EEGFT')
            
            X=reshape(cell2mat(EEGFT.trial),size(EEGFT.trial{1},1),size(EEGFT.trial{1},2),length(EEGFT.trial));
%             cond=unique(EEGFT.Conditions)';
%             for indCond=1:length(cond)
%                 [P(:,:,:,indCond,indGroup,indPlayer) Class]=EnsembleAverage(X(:,:,EEGFT.Conditions==cond(indCond)),EEGFT.sampleinfo(EEGFT.Conditions==cond(indCond)));
%                 LengthCont=max(EEGFT.trialinfo(:,2))+EEGFT.fsample*2;
%                 Trigger=Position2Trigger(EEGFT.trialinfo(EEGFT.Conditions==cond(indCond),1),LengthCont);
%                 [Xcont TriggerC]=epoch2continuous(X(:,:,EEGFT.Conditions==cond(indCond)),Trigger,'compact');
%                 [Pover(:,:,:,indCond,indGroup,indPlayer)]=meanOverlap(Xcont,TriggerC,size(X,2),EEGFT.sampleinfo(EEGFT.Conditions==cond(indCond)));
% 
%             end
            tic
            for indK=1:size(X,3)
                Xf(:,:,indK)= preprocessingEEG(X(:,:,indK)',EEGFT.fsample,[1 20 2 0])';
            end
            
%             for indCond=1:length(cond)
%                 [Pf(:,:,:,indCond,indGroup,indPlayer) Class]=EnsembleAverage(Xf(:,:,EEGFT.Conditions==cond(indCond)),EEGFT.sampleinfo(EEGFT.Conditions==cond(indCond)));
%                 [Xcont TriggerC]=epoch2continuous(Xf(:,:,EEGFT.Conditions==cond(indCond)),Trigger,'compact');
%                 [Pfover(:,:,:,indCond,indGroup,indPlayer)]=meanOverlap(Xcont,TriggerC,size(X,2),EEGFT.sampleinfo(EEGFT.Conditions==cond(indCond)));
% 
%             end
            time=0:1/EEGFT.fsample:3-1/EEGFT.fsample
            prestim=0.5;poststim=2;
            time=find(time>=prestim,1): find(time<poststim,1,'last');
            Xte{indGroup,indPlayer}=Xf(:,time,:);
            Cond{indGroup,indPlayer}=EEGFT.Conditions;
            Yte{indGroup,indPlayer}=EEGFT.sampleinfo
            dimord='chan_time_class_cond_group_subj'
%             save([Directory 'mat\epoch\EnsembleAverage.mat'],'P','Pover','Pf','Pfover','cond','GroupsName','Class','dimord','LABELS')
            save([Directory 'mat\epoch\MARTI_epochs.mat'],'Xte','Yte','Cond','GroupsName','dimord','LABELS','prestim','poststim')

            toc
        catch e
            warning(['Group' num2str(indGroup)])
        end
    end
end
% 
size(P)
%% show the data
clear all;close all;
Directory='D:\data\Hyperscanning\MARTI\Groups\' %office
% Directory='E:\data\MARTI\' %home

load(['MARTI_GroupsName.mat'])
LABELS=getfield(load(['MARTI_ElectrodesName.mat']),'MARTImapping');
load([Directory 'mat\epoch\EnsembleAverage.mat'])

Col=[10 5 5;1 5 1;5 10 10; 5 3 10]/10;
indClass=2
     figure
     
% pp = spline(1:384,Pover(28,:,2,1,1,1))

 for condition=1:4

 subj2take=[1 2 3 4 5 6 7 8 11 12 13 14 16 18 19 20 21];
 elc2take=1:32;
 elc2take(8)=[];
 elc2take=15

 Pbas=P-repmat(mean(P(:,64:128,:,:,:,:),2),[1,size(P,2),1,1,1,1]);
 Pbasover=Pover-repmat(mean(Pover(:,64:128,:,:,:,:),2),[1,size(Pover,2),1,1,1,1]);
  
   subplot(231);hold all
plotEEG(mean(mean(P(elc2take,:,indClass,condition,subj2take,:),6),5),4,128,LABELS(elc2take),12,Col(condition,:));
title('AEA')


 subplot(234);hold all
plotEEG(mean(mean(Pover(elc2take,:,indClass,condition,subj2take,:),6),5),4,128,LABELS(elc2take),12,Col(condition,:));
title('regression AEA')

subplot(232);hold all
plotEEG(mean(mean(Pbas(elc2take,:,indClass,condition,subj2take,:),6),5),4,128,LABELS(elc2take),12,Col(condition,:));
title('AEA + baseline')

 subplot(235);hold all
plotEEG(mean(mean(Pbasover(elc2take,:,indClass,condition,subj2take,:),6),5),4,128,LABELS(elc2take),12,Col(condition,:));
title('regression AEA + baseline')

subplot(233);hold all
plotEEG(mean(mean(Pf(elc2take,:,indClass,condition,subj2take,:),6),5),4,128,LABELS(elc2take),12,Col(condition,:));
title('AEA (1-20Hz)')

subplot(236);hold all
plotEEG(mean(mean(Pfover(elc2take,:,indClass,condition,subj2take,:),6),5),4,128,LABELS(elc2take),12,Col(condition,:));
title('AEA (1-20Hz) + baseline')

 end
 hold off

