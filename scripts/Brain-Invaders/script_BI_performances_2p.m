%% LOAD FILES WITH MAT
clear all;
close all
%%for 1 players from the training data
nb_players=71;
AUsers = Generate_Users_Numbers(1:nb_players);
%Users = {'11','15','18','21'}; % selected players
Directory= 'D:\data\Hyperscanning\BI-multiplayers\'
load([Directory 'Training\ALLsolo.mat'])
load([Directory 'BIpairs.mat'])
RNDtr={};RNDte={};

%% prepare data
    Xte=ALLdata.Xte;
    Yte=ALLdata.Yte;
 nusers=size(allcombi,2);
% allcombi=nchoosek(AUsers,nusers);
auc=[];
se=1; %selection session for users data 
%% generate multiplayers set
for Trial=1:10
for i=1:size(allcombi,1)
IND=find(ismember(AUsers,allcombi(i,:)))
%%generate multiplayers set
RND=[];P=0.25;
        Xte=ALLdata.Xte(IND);Yte=ALLdata.Yte(IND);
        [Xout, Yout,RND]=Multisubject_generator(ALLdata.Xte(IND),ALLdata.Yte(IND),RND);
        [X, Y]=Randomize_trials(Xout, Yout);
        [Xtraining Ytraining Xtest Ytest]=Generate_Training_Test_Set(X,Y,P);
        
% end
%% classification
P300_ref_orientation='multiP1'
% Estimation des matrices de covariance spéciale P300
[COVtr, P1] = covariances_p300_hyper(Xtraining,Ytraining,nusers,P300_ref_orientation);
COVte = covariances_p300_hyper(Xtest,P1,nusers,P300_ref_orientation);
Stats={'all', 'intra'}
for StatsIND=1:length(Stats)
    Stats{StatsIND}
%%% Classification par MDM riemannienne
[Yestimated, distances,C] = mdm_hyper(COVte,COVtr,Ytraining,nusers,Stats{StatsIND},P300_ref_orientation,'ld','ld'); 
%[Yestimated, distances,C] = mdm(COVte,COVtr,Ytraining); 

Perfclassif=length(find(Yestimated==Ytest))/length(Ytest);
disp(['Performance classification using closest reference :' num2str(Perfclassif) '%']);
% matrice de confusion
disp('matrice de confusion');
disp(confusionmat(Ytest,Yestimated));
% courbe roc 
scores = -diff(distances')';
 [TPR,FPR,TH] =roc(Ytest',scores');
%figure;plotroc(Ytest',scores');

% Area under curve (higer is better, max=1);
[X,Y,~,AUC{i,Trial,StatsIND},OPTROCPT] = perfcurve(Ytest',scores',1);
AUC{i,Trial,StatsIND}
figure;plot(X,Y);
       xlabel('False positive rate'); ylabel('True positive rate')
       title(['ROC for classification by logistic regression Couple=' cat(2,allcombi{i,:})])
               text(0.5,0.5,['AUC=' num2str(AUC{i,Trial,StatsIND},4)],'FontSize',12)
               close all
end
%auc = area_under_curve(scores',Ytest')
%disp(['Area Under Curve : ' num2str(auc)]);
end
end
save( [Directory 'results_AUC_STATS_all_VS_intra.mat'],'AUC')
%% load all data
clear all
Stats={'all', 'intra'}

nb_players=71;
AUsers = Generate_Users_Numbers(1:nb_players);
Directory='D:\data\Hyperscanning\BI-multiplayers\';
load([Directory 'BIpairs.mat']);
load( [Directory 'results_AUC_STATS_all_VS_intra.mat'])
Directory='D:\data\Hyperscanning\BI-multiplayers\Training\';
load([Directory 'ALLp1.mat'])
load([Directory 'ALLsolo.mat'])
[P1 r m]=zscore(P1,2);
for Col=1:length(allcombi) % compute SYNC
    try
        Users=allcombi(Col,:);
        IND=find(ismember(AUsers,allcombi(Col,:)));
        P1_riem_dist(Col)=distance(cov(P1(:,:,IND(1))'),cov(P1(:,:,IND(2))'),'riemman');
        P1toP2_riem_dist(Col)=SYNC_Riem(P1(:,:,IND(1)),P1(:,:,IND(2)));
    catch e
        disp(['error occurs at ' num2str(Col)])
        e
    end
end
%% analysis
close all
Groups=[];
for i=1:size(allcombi,1)
Groups{i}=cat(2,allcombi{i,:});
end
AUCn=cell2mat(AUC);
AUCmean=mean(AUCn,2);
figure;plot(squeeze(AUCmean(:,1,1)))
hold all;plot(squeeze(AUCmean(:,:,2)));hold off;
set(gca, 'XTick', 1:length(Groups),'XTickLabel', Groups);
legend(Stats)
ylabel('ROC')
title('Overview of group performances')

figure
diffAUC=AUCmean(:,1,1)-AUCmean(:,1,2);
plot(P1_riem_dist,P1toP2_riem_dist,'*')