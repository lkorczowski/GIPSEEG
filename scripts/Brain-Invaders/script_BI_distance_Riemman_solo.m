%% LOAD FILES WITH MAT
clear all;
close all
% select the players in the analysis
nb_players=71;
AUsers=Generate_Users_Numbers(1:nb_players);
%AUsers = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
%Users = {'11','15','18','21'}; % selected players
Directory='D:\data\Hyperscanning\BI-multiplayers\Training\';
load([Directory 'ALLsolo.mat'])
%% prepare data
nusers=1;
allcombi=nchoosek(AUsers,nusers);
%
auc=[];
P=0;
%% compute P1 and COVP1 for all users and save it 
for Col=1:length(allcombi) 
    try
        Users=allcombi(Col,:);
        IND=find(ismember(AUsers,allcombi(Col,:)));
        RND=[];
        Xte=ALLdata.Xte{IND};Yte=ALLdata.Yte{IND};
        P1(:,:,Col)= mean(Xte(:,:,Yte==1),3);
        COVp1(:,:,Col) = cov(P1(:,:,Col)');
        
    catch e
        disp(['error occurs at ' num2str(Col)])
        e
    end
end
save([Directory 'ALLp1.mat'],'P1','COVp1')
%% compute Riemmanien distance between P1
clear all
close all
Directory='D:\data\Hyperscanning\BI-multiplayers\Training\';
load([Directory 'ALLp1.mat'])
nb_players=71;
AUsers=Generate_Users_Numbers(1:nb_players);
nusers=2;
allcombi=nchoosek(AUsers,nusers);
[P1 r m]=zscore(P1,2);
for Col=1:length(allcombi) 
    try
        Users=allcombi(Col,:);
        IND=find(ismember(AUsers,allcombi(Col,:)));
        RND=[];
        P1_riem_dist(Col)=distance(cov(P1(:,:,IND(1))'),cov(P1(:,:,IND(2))'),'riemman');
       	Cref=blkdiag(cov(P1(:,:,IND(1))'),cov(P1(:,:,IND(2))'));
        C=cov([P1(:,:,IND(1));P1(:,:,IND(2))]');
        P1toP2_riem_dist(Col)=distance(Cref,C,'riemman');
    catch e
        disp(['error occurs at ' num2str(Col)])
        e
    end
end
%save([Directory 'P1_riem_dist.mat'],'P1_riem_dist','allcombi','P1toP2_riem_dist')
%% compute analysize
close all
for i=1:length(AUsers)
[IND,INDc]=find(ismember(allcombi,AUsers{i}));
AN(i,1)=mean(P1toP2_riem_dist(IND)); %mean of each individu
AN(i,2)=var(P1toP2_riem_dist(IND)); %variance of the individual performance
%AN(Col,3)=AUC(Col); %group performance

end
figure
subplot(211);plot(AN(:,1)');subplot(212);plot(AN(:,2)');
figure
plot(P1toP2_riem_dist,'*')
%%
results.AUC=AUC;
results.allcombi=allcombi;
save([Directory 'DistanceR.mat'],'results')
%}
