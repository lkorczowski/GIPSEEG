%% LOAD FILES WITH MAT
clear all;
close all
%%for 1 players from the training data
nb_players=71;
AUsers = Generate_Users_Numbers(1:nb_players);
%Users = {'11','15','18','21'}; % selected players

Directory= 'D:\data\Hyperscanning\BI-multiplayers\Groups\'
SAVEPATH=[Directory 'Rejected_parts\Artefacts.mat'];
load(SAVEPATH)
load([Directory 'ALLgroups.mat'])
load([Directory 'ALLrawGroupsData.mat'])
load([Directory 'Groups.mat'])
RNDtr={};RNDte={};
method_mean = 'ld';
    method_dist = 'riemann';
%%prepare data
    Xte=ALLdata.Xte;
    Yte=ALLdata.Yte;
 nusers=2;
% allcombi=nchoosek(AUsers,nusers);
auc=[];
se=1; %selection session for users data 
%% generate multiplayers set loop and save
BImapping=[1,2,5,3,6,12,14,16,21,22,23,24,25,27,28,29];
BImap={[1,2,5,3,6,12,14,16,21,22,23,24,25,27,28,29],... %full BI mapping (16)
    [5,3,6,12,14,16,21,22,23,24,25,27,28,29],... %remove Fp1,Fp2 (14)
    [5,6,12,14,16,21,22,23,24,25,27,28,29],... %remove Fp1,Fp2, AFz (13)
    [5,6,14,21,22,23,24,25,27,28,29],...%remove Fp1,Fp2, AFz, T7,T8 (11)
    [5,6,14,22,23,24,27,28,29],... %remove Fp1,Fp2, AFz, T7,T8,P7,P8 (9)
    [5,6,13,14,15,22,23,24,27,28,29],... %remove Fp1,Fp2, AFz, T7,T8,P7,P8, add C3 C4 (11)
    [22,23,24,27,28,29]};% just P3,Pz,P4,01,0z,O2
for MAP=1:length(BImap)
BImapping=[BImap{MAP} BImap{MAP}+32];
Stats={'all', 'intra'}
for i=1:size(ALLgroups,1)
%%generate multiplayers set
PourTest=25;
P=PourTest/100;X=ALLdata.Xte{i}(BImapping,:,:);Y=ALLdata.Yte{i};
NoArte=~ALLdata.isBad{i};
        [Xtraining Ytraining Xtest Ytest]=Generate_Training_Test_Set(X(:,:,NoArte),Y(NoArte),P);
        %      [Xtraining Ytraining Xtest Ytest]=Generate_Training_Test_Set(X,Y,P);
% end
%% classification
P300_ref_orientation='multiP1'
% Estimation des matrices de covariance spéciale P300
[COVtr, P1] = covariances_p300_hyper(Xtraining,Ytraining,nusers,P300_ref_orientation);
%eegplot([P1(:,:,1);P1(:,:,2)],'srate',128,'winlength',1)
    %return
COVte = covariances_p300_hyper(Xtest,P1,nusers,P300_ref_orientation);
for StatsIND=1:length(Stats)
    disp(['Classification initialized on Group ' ALLgroups{i} ' with stastistics='])
    disp(Stats{StatsIND})
%%% Classification par MDM riemannienne
[Yestimated, distances,C,COV] = mdm_hyper(COVte,COVtr,Ytraining,nusers,Stats{StatsIND},P300_ref_orientation,method_mean,method_dist); 
%[Yestimated, distances,C] = mdm(COVte,COVtr,Ytraining); 

Perfclassif=length(find(Yestimated==Ytest))/length(Ytest);
disp(['Performance classification using closest reference :' num2str(Perfclassif) '%']);
% matrice de confusion
disp('matrice de confusion');
disp(confusionmat(Ytest,Yestimated));
% courbe roc 
scores{i,StatsIND,MAP} = -diff(distances')';
 [TPR,FPR,TH] =roc(Ytest',scores{i,StatsIND,MAP}');

% Area under curve (higer is better, max=1);
[PerfX,PerfY,~,AUC{i,StatsIND,MAP},OPTROCPT] = perfcurve(Ytest',scores{i,StatsIND,MAP}',1);
AUC{i,StatsIND,MAP}
figure;plot(PerfX,PerfY);
       xlabel('False positive rate'); ylabel('True positive rate')
       title(['ROC for classification by logistic regression Couple=' cat(2,ALLgroups{i,:})])
               text(0.5,0.5,['AUC=' num2str(AUC{i,StatsIND},4)],'FontSize',12)
               close all
end
%auc = area_under_curve(scores',Ytest')
%disp(['Area Under Curve : ' num2str(auc)]);

% compute SYNC variable
SYNC.riem(i,MAP)=SYNC_Riem(supressMean(P1(:,:,1)),supressMean(P1(:,:,2)));
SYNC.P1_riem_dist(i,MAP)=distance(cov(P1(:,:,1)'),cov(P1(:,:,2)'),'riemann');
SYNC.FroNormP1P2(i,MAP)=distance(blkdiag(cov(P1(:,:,1)'),cov(P1(:,:,2)')),cov([P1(:,:,1);P1(:,:,2)]'),'fro');%fro norm
end
end
save( [Directory 'results_AUC_STATS_all_VS_intra_vsMAPPING_P0_' num2str(PourTest) '.mat'],'AUC','scores','SYNC','BImap')
%% load all data
clear all
close all
Stats={'all', 'intra'}
Directory= 'D:\data\Hyperscanning\BI-multiplayers\Groups\'

load([Directory 'Groups.mat']);
load([Directory 'GroupsSYNC.mat']);
%load( [Directory 'results_AUC_STATS_all_VS_intra_P0_25.mat'])
load( [Directory 'results_AUC_STATS_all_VS_intra_vsMAPPING_P0_25.mat'])

%GroupSYNC
%%analysis
%close all
Groups=[];
Groups=ALLgroups;
AUCn=cell2mat(AUC);
AUCmean=AUCn;
%% Analyse mapping
plot(squeeze(mean(AUCmean(:,1,:),1)))

%% Analyse intra VS inter
figure;subplot(211);bar(squeeze(AUCmean(:,1)),1,'g')
axis([0 19 min(min(AUCmean)) 1])
hold all;bar(squeeze(AUCmean(:,2)),0.5,'r');hold off;
set(gca, 'XTick', 1:length(Groups),'XTickLabel', Groups);
legend(Stats)
ylabel('AUC ROC')
title({'Overview of group performances' ...
    'If Green>Red, the InterStats increase perf'...
    'If Red>Green, the InterStats decrease perf'})


diffAUC=AUCmean(:,1)-AUCmean(:,2);
x=SYNC.riem;y=diffAUC;
subplot(212);figure; plot(x,y,'*')
title(['Correlation Coefficient=' num2str(corrcoef(x',y))])
xlabel('RiemSYNC');ylabel('Diff Perf AUCinter-AUCintra')

%%permutation test
%--- set some parameters
alpha   = 0.05;     % significance level
nPerm   = 50000;    % number of permutations (i.e., size of surrogate data)
nObs    = length(AUCmean(:,1));     	% number of observations
mean_1  = 10;       % mean of first vector of observations
mean_2  = 12.5;       % mean of second vector of observations
var_1   = 5;        % variance of first vector of observations
var_2   = 5;        % variance of second vector of observations


% --- simulate two vectors of paired obervation
vec_OBS_1   = mean_1 + var_1*randn(1,nObs);
vec_OBS_2   = mean_2 + var_2*randn(1,nObs);

vec_OBS_1   = AUCmean(:,1);
vec_OBS_2   = AUCmean(:,2);

% --- perform permutation tests
if alpha < (1/nPerm)
   disp(['Not enough permutations for this significance level (' num2str(alpha) ')']);
   alpha = 1/nPerm;
   disp(['--> Significance level set to minimum possible (' num2str(alpha) ')']);
end
schuffle_index = uniqueShuffle2(nPerm, nObs, 1); % permutation matrix used to shuffle values between observation sets
diff_PERM = zeros(1,nPerm);     
% compute OBSERVATION (reference/real value) 
diff_OBS = mean(vec_OBS_2 - vec_OBS_1);  % here: average paired differences (CAN BE A SIMPLE DIFFERENCE, OR ANY CALCULATION!!!)
% compute PERMUTATION (surrogate values) 
vec_PERM_1 = zeros(1,nObs);
vec_PERM_2 = zeros(1,nObs);
for perm_ix = 1:nPerm
    % draw specific permutation using 'schuffle_index' logical indexes
	vec_PERM_1(schuffle_index(perm_ix,:))   = vec_OBS_1(schuffle_index(perm_ix,:));
    vec_PERM_1(~schuffle_index(perm_ix,:))  = vec_OBS_2(~schuffle_index(perm_ix,:));
    vec_PERM_2(~schuffle_index(perm_ix,:))  = vec_OBS_1(~schuffle_index(perm_ix,:));
    vec_PERM_2(schuffle_index(perm_ix,:))   = vec_OBS_2(schuffle_index(perm_ix,:));
    % compute surrogate surrogate value for this specific permutation
    diff_PERM(perm_ix) =  mean(vec_PERM_2 - vec_PERM_1);
end

% compute statistical significance and plot results
p_val = (1/nPerm)*(1+length(find(diff_PERM > diff_OBS)));
h_val = p_val < alpha;  
figure;dispBootstrap(diff_OBS, diff_PERM)    
