% this script contains plenty of different steps. All the step are not
% mandatory but some are important.
%
% (1) this script load the data and compute a general parameter file for the
% different conditions to be tested
% (2) This is the computation of all different, it uses
% BCI_chain. Please consider to modify this function to add new
% conditions.
% (3)-(N) are some script to plot and check the results depending of
% different conditions. Please note that theses scripts could not work if
% the conditions which are tested have not been computed in (2).
%
%
% *** History: 2015-12-16
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015 louis[dot]korczowski[at]gmail[dot]com

%% 0 LOAD THE DATA
clear all
Directory='D:\data\Hyperscanning\EKATE\Groups\'
% Directory= 'D:\data\Hyperscanning\MARTI\Groups\'
% Directory ='F:\' %usb stick at home
% load([Directory 'mat\epoch\ALL_epochs_107.mat'])
load([Directory 'ALLgroups.mat']);load([Directory 'Groups.mat']);ALLdata.subject=ALLgroups;

load('BImap.mat')
%% 1 GENERATE ALL THE PARAMETERS
% if any(any(~[ALLdata.Yte{:,1}]==[ALLdata.Yte{:,2}])), error('Yte are not the same');else ALLdata.Yte=ALLdata.Yte(:,1); end
clear Parameters
Parameters.method_mean = {'riemann','congedo','opttransp','ld'}; % PARAMETER ld for logdet alphadivergence (Chebbi & Moakher, 2012)
Parameters.method_dist = {'riemann'}; %PARAMETER
Parameters.nusers=2; %to replace with Parameters.nusers={[1],[2],[1,2]}
Parameters.users={[1,2],1,2}%,1,2};
Parameters.elec=BImap(1); % PARAMETER #3 use 1 for the usual 16 BI electrods, 3 for the optimal subset of 13 electrodes
Parameters.samples={[1:128]}; % PARAMETER #3 use 1 for the usual 16 BI electrods, 3 for the optimal subset of 13 electrodes
Parameters.nelec=32; % number of electrodes per subject before reduction
% Parameters.Stats={'MDM-hyper', 'MDM-solo','MDM','common','SWLDA-hyper','SWLDA-multi','SWLDA-solo'};% PARAMETER #2
Parameters.classifier={'MDM'}%,'MDM','MDM-LR'}%, 'MDM-multi','MDM'};% PARAMETER #2

% Parameters.Stats={'all','SWLDA-hyper'};% PARAMETER #2

Parameters.GroupsName=ALLdata.subject;
Parameters.EIG={0}%,1,2,3,4,5,6,7,8,9,10,1:4,2:5,3:6,4:8,5:9,6:10}
Parameters.P300_ref_orientation={'multiP1'}; % PARAMETER
Parameters.RND=[]; %random seed for training/test set
% Parameters.TestSetRatio={[5 50],[10 50],[15 50],[20 50],[25 50],[30 50],[35 50],[40 50],[45 50],[50,50]};
% Parameters.TestSetRatio={0.1 0.2 0.3 0.4 0.5};
Parameters.TestSetRatio={0.05 0.15 0.25};

Parameters.Trials=1;
Parameters.ShuffleCouple={0};
Parameters.SHUF=[]; %random seed for shuffling
Parameters.SetBalance=0;
Parameters.NbMax=720;

%%generate single trial parameters
%{
             method_mean: {'ld'}
             method_dist: {'riemann'}
                  nusers: 2
                   BImap: {1x9 cell}
                   Stats: {'all'  'intra'}
              GroupsName: {19x1 cell}
                     EIG: {'all'  'limited'}
    P300_ref_orientation: {'multiP1'}
                     RND: []
%}
indP=1;
clear P
%for TrialInd=1:25
Parameters
groups2take=1:19
for MeanInd=1:length(Parameters.method_mean)
for DistInd=1:length(Parameters.method_dist)
    for GroupsInd=groups2take%1:length(Parameters.GroupsName)
        for P300ind=1:length(Parameters.P300_ref_orientation)
            for StatsInd=1:length(Parameters.classifier)
                for EIGInd=1:length(Parameters.EIG)
                    for TestTrainInd=1:length(Parameters.TestSetRatio)
                        for SchInd=1:length(Parameters.ShuffleCouple)
                            for smpInd=1:length(Parameters.samples)
                                for indS=1:length(Parameters.users)
                                    
                                    INDEX={'method_mean',MeanInd;
                                        'method_dist',DistInd;
                                        'nusers',1;
                                        'elec',1;%1 for BI normal, 3 for optimal (remove Fp1,Fp2)
                                        'classifier',StatsInd;
                                        'GroupsName',GroupsInd;
                                        'EIG',EIGInd;
                                        'P300_ref_orientation',P300ind;
                                        'RND',1;
                                        'SHUF',1;
                                        'TestSetRatio',TestTrainInd;
                                        'ShuffleCouple',SchInd;
                                        'samples',smpInd;
                                        'users',indS};
                                    tmp=generateParam(Parameters,INDEX);
                                    BOOL=1;
                                    if exist('P')
                                        for tmpIND=1:length(P)
                                            if isequal(P(tmpIND),tmp)
                                                BOOL=0;
                                            end
                                        end
                                    end
                                    if BOOL
                                        P(indP)=tmp;
                                        %                             fprintf('added ')
                                    end
                                    indP=indP+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end
fprintf('\n')
length(P)
% P=P(1);disp('warning, debug mode')
%load('D:\data\Hyperscanning\BI-multiplayers\results_Groups.mat','R')

%% 2 (1 REQUIRED) COMPUTE THE RESULTS
%close all
%P(1).RND=[];
maxTest=length(P)
R.P=P;
R.Parameters=Parameters;
MaxTrial=25;
NbGroups=length(Parameters.GroupsName);

%only for generate new seed
%if ~isfield(R,'RNDseed')
R.RNDseed=cell(MaxTrial,NbGroups);
%end
%if ~isfield(R,'SHUF')
R.SHUF=cell(MaxTrial,NbGroups);
%end
tic
t1=0;
for TrialInd=1:MaxTrial
    disp(['CV (' num2str(TrialInd) '/' num2str(MaxTrial) '). last call in ' num2str(round(toc-t1)) 's'])
                t1=toc;

    for TestInd=1:maxTest
        %         try
        P(TestInd).RND=R.RNDseed{TrialInd,strcmp(Parameters.GroupsName,P(TestInd).GroupsName)};
        P(TestInd).SHUF=R.SHUF{TrialInd,strcmp(Parameters.GroupsName,P(TestInd).GroupsName)};
        %[1AUC, 2scores,3distances, 4ConfM, 5Perfclassif,6Ytest,7RND, 8P1,9SHUF]=mdm_chain_hyper(ALLdata,Parameters,P)
        
        [   R.R(TrialInd,TestInd)...%1
            R.AUCall{TrialInd,TestInd}... %2
            R.RNDseed{TrialInd,strcmp(Parameters.GroupsName,P(TestInd).GroupsName)}... %3
            R.SHUF{TrialInd,strcmp(Parameters.GroupsName,P(TestInd).GroupsName)}]=... %4
            BCI_chain(ALLdata,Parameters,P(TestInd));
        R.P(TestInd).RND=R.RNDseed{TrialInd,strcmp(Parameters.GroupsName,P(TestInd).GroupsName)};
        R.P(TestInd).SHUF=R.SHUF{TrialInd,strcmp(Parameters.GroupsName,P(TestInd).GroupsName)};
        %         catch
        %             disp(['error occured on BCI_chain for P(' num2str(TestInd) ')'])
        %         end
        
    end
end
save([Directory 'mat\EKATE_fisher_ws_' datestr(now,30) '.mat'],'R')

toc
%% quick check 1 (figure)
clear Pall CellRes Param auc aucstd toPlot allParam
File='D:\data\Hyperscanning\EKATE\Groups\mat\EKATE_fisher_ws_20151225T111228.mat';load(File)
% File='D:\data\Hyperscanning\MARTI\Groups\mat\MARTI_fisher_ws_20151220T043219.mat';load(File)
%%
compt=0;
allParam=R.Parameters;
% create a structure with only the following parameters to plot
toPlot = {'method_mean','users','TestSetRatio','method_dist'};
for inpName=toPlot
    for indP=1:length(allParam)
        % read the acceptable names
        Param(indP).(inpName{:}) = allParam(indP).(inpName{:});
    end
end

% add the results in the P structure
auc=mean(cellfun(@mean,R.AUCall),1);
aucstd=std(cellfun(@mean,R.AUCall),[],1);
Pall=R.P;
for indP=1:length(Pall)
   Pall(indP).AUC=auc(indP);
      Pall(indP).STD=aucstd(indP);
end
figure
Res=[]
for indPl=1:length(allParam.users)  %players
    for testR=1:length(allParam.TestSetRatio) %training set ratio
        [indicesP, namesP, cvec]=struc2res(Pall,Param,{1:4,indPl,testR,2}) %struc2res extracts specific results from the structure of parameters Pall
        for indI=1:length(indicesP)
            AUC{indI,testR}=[Pall(indicesP{indI}).AUC];
            STD{indI,testR}=[Pall(indicesP{indI}).STD];
        end
    end
    subplot(1,3,indPl)
    errorbar(repmat([allParam.TestSetRatio{:}]',1,4),cellfun(@mean,AUC'),cellfun(@mean,STD')/3);legend(Param.method_mean,'location','south')
    xlabel('training set (%)')
    ylabel('AUC')
    title(['users' num2str(Param.users{indPl})])
    axis([0 0.6 0.5 1])
    
    
    
end
    set(gcf, 'PaperPosition', [0 0 30 12.5],'units','normalized','outerposition',[0.4 0.1 0.5 .9]),set(gcf,'color',[1 1 1])

        print(gcf, ['meanAUC_ratios'],'-dtiff','-r450')

%% quick check 2 (table in console)
auc=mean(cellfun(@mean,R.AUCall),1);
aucstd=std(cellfun(@mean,R.AUCall),[],1);
Param=R.P;
for indP=1:length(Param)
   Param(indP).AUC=auc(indP);
      Param(indP).STD=aucstd(indP);
end

ShowRes=rmfield(Param,{'Trials','samples','elec','EIG','P300_ref_orientation',...
    'RND','ShuffleCouple','SHUF','SetBalance','NbMax','nelec','nusers',...
    'classifier','method_dist'});struct2table(ShowRes)

%% 3 (REQUIRED FOR 4-N) LOAD RESULTS (R structure), can be long -OBSOLETE-
clear all
Directory= 'D:\data\Hyperscanning\MARTI\Groups\results\BCI\'

load([Directory 'IEEEtrans_20151029T115207.mat'])

%% 4 (3 REQUIRED) Show all the conditions mean and variances -OBSOLETE-
% we compute the mean of each condition as well as the
% close all
clear INDratio CONDITIONmean CONDITIONvar conditions_R indSTATS
P=R.P;
ALLgroups=[1,2,3,4,5,6,7,8,9,10,11,14,15,16,17,18,19];

figfolder=['D:\Mes Documents GIPSA\MATLAB\figures\BI-multiplayer\IEEEtrans15\'];
GroupsIND=[1 2 3 4 5 6 7 8 9 10 11 14 15 16 17 18 19];
Parameters=R.Parameters;
% Parameters.Stats={'all','intra','solo','SWLDA-hyper','SWLDA-multi','SWLDA-solo'}
%GROUPS=ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{7});% take only group1
%EKATE groups ton remove
GROUPS=~(ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{12})|ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{13}));
%MARTI groups to remove
GROUPS=~(ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{4})|ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{10}));

nb_conditions=length(Parameters.Stats);
nb_groups=length(ALLgroups);
nb_ratios=length(Parameters.TestSetRatio);
nb_sets=size(R.AUCall,1);
%%%%%%%%%%%%%%%%%%%%%% EACH CONDITION mean var %%%%%%%%%%%%%%%%%%%%%%%%%%
for indR=1:size(R.AUCall,2) %each condition
    conditions_R=[R.AUCall{:,indR}];
    CONDITIONmean(indR)=mean(conditions_R);
    CONDITIONvar(indR)=var(conditions_R);
end

%%%%%%%%%%%%%%%%%%%%%% EACH CONDITION indices %%%%%%%%%%%%%%%%%%%%%%%%%%
for indR=1:length(Parameters.Stats)
    indSTATS{indR}=ParametersIndex(P,'Stats',Parameters.Stats{indR});
end
%%%%%%%%%%%%%%%%%%%%%% EACH RATION indices %%%%%%%%%%%%%%%%%%%%%%%%%%
for indR=1:length(R.Parameters.TestSetRatio)
    INDratio{indR}=ParametersIndex(P,'TestSetRatio',Parameters.TestSetRatio{indR});
end

%%%%%%%%%%%%%%%%%%%%%% EACH RATION indices %%%%%%%%%%%%%%%%%%%%%%%%%%
AUCratio=[]
AUCratio_var=[]
for indR=1:length(INDratio)
    tmpmean=[];
    tmpvar=[];
    for indCond=1:length(indSTATS)
        tmpmean=[tmpmean;mean([R.AUCall{:,indSTATS{indCond}&GROUPS&INDratio{indR}}])];
        tmpvar=[tmpvar;var([R.AUCall{:,indSTATS{indCond}&GROUPS&INDratio{indR}}])];
        
    end
    AUCratio=[AUCratio,tmpmean];
    AUCratio_var=[AUCratio_var,tmpvar];
end

AUCgroups=[]
AUCgroups_v=[]
for indR=ALLgroups
    AUCgroups=[]
    AUCgroups_v=[]
    tmpmean=[];
    tmpvar=[];
    for indCond=1:length(indSTATS)
        GROUP=ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{indR});
        tmp=[];
        %tmp=[R.AUCall{:,indSTATS{indCond}&GROUP&INDratio{4}}];%specific
        %training ratio
        tmp=[R.AUCall{:,indSTATS{indCond}&GROUP&INDratio{4}}];
        if length(tmp)>nb_sets
            tmp=reshape(tmp,2,size(R.AUCall,1));
        end
        tmpmean=[tmpmean;max(tmp,[],1)];
        tmpvar=[tmpvar;var(tmp,[],2)];
        
    end
    AUCgroups=[AUCgroups,tmpmean];
    AUCgroups_v=[AUCgroups_v,tmpvar];
    lol6(indR)=testt(AUCgroups(1,:),AUCgroups(3,:))
end
figure;
Markers={'o','x','*','^','+','>','<','<'}
FontSize=24

hold all
for indR=1:size(AUCgroups,1)
    ebplot=plot(AUCgroups(indR,:)');%,sqrt(AUCvar(indR,:))');%,sqrt(VARgroups(indR,:))');
    
    %consider errorbar
    set(ebplot, 'Marker',Markers{indR},'MarkerSize',12,'LineWidth',2);
end
set(gca,'XTick', 1:nb_groups);set(gca,'XTickLabel', ALLgroups,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
hold off


%close all
figure;

set(gcf, 'color', [1 1 1])

set(gcf, 'PaperPosition', [0 0 20 16],'units','normalized','outerposition',[0.6 0.1 0.375 .9])
hold all
LEGENDS={'MDM-hyper','MDM-multi','MDM-mean','SWLDA-hyper','SWLDA-multi','SWLDA-mean'}
for indR=1:size(AUCratio,1)
    ebplot=plot((AUCratio(indR,:))','Color',[.1 .1 .1]);%,sqrt(AUCvar(indR,:))');%,sqrt(VARgroups(indR,:))');
    
    %consider errorbar
    set(ebplot, 'Marker',Markers{indR},'MarkerSize',12,'LineWidth',2);
end
%end
set(gca,'XTick', 1:nb_ratios);set(gca,'XTickLabel', {'5','10','15','20','25','30','35','40','45','50'},'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
xlabel('Number of trials per class used for training')
ylabel('mean AUC')
legend(LEGENDS,'location','SouthEast','fontsize',FontSize-10)
hold off
colormap('gray')

saveas(gcf,[figfolder 'meanAUC_ratios.tiff'])
figure
hist(sqrt(AUCgroups(:)),10)

%% 6 (5 REQUIRED) AFTER 5 COMPUTE THE BAR PLOT TO TEST BETWEEN RATIO -OBSOLETE-
close all
figure;
hold all
Markers={'o','x','*','^','+'}
for indR=1:size(MEANgroups,1)
    ebplot=plot(MEANgroups(indR,:)')%,sqrt(VARgroups(indR,:))')
    %consider errorbar
    set(ebplot, 'Marker',Markers{indR},'MarkerSize',10,'LineWidth',2)
end
set(gca,'XTick', 1:nb_ratios);set(gca,'XTickLabel', {'5','10','15','20','25','30','35','40','45','50'},'fontsize',FontSize-10,'fontname','times new roman','FontAngle','italic')
xlabel('Number of repetition per class used for training')
legend({'MDM1+MDM2','MDM-multi','MDM-solo'},'location','best')
hold off
%close all