% this script is made to load and evaluate a data set of EEG with several
% users.
% All the results should be saved in the structure 'Results'
% All the parameters can be saved in the structure 'Parameters'
%% generate global parameters TO DO BEFORE ANY ANALYSIS
clear all
Directory= 'D:\data\Hyperscanning\EKATE\Groups\'
load( [Directory 'results_AUC_STATS_all_VS_intra_vsMAPPING_P0_25.mat'],'BImap')
load([Directory 'Groups.mat'])
load([Directory 'ALLgroups.mat'])
load([Directory 'Flash.mat'])
clear Parameters
Parameters.method_mean = {'ld'}; % PARAMETER
Parameters.method_dist = {'riemann'}; %PARAMETER
Parameters.nusers=2;
Parameters.BImap=BImap; % PARAMETER #3
Parameters.Stats={'all', 'intra'};% PARAMETER #2
Parameters.GroupsName=ALLgroups;
Parameters.EIG={0}%,1,2,3,4,5,6,7,8,9,10,1:4,2:5,3:6,4:8,5:9,6:10}
Parameters.P300_ref_orientation={'multiP1'}; % PARAMETER
%Parameters.RND=[]; %random seed for training/test set
Parameters.TestSetRatio={0.875,0.75,0.625,0.5,0.25};
Parameters.Trials=1;
%Parameters.ShuffleCouple={1,0};
%Parameters.SHUF=[]; %random seed for shuffling
%Parameters.SetBalance=1;
%Parameters.NbMax=456;

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
for GroupsInd=1:length(Parameters.GroupsName)
    for P300ind=1:length(Parameters.P300_ref_orientation)
        for StatsInd=1:length(Parameters.Stats)
            for EIGInd=1:length(Parameters.EIG)
                for TestTrainInd=1:length(Parameters.TestSetRatio)
                    %     for SchInd=1:length(Parameters.ShuffleCouple)
                    INDEX={'method_mean',1;
                        'method_dist',1;
                        'nusers',1;
                        'BImap',3;
                        'Stats',StatsInd;
                        'GroupsName',GroupsInd;
                        'EIG',EIGInd;
                        'P300_ref_orientation',P300ind;
                        %'RND',1;
                        %'SHUF',1;
                        'TestSetRatio',TestTrainInd;
                        %'ShuffleCouple',SchInd
                        };
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
                            fprintf('added ')
                        end
                        indP=indP+1;
                end
            end
        end
    end
end
%end
fprintf('\n')
length(P)
%load('D:\data\Hyperscanning\BI-multiplayers\results_Groups.mat','R')

%% COMPUTE EVERYTHING (do not run if you plan to just load the results)
%close all
%P(1).RND=[];
maxTest=length(P)
R.P=P;
R.Parameters=Parameters;
MaxTrial=10;
NbGroups=length(Parameters.GroupsName);

%only for generate new seed
%if ~isfield(R,'RNDseed')
R.RNDseed=cell(MaxTrial,NbGroups);
%end
%if ~isfield(R,'SHUF')
R.SHUF=cell(MaxTrial,NbGroups);
%end
tic
for TrialInd=1:1
    
    for TestInd=1:maxTest
        
        %P(TestInd).RND=R.RNDseed{TrialInd,strcmp(Parameters.GroupsName,P(TestInd).GroupsName)};
        %P(TestInd).SHUF=R.SHUF{TrialInd,strcmp(Parameters.GroupsName,P(TestInd).GroupsName)};
        [R.AUCall{TrialInd,TestInd}...
            R.Scores{TrialInd,TestInd}...
            R.Distances{TrialInd,TestInd}...
            R.ConfM{TrialInd,TestInd}...
            R.Perf{TrialInd,TestInd}...
            R.Ytest{TrialInd,TestInd}]=...
            mdm_chain_hyper_chrono(ALLdata,Parameters,P(TestInd));
        %R.P(TestInd).RND=R.RNDseed{TrialInd,strcmp(Parameters.GroupsName,P(TestInd).GroupsName)};
        %R.P(TestInd).SHUF=R.SHUF{TrialInd,strcmp(Parameters.GroupsName,P(TestInd).GroupsName)};
        
        
    end
end
directory='D:\data\Hyperscanning\BI-multiplayers\Groups\Results\';
save([directory 'results_GroupsChrono_' datestr(now,30) '.mat'],'R')

toc

figure
AUC=cellfun(@mean,[R.AUCall]);
hist(AUC(:),20)
figure
Perf=cellfun(@mean,[R.Perf]);
hist(Perf(:),20)

%% chronological analysis of performance during the BI session VERSUS OFFLINE
% this part compares the offline classification contained in R.Scores to
% the real online performance.
% you need to load the parameters as well as the results contained in R to
% proceed
close all
for indF=1:length(Flash)
    Serie{indF}=Flash2Serie(Flash{indF});%generate tag
    %figure;semilogy(Serie{indF})
end

close all
nbPinit=4
figure
for Geo=2%[1 2 3 4 5 6 7 8 9 10 11 14 16 17 18 19]
    nbP=nbPinit+10*(Geo-1);
Perf=(R.Scores{nbP}>0)==R.Ytest{nbP};
N=24;
CONVO{Geo}=conv(Perf+0,ones(1,N),'valid')./N;
Scale{Geo}=(1:length(CONVO{Geo}))./length(CONVO{Geo});
figure;
subplot(211);plot(CONVO{Geo});legend('Mean Perf Classification (sliding window 24 flashs)','Location','Best');title(num2str(Geo));
subplot(212);semilogy(Serie{Geo}(end-length(CONVO{Geo})+1:end));xlabel('Flash #');legend('Label (10^3=break) (10^2=new level) (10^1=new trial)','Location','Best')
end

%% chronological analysis of performance during the BI session BEHAVIOURAL ONLY (only need Flash.mat)
% this code just plot the real online performance.
% you need to load the FLASHES. The performance is computed thanks to the
% time signature of the game to differentiate trial/levels/breaks etc...
close all
for indS=1:size(Flash,2) % for each session
for indF=1:size(Flash,1)%for each group
    Serie{indF,indS}=Flash2Serie(Flash{indF,indS});%generate tag
    OnlineResults(indF,indS).fails=count(Serie{indF,indS}==10);
    OnlineResults(indF,indS).levels=count(Serie{indF,indS}==100)+count(Serie{indF,indS}==1000);
    OnlineResults(indF,indS).succes=count(Serie{indF,indS}==100)+count(Serie{indF,indS}==1000);
    OnlineResults(indF,indS).games=count(Serie{indF,indS}==1000);
    PerfG(indF,indS)=OnlineResults(indF,indS).succes/(OnlineResults(indF,indS).succes+OnlineResults(indF,indS).fails)
end
end
[OnlineResults(:,1).levels]'
[OnlineResults(:,3).fails]'+([OnlineResults(:,3).levels]'-36)
PerfSOLO=PerfG(:,[1:2]);
PerfSOLO=PerfSOLO(:);

close all
nbPinit=4
for Geo=2%[1 2 3 4 5 6 7 8 9 10 11 14 16 17 18 19]
figure;
semilogy(Serie{Geo});xlabel('Flash #');legend('Label (10^3=break) (10^2=new level) (10^1=new trial)','Location','Best')
end

    hfig=figure(1)
    subplot(1,4,1:4)
        xticklabels = 0:.1:1;
    FontSize=24
    set(hfig, 'PaperPosition', [0 0 30 20],'units','normalized','outerposition',[0.6 0.1 0.375 .9])
        set(gcf, 'color', [1 1 1])
        
        %show histogram solo versus multi
[n1, xout1] = hist(PerfG(:,3),20);
%bar(xout1,n1/19,.5); grid; hold all
h1=plot(xout1*100,smooth(n1/19),'m','linewidth',4); grid; hold all
x=[mean(PerfG(:,3))*100,mean(PerfG(:,3))*100]; %the mean performance of the multiplayer collaborative performance
y=[0,max(smooth(n1/19))*1.1];
plot(x,y,'m--','linewidth',4)
plot(x,y,'m')
text(mean(PerfG(:,3))*101,max(smooth(n1/19))*1.1,[num2str(mean(PerfG(:,3))*100,3) '%'],'fontsize',FontSize,'fontname','times new roman','FontAngle','italic','color','m') 
[n2, xout2] = hist(PerfSOLO,20);
%bar(xout2,n2/38,.5,'g');
h2=plot(xout2*100,smooth(n2/38),'b','linewidth',4);
x=[mean(PerfSOLO)*100,mean(PerfSOLO)*100];%the mean performance of the solo performance
y=[0,max(smooth(n1/19))*1.1];
plot(x,y,'b--','linewidth',4)
plot(x,y,'b')
text(mean(PerfSOLO)*101,max(smooth(n1/19))*1.1,[num2str(mean(PerfSOLO)*100,3) '%'],'fontsize',FontSize,'fontname','times new roman','FontAngle','italic','color','b') 
%bar(xout1(1:4),n1(1:4)/19,.5); grid; hold all
hold off;


xlabel('Performance (%)','fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
ylabel('Probability','fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
title('Online Classification','fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
    set(gca, 'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
    legend([h1 h2],{'MDM-hyper','MDM-solo'},'location','best','FontSize',FontSize-4)
        %set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
        set(gca,'YtickLabel',[])
        saveas(hfig,['D:\Mes Documents GIPSA\MATLAB\figures\BI-multiplayer\Online\online_MDMhyper_2players_perf2.tiff'])
saveas(hfig,['D:\Mes Documents GIPSA\MATLAB\figures\BI-multiplayer\Online\online_MDMhyper_2players_perf2.fig'])

%% performances analysis
AUCall=cellfun(@(x) x,[R.AUCall(11:50,:)]);
mAUCall(1,:)=mean(AUCall,1)
vAUCall(1,:)=var(AUCall,1)

Perf=cellfun(@(x) x,[R.Perf(11:50,:)]);
mPerf(1,:)=mean(Perf,1)
vPerf(1,:)=var(Perf,1)

figure
subplot(211)
plot(mAUCall)
hold on
plot(mAUCall+sqrt(vAUCall),'--','Color',[0.5 0.5 0.5])
plot(mAUCall-sqrt(vAUCall),'--','Color',[0.5 0.5 0.5])
hold off
subplot(212)
plot(mPerf)
hold on
plot(mPerf+sqrt(vPerf),'--','Color',[0.5 0.5 0.5])
plot(mPerf-sqrt(vPerf),'--','Color',[0.5 0.5 0.5])
hold off



%% plot analysis
nbsections=20;
figure
autoSubplot(2,hist(R.Scores{1},nbsections),hist(R.Scores{2},nbsections));
figure
autoSubplot(2,hist(R.Distances{1}(:,1),nbsections),hist(R.Distances{1}(:,2),nbsections),hist(R.Distances{2}(:,1),nbsections),hist(R.Distances{2}(:,2),nbsections));
figure
autoSubplot(2,hist(R.EIG{INDICEALL},nbsections),hist(R.EIG{INDICEALL}(:,2),nbsections));

%% EIGvalues analysis

for EIGind=1:length(R.Parameters.EIG)
    EIGindex=ParametersIndex(P,'EIG',R.Parameters.EIG{EIGind});
    INTRA=ParametersIndex(P,'Stats',R.Parameters.Stats{2});
    INTER=ParametersIndex(P,'Stats',R.Parameters.Stats{1});
    GROUPS=~ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{12})|ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{18});
    GROUP=ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{1});
    
    PerfINTRA(EIGind)=mean(mean([R.AUCall{:,EIGindex&INTRA&GROUPS}]));
    PerfINTER(EIGind)=mean(mean([R.AUCall{:,EIGindex&INTER&GROUPS}]));
    
end

figure;plot(PerfINTER);hold all;plot(PerfINTRA);legend(Parameters.Stats)
set(gca,'XTick',1:17,'XTickLabel',cellfun(@num2str,Parameters.EIG,'UniformOutput',0))
rotateticklabel(gca)
ylabel('AUC ROC')

%% Performances analysis (shuffle VS not_shuffle vs inter-intra vs test set ratio) ALL TEST
close all
%f1=figure
%title('lol')
%subplot(121)
load(['D:\data\Hyperscanning\BI-multiplayers\Groups\Results\' 'results_GroupsChrono_20140929T152419.mat'],'R')

INTRA=ParametersIndex(P,'Stats',R.Parameters.Stats{2});
INTER=ParametersIndex(P,'Stats',R.Parameters.Stats{1});


%NOP1=ParametersIndex(P,'P300_ref_orientation',R.Parameters.P300_ref_orientation{2});
%MULTIP1=ParametersIndex(P,'P300_ref_orientation',R.Parameters.P300_ref_orientation{1});
GROUPS=ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{1});%
%GROUPS=~(ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{12})|ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{18}));%remove group 12 and 18

R87=ParametersIndex(P,'TestSetRatio',0.8750);
R75=ParametersIndex(P,'TestSetRatio',0.75);
R62=ParametersIndex(P,'TestSetRatio',0.625);
R50=ParametersIndex(P,'TestSetRatio',0.50);
R25=ParametersIndex(P,'TestSetRatio',0.25);

PerfINTRA=[R.AUCall{:,INTRA&GROUPS}];
PerfINTER=[R.AUCall{:,INTER&GROUPS}];


ALLr={...
    [R.AUCall{:,INTRA&GROUPS&R87}]...
    [R.AUCall{:,INTRA&GROUPS&R75}]...
    [R.AUCall{:,INTRA&GROUPS&R62}]...
    [R.AUCall{:,INTRA&GROUPS&R50}]...
    [R.AUCall{:,INTRA&GROUPS&R25}];...
    [R.AUCall{:,INTER&GROUPS&R87}]...
    [R.AUCall{:,INTER&GROUPS&R75}]...
    [R.AUCall{:,INTER&GROUPS&R62}]...
    [R.AUCall{:,INTER&GROUPS&R50}]...
    [R.AUCall{:,INTER&GROUPS&R25}];...
    }
%    figure
% [yaa1 xaa1]=hist(ALLr{1,1},10,'r')
% bar(xaa1,yaa1)
% hold on
% [yaa1 xaa1]=hist(ALLr{4,1},10,'FaceColor',[0,0.7,0.7])
% bar(xaa1,yaa1)
% hold off
%figure
%boxplot([PerfINTRA' PerfINTER'])
%testSTATS=PermTest(PerfINTRA, PerfINTER)

%figure
%boxplot([PerfSHUFF' PerfnoSHUFF'])
%testSHUFF=PermTest(PerfSHUFF, PerfnoSHUFF)

%test=PermTest(ALLr{1,5},ALLr{1,4},10000)


%% P300 analysis
close all
INTER=ParametersIndex(P,'Stats',R.Parameters.Stats{1});
R25=ParametersIndex(P,'TestSetRatio',0.25);
GROUPS=ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{1});%
noSHUFF=ParametersIndex(P,'ShuffleCouple',R.Parameters.ShuffleCouple{2});
figure
subplot(221);plot(P1(:,:,1)');title('X');subplot(223);plot(P1(:,:,2)');title('Y')
P1=R.P1{:,INTER&R25&noSHUFF&GROUPS}
X=P1(:,:,1);
Y=P1(:,:,2);
[Ux Wx Vx]=svd(squeeze(X));
[Uy Wy Vy]=svd(squeeze(Y));
%subplot(232);plot((Ux*Wx*Vx')');subplot(235);plot((Uy*Wy*Vy')')


Y_x=pinv((sqrt(pinv(Wx))*Ux'))*(sqrt(pinv(Wy))*Uy')*P1(:,:,2)*(Vy*sqrt(pinv(Wy)))*pinv((Vx*sqrt(pinv(Wx))))
subplot(222);plot(P1(:,:,1)');title('X');subplot(224);plot(Y_x');
error=sum(sum((X-Y_x).^2))
title(['Y_x, error=' num2str(error)])
%% Performances analysis (shuffle VS not_shuffle vs inter-intra vs test set ratio) ITERATIVE TEST
clear SHtest STATtest
for indGroup=1:19
    INTRA=ParametersIndex(P,'Stats',R.Parameters.Stats{2});
    INTER=ParametersIndex(P,'Stats',R.Parameters.Stats{1});
    
    SHUFF=ParametersIndex(P,'ShuffleCouple',R.Parameters.ShuffleCouple{1});
    noSHUFF=ParametersIndex(P,'ShuffleCouple',R.Parameters.ShuffleCouple{2});
    
    %NOP1=ParametersIndex(P,'P300_ref_orientation',R.Parameters.P300_ref_orientation{2});
    %MULTIP1=ParametersIndex(P,'P300_ref_orientation',R.Parameters.P300_ref_orientation{1});
    GROUPS=ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{indGroup});%
    %GROUPS=~(ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{12})|ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{18}));%remove group 12 and 18
    
    R87=ParametersIndex(P,'TestSetRatio',0.8750);
    R75=ParametersIndex(P,'TestSetRatio',0.75);
    R62=ParametersIndex(P,'TestSetRatio',0.625);
    R50=ParametersIndex(P,'TestSetRatio',0.50);
    R25=ParametersIndex(P,'TestSetRatio',0.25);
    
    PerfINTRA=[R.AUCall{:,INTRA&GROUPS}];
    PerfINTER=[R.AUCall{:,INTER&GROUPS}];
    
    PerfSHUFF=[R.AUCall{:,GROUPS&SHUFF}];
    PerfnoSHUFF=[R.AUCall{:,GROUPS&noSHUFF}];
    
    ALLr={...
        [R.AUCall{:,INTRA&GROUPS&R87&SHUFF}]...
        [R.AUCall{:,INTRA&GROUPS&R75&SHUFF}]...
        [R.AUCall{:,INTRA&GROUPS&R62&SHUFF}]...
        [R.AUCall{:,INTRA&GROUPS&R50&SHUFF}]...
        [R.AUCall{:,INTRA&GROUPS&R25&SHUFF}];...
        [R.AUCall{:,INTER&GROUPS&R87&SHUFF}]...
        [R.AUCall{:,INTER&GROUPS&R75&SHUFF}]...
        [R.AUCall{:,INTER&GROUPS&R62&SHUFF}]...
        [R.AUCall{:,INTER&GROUPS&R50&SHUFF}]...
        [R.AUCall{:,INTER&GROUPS&R25&SHUFF}];...
        [R.AUCall{:,INTRA&GROUPS&R87&noSHUFF}]...
        [R.AUCall{:,INTRA&GROUPS&R75&noSHUFF}]...
        [R.AUCall{:,INTRA&GROUPS&R62&noSHUFF}]...
        [R.AUCall{:,INTRA&GROUPS&R50&noSHUFF}]...
        [R.AUCall{:,INTRA&GROUPS&R25&noSHUFF}];...
        [R.AUCall{:,INTER&GROUPS&R87&noSHUFF}]...
        [R.AUCall{:,INTER&GROUPS&R75&noSHUFF}]...
        [R.AUCall{:,INTER&GROUPS&R62&noSHUFF}]...
        [R.AUCall{:,INTER&GROUPS&R50&noSHUFF}]...
        [R.AUCall{:,INTER&GROUPS&R25&noSHUFF}];...
        ...
        };
    R.Parameters.GroupsName{indGroup}
    %verify Schuffle for each group
    nbPermu=1000
    for indRatio=1:length(R.Parameters.TestSetRatio)
        indRatio
        SHtest{indRatio}(indGroup,1)=PermTest(ALLr{3,indRatio},ALLr{1,indRatio},nbPermu);
        SHtest{indRatio}(indGroup,2)=PermTest(ALLr{4,indRatio},ALLr{2,indRatio},nbPermu);
        SHtest{indRatio}(indGroup,3)=mean(ALLr{3,indRatio});
        SHtest{indRatio}(indGroup,4)=mean(ALLr{1,indRatio});
        SHtest{indRatio}(indGroup,5)=mean(ALLr{4,indRatio});
        SHtest{indRatio}(indGroup,6)=mean(ALLr{2,indRatio});
        
        STATtest{indRatio}(indGroup,1)=PermTest(ALLr{2,indRatio},ALLr{1,indRatio},nbPermu);
        STATtest{indRatio}(indGroup,2)=PermTest(ALLr{4,indRatio},ALLr{3,indRatio},nbPermu);
        STATtest{indRatio}(indGroup,3)=mean(ALLr{2,indRatio});
        STATtest{indRatio}(indGroup,4)=mean(ALLr{1,indRatio});
        STATtest{indRatio}(indGroup,5)=mean(ALLr{4,indRatio});
        STATtest{indRatio}(indGroup,6)=mean(ALLr{3,indRatio});
        
    end
end
%% Performances analysis interative for all groups
close all
%f1=figure
%title('lol')
%subplot(121)

INTRA=ParametersIndex(P,'Stats',R.Parameters.Stats{2});
INTER=ParametersIndex(P,'Stats',R.Parameters.Stats{1});

SHUFF=ParametersIndex(P,'ShuffleCouple',R.Parameters.ShuffleCouple{1});
noSHUFF=ParametersIndex(P,'ShuffleCouple',R.Parameters.ShuffleCouple{2});

%NOP1=ParametersIndex(P,'P300_ref_orientation',R.Parameters.P300_ref_orientation{2});
%MULTIP1=ParametersIndex(P,'P300_ref_orientation',R.Parameters.P300_ref_orientation{1});
GROUPS=ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{1});%
%GROUPS=~(ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{12})|ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{18}));%remove group 12 and 18

R87=ParametersIndex(P,'TestSetRatio',0.8750);
R75=ParametersIndex(P,'TestSetRatio',0.75);
R62=ParametersIndex(P,'TestSetRatio',0.625);
R50=ParametersIndex(P,'TestSetRatio',0.50);
R25=ParametersIndex(P,'TestSetRatio',0.25);

PerfINTRA=[R.AUCall{:,INTRA&GROUPS}];
PerfINTER=[R.AUCall{:,INTER&GROUPS}];

PerfSHUFF=[R.AUCall{:,GROUPS&SHUFF}];
PerfnoSHUFF=[R.AUCall{:,GROUPS&noSHUFF}];

ALLr={...
    [R.AUCall{:,INTRA&GROUPS&R87&SHUFF}]...
    [R.AUCall{:,INTRA&GROUPS&R75&SHUFF}]...
    
    %%
    disp(GroupInd)
    set(gca,'XTick', 1:4);set(gca,'XTickLabel', {'P1' 'P2' 'INTRA' 'INTRA+INTER'})
    ylabel('ROC AUC score')
    subplot(122)
    
    
    
    boxplot([PerfSOLO1(:) PerfSOLO2(:) PerfINTRA' PerfINTER'])
    hold all
    disp(GroupInd)
    set(gca,'XTick', 1:4);set(gca,'XTickLabel', {'P1' 'P2' 'INTRA' 'INTRA+INTER'})
    saveas(f1,[directory 'figures\' 'ROC' 'all' '.jpeg'])
    
    %% Performances Analysis (comparison P_300 orientation)
    close all
    f1=figure
    title('lol')
    subplot(121)
    %for GroupInd=[1 2 3 4 5 6 7 8 9 10 11 13 14 15 16 17 19]
    
    SOLO=ParametersIndex(P,'Stats',R.Parameters.Stats{3});
    INTRA=ParametersIndex(P,'Stats',R.Parameters.Stats{2});
    INTER=ParametersIndex(P,'Stats',R.Parameters.Stats{1});
    
    NOP1=ParametersIndex(P,'P300_ref_orientation',R.Parameters.P300_ref_orientation{2});
    MULTIP1=ParametersIndex(P,'P300_ref_orientation',R.Parameters.P300_ref_orientation{1});
    
    GROUPS=~(ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{12})|ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{18}));
    
    PerfINTRA=[R.AUCall{:,INTRA&~MULTIP1}];
    PerfINTER=[R.AUCall{:,INTER&~MULTIP1}];
    PerfSOLO=[R.AUCall(:,SOLO&~MULTIP1)];
    PerfSOLO1=cellfun(@(x) x{1},PerfSOLO);
    PerfSOLO2=cellfun(@(x) x{2},PerfSOLO);
    
    
    boxplot([PerfSOLO1(:) PerfSOLO2(:) PerfINTRA' PerfINTER'])
    disp(GroupInd)
    set(gca,'XTick', 1:4);set(gca,'XTickLabel', {'P1' 'P2' 'INTRA' 'INTRA+INTER'})
    ylabel('ROC AUC score')
    subplot(122)
    %for GroupInd=[1 2 3 4 5 6 7 8 9 10 11 13 14 15 16 17 19]
    
    SOLO=ParametersIndex(P,'Stats',R.Parameters.Stats{3});
    INTRA=ParametersIndex(P,'Stats',R.Parameters.Stats{2});
    INTER=ParametersIndex(P,'Stats',R.Parameters.Stats{1});
    
    NOP1=ParametersIndex(P,'P300_ref_orientation',R.Parameters.P300_ref_orientation{2});
    MULTIP1=ParametersIndex(P,'P300_ref_orientation',R.Parameters.P300_ref_orientation{1});
    
    GROUPS=~(ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{12})|ParametersIndex(P,'GroupsName',R.Parameters.GroupsName{18}));
    
    PerfINTRA=[R.AUCall{:,INTRA&MULTIP1}];
    PerfINTER=[R.AUCall{:,INTER&MULTIP1}];
    PerfSOLO=[R.AUCall(:,SOLO&MULTIP1)];
    PerfSOLO1=cellfun(@(x) x{1},PerfSOLO);
    PerfSOLO2=cellfun(@(x) x{2},PerfSOLO);
    
    
    boxplot([PerfSOLO1(:) PerfSOLO2(:) PerfINTRA' PerfINTER'])
    hold all
    disp(GroupInd)
    set(gca,'XTick', 1:4);set(gca,'XTickLabel', {'P1' 'P2' 'INTRA' 'INTRA+INTER'})
    saveas(f1,[directory 'figures\' 'ROC' 'all' '.jpeg'])
    
    %end
    %%
    mean(mean(xcorr(R.P1{1,1}(:,:,1),R.P1{2,1}(:,:,1))))
    %plot(R.P1{:,1}(:,:,1)')
    %%
    figure;plot(PerfINTER);hold all;plot(PerfINTRA);legend(Parameters.Stats)
    set(gca,'XTick',1:17,'XTickLabel',cellfun(@num2str,Parameters.EIG,'UniformOutput',0))
    rotateticklabel(gca)
    ylabel('AUC ROC')
    
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
    
    selectedGroups=[1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,19];
    %%
    selectedGroups=[1];
    
    Groups=ALLgroups(selectedGroups);
    AUCn=cell2mat(AUC(selectedGroups,:,SelectedMap,1:2,:));
    AUCn=squeeze(AUCn)
    subplot(211);hist(scores{1,strcmp(Stats,'intra'),3,Trial,1}(Ytest==1));hold all;scores{1,strcmp(Stats,'intra'),3,Trial,1}(Ytest==2)
    scores{1,strcmp(Stats,'all'),3,Trial,1}
    scores{1,strcmp(Stats,'intra'),3,Trial,2}
    scores{1,strcmp(Stats,'all'),3,Trial,2}
    
    %% Analyse mapping
    close all
    MEANmapping=squeeze(mean(mean(AUCn,2),1));
    VARmapping=squeeze(mean(var(AUCn,[],4),1));
    figure
    subplot(211);plot(MEANmapping');legend(Stats);ylabel('Mean ROC Perf');xlabel('Set of electrods')
    [MAXperf INDmap]=max(MEANmapping,[],2);
    hold on;plot(INDmap,MEANmapping(:,INDmap),'or')
    title('Test on different set of electrods')
    
    subplot(212);plot(VARmapping');legend(Stats);ylabel('Var ROC Perf');xlabel('Set of electrods')
    [MAXperf INDmap]=min(VARmapping,[],2);
    hold on;plot(INDmap,VARmapping(:,INDmap),'or')
    
    size(AUCn)
    %%
    SelectedMap=3
    figure
    
    plot(MEANmapping(:,SelectedMap)');legend(Stats);ylabel('Mean ROC Perf');xlabel('Set of electrods')
    [MAXperf INDmap]=max(MEANmapping,[],2);
    hold on;plot(INDmap,MEANmapping(SelectedMap)-VARmapping(SelectedMap,INDmap),'or')
    title('Test on different set of electrods')
    plot(VARmapping');legend(Stats);ylabel('Var ROC Perf');xlabel('Set of electrods')
    [MAXperf INDmap]=min(VARmapping,[],2);
    hold on;plot(INDmap,VARmapping(SelectedMap,INDmap),'or')
    
    %% Analyse intra VS inter
    close all
    SelectedMap=3
    AUCmean=squeeze(mean(AUCn(:,:,SelectedMap,:),4));
    figure;bar(squeeze(AUCmean(:,1)),1,'g')
    axis([0 19 min(min(AUCmean)) 1])
    hold all;bar(squeeze(AUCmean(:,2)),0.5,'r');hold off;
    set(gca, 'XTick', 1:length(Groups),'XTickLabel', Groups);
    legend(Stats)
    ylabel('AUC ROC')
    title({'Overview of group performances' ...
    'If Green>Red, the InterStats increase perf'...
    'If Red>Green, the InterStats decrease perf'})
    %%
    close all
    SelectedMap=3
    AUCmean=squeeze(mean(AUCn(:,:,SelectedMap,:),4));
    figure;bar(squeeze(AUCmean(:,1)),1,'g')
    axis([0 19 min(AUCmean) 1])
    hold all;bar(squeeze(AUCmean(:,2)),0.5,'r');hold off;
    set(gca, 'XTick', 1:length(Groups),'XTickLabel', Groups);
    legend(Stats)
    ylabel('AUC ROC')
    title({'Overview of group performances' ...
    'If Green>Red, the InterStats increase perf'...
    'If Red>Green, the InterStats decrease perf'})
    %%
    close all
    SelectedMap=3
    AUCmean=squeeze(mean(AUCn(:,:,SelectedMap,:),4));
    figure;subplot(311);bar(squeeze(AUCmean(:,1)),1,'g')
    axis([0 19 0.5 1])
    hold all;bar(squeeze(AUCmean(:,2)),0.5,'r');hold off;
    set(gca, 'XTick', 1:length(Groups),'XTickLabel', Groups);
    legend(Stats)
    ylabel('AUC ROC')
    title({'Overview of group performances' ...
    'If Green>Red, the InterStats increase perf'...
    'If Red>Green, the InterStats decrease perf'})
    
    diffAUC=AUCmean(:,1)-AUCmean(:,2);
    x=squeeze(mean(SYNC.riem(selectedGroups,SelectedMap,:),3));y=diffAUC;
    subplot(312); plot(x,y,'*')
    title(['Correlation Coefficient=' num2str(corrcoef(x,y))])
    xlabel('RiemSYNC, d_R^2(C^{-1/2}*C_{ref}*C^{-1/2} <-> \itI)');ylabel('Diff Perf AUCinter-AUCintra')
    x=squeeze(mean(SYNC.FroNormP1P2(selectedGroups,SelectedMap,:),3));y=diffAUC;
    subplot(313); plot(x,y,'*')
    title(['Correlation Coefficient=' num2str(corrcoef(x,y))])
    xlabel('FroNorm, ||(C^{-1/2}*C_{ref}*C^{-1/2} - \itI||_F^2');ylabel('Diff Perf AUCinter-AUCintra')
    
    %% permutation test
    close all
    %--- set some parameters
    alpha   = 0.05;     % significance level
    nPerm   = 10000;    % number of permutations (i.e., size of surrogate data)
    nObs    = length(ALLr{1,1});     	% number of observations
    mean_1  = 10;       % mean of first vector of observations
    mean_2  = 12.5;       % mean of second vector of observations
    var_1   = 5;        % variance of first vector of observations
    var_2   = 5;        % variance of second vector of observations
    
    
    % --- simulate two vectors of paired obervation
    vec_OBS_1   = mean_1 + var_1*randn(1,nObs);
    vec_OBS_2   = mean_2 + var_2*randn(1,nObs);
    
    vec_OBS_1   = ALLr{3,5};
    vec_OBS_2   = ALLr{4,5};
    
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
        % compute surrogate surrogate value for this specific permutationLrL
        diff_PERM(perm_ix) =  mean(vec_PERM_2 - vec_PERM_1);
    end
    
    % compute statistical significance and plot results
    p_val = (1/nPerm)*(1+length(find(diff_PERM > diff_OBS)))
    h_val = p_val < alpha;
    %figure;dispBootstrap(diff_OBS, diff_PERM)
