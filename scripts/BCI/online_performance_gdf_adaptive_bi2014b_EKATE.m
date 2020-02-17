%% COMPUTE
clear all
load('EKATE_GroupsName.mat');users=ALLgroups;
% load('EKATE_GroupsName.mat');users=ALLgroups;
% users=Generate_Users_Numbers([1:71]);
Directory='F:\data\Hyperscanning\EKATE\Groups'
nsession=3;
% dbstop in load_EEG_data at 30 if indU==4

for indU=1:length(users)
    %     try
    %bi2015a and earlier
    data = load_EEG_data(Directory,['' users{indU}]);
    for indS=1:nsession
        %             try
        h=data.session(indS).h;
        T=length(data.session(indS).Trigger);
        StimCode=h.EVENT.TYP;
        nrep2(indU,indS)=count(StimCode==33040 | StimCode==1409);
        nrep(indU,indS)=count(StimCode==33279);%33279 : end of rep
        nlvl(indU,indS)=count(StimCode==786);%786 : end of LEVEL
        nlvl2(indU,indS)=count(StimCode==781);%786 : give positive feedback
        nNT(indU,indS)=count(StimCode==33286);
        nTA(indU,indS)=count(StimCode==33285);
        nstart(indU,indS)=count(StimCode==100);
        nstop(indU,indS)=count(StimCode==101);
        
        
        % ONLY bi2015b BEGIN (65th channel is stimulation code)
        if size(data.session(indS).Channels,2)==65,
            StimCode2=findpeaks(data.session(indS).Channels(:,65));
            adaptive(indU,indS)=triggers_BI2015b_extractONLINE(StimCode2);
        else
            %Begin of Level
            lvlbegPos =  h.EVENT.POS(h.EVENT.TYP==781);
            adaptive(indU,indS).lvlbeg = zeros(T,1);
            adaptive(indU,indS).lvlbeg(lvlbegPos)=1;
            
            %End of Level
            lvlendPos =  h.EVENT.POS(h.EVENT.TYP==781);
            adaptive(indU,indS).lvlend = zeros(T,1);
            adaptive(indU,indS).lvlend(lvlendPos)=1;
            
            %End of Repetition
            endrep =  h.EVENT.POS(h.EVENT.TYP==33279);
            adaptive(indU,indS).endrep = zeros(T,1);
            adaptive(indU,indS).endrep(endrep) = 1;
            
            %Repetition online success (1 : successed, 0 : failed)
            adaptive(indU,indS).OnlineSuccessREP=zeros(length(endrep),1);
            for indEND=find(adaptive(indU,indS).lvlend)'
                tmpREP=indEND-endrep;
                adaptive(indU,indS).OnlineSuccessREP=adaptive(indU,indS).OnlineSuccessREP+(min(tmpREP(tmpREP>0))==tmpREP);
                %=TargetPos;
            end         
            adaptive(indU,indS).StimCode=StimCode;                       
        end
        %             catch
        %                 disp(['error user' num2str(indU) 'session' num2str(indS)])
        
        %             end
    end
    %     catch
    %         disp(['error user' num2str(indU)])
    %     end
end
%% SAVE
% outputfolder='D:\Mes Documents GIPSA\EEG Hyperscanning\Hyperscanning, common documents\2016 JOURNAL BCI\data\';
outputfolder='D:\GoogleDrive\Présentations & Communications\2016-05 BCI Meeting 2016\data\';
load([outputfolder 'results.mat'],'bi2014b')
bi2014b.adaptive=adaptive;bi2014b.nrep=nrep;bi2014b.nlvl=nlvl;bi2014b.nlvl2=nlvl2;
bi2014b.nTA=nTA;bi2014b.nNT=nNT;bi2014b.nrep2=nrep2;
save([outputfolder 'results.mat'],'bi2014b','-append')
%% ANALYSIS
clear all;clc;close all
inputfolder='D:\GoogleDrive\Présentations & Communications\2016-05 BCI Meeting 2016\data\';
outputfolder=[inputfolder '\bi2014\figures\']

load([inputfolder 'results.mat'],'bi2014b')
FontSize=12;
Scaling=20;
PaperFactor=1
maxAdaptation=max(cellfun(@length,{bi2014b.adaptive(:,1:2).OnlineSuccessREP}));
maxAdaptation2=max(cellfun(@length,{bi2014b.adaptive(:,3).OnlineSuccessREP}));
mean(cellfun(@mean,{bi2014b.adaptive.OnlineSuccessREP}))
ALLsessions1=NaN(size(bi2014b.adaptive,1),2,maxAdaptation);
ALLsessions2=NaN(size(bi2014b.adaptive,1),1,maxAdaptation2);

ACC=[];
for indGroup=1:size(bi2014b.adaptive,1)
    for indSession=1:size(bi2014b.adaptive,2)
        if indSession<3
        ALLsessions1(indGroup,indSession,1:length(bi2014b.adaptive(indGroup,indSession).OnlineSuccessREP))=bi2014b.adaptive(indGroup,indSession).OnlineSuccessREP;
        else
        ALLsessions2(indGroup,1,1:length(bi2014b.adaptive(indGroup,indSession).OnlineSuccessREP))=bi2014b.adaptive(indGroup,indSession).OnlineSuccessREP;
        end
        ACC(indGroup,indSession)=mean(bi2014b.adaptive(indGroup,indSession).OnlineSuccessREP);
    end
end
perfAdaptive1=squeeze(nanmean(ALLsessions1,1))';
perfAdaptive2=squeeze(nanmean(ALLsessions2,1));
ACC2=[];
ACC2(:,1)=mean(mean(ACC(:,1:2),1));ACC2(:,2)=mean(mean(ACC(:,3),1));
ACCplot=repmat(ACC2,maxAdaptation2,1);

figure
plot([smooth(mean(perfAdaptive1(1:20,:),2),'lowess')],'linewidth',1.5);hold all;
plot([smooth(perfAdaptive2(1:20,:),'lowess')],'linewidth',1.5);
plot([ACCplot],'--k','linewidth',0.75);hold off;
xlim([0 20])

xlabel('Nb Répétitions','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
ylabel('Taux Réussite','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
set(gcf, 'PaperPosition', [0 0 12.5 6]*PaperFactor,'units','normalized','outerposition',[0.4 0.1 0.5 .9])
legend('SOLO','COLLABORATIVE','location','northeastoutside')
legend boxoff                   % Hides the legend's axes (legend border and background)
ylim([0 1])

print([outputfolder 'adaptation_full'], '-dpng', '-r450');
print([outputfolder 'adaptation_full_FR'], '-dpdf', '-r450');