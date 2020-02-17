%% COMPUTE
clear all
load('MARTI_GroupsName.mat');users=GroupsName;
% load('EKATE_GroupsName.mat');users=ALLgroups;
% users=Generate_Users_Numbers([1:71]);
Directory='D:\data\Hyperscanning\MARTI\Groups'
nsession=4
% dbstop in load_EEG_data at 30 if indU==4

for indU=1:length(users)
    %     try
    %bi2015a and earlier
    data = load_EEG_data(Directory,['' users{indU}]);
    for indS=1:nsession
        %             try
        StimCode=data.session(indS).h.EVENT.TYP;
        nrep2(indU,indS)=count(StimCode==33040 | StimCode==1409);
        nrep(indU,indS)=count(StimCode==33279);
        nlvl(indU,indS)=count(StimCode==786);
        nlvl2(indU,indS)=count(StimCode==781);
        nNT(indU,indS)=count(StimCode==33286);
        nTA(indU,indS)=count(StimCode==33285); 
        nstart(indU,indS)=count(StimCode==100);
        nstop(indU,indS)=count(StimCode==101);
        
        %TARGET SUCCESS ?
        StimCode==33285
                BEGlvl=find(StimCode==786) 
                ENDlvl=find(StimCode==781)
                
                

        % ONLY bi2015b BEGIN (65th channel is stimulation code)
        if size(data.session(indS).Channels,2)==65,
            StimCode2=findpeaks(data.session(indS).Channels(:,65));
            adaptive(indU,indS)=triggers_BI2015b_extractONLINE(StimCode2);
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

load([outputfolder 'results.mat'],'bi2015b')
bi2015b.adaptive=adaptive;
save([outputfolder 'results.mat'],'bi2015b','-append')
%% ANALYSIS
clear all;clc;close all
% inputfolder='D:\Mes Documents GIPSA\EEG Hyperscanning\Hyperscanning, common documents\2016 JOURNAL BCI\data\';
inputfolder='D:\GoogleDrive\Présentations & Communications\2016-05 BCI Meeting 2016\data\'
outputfolder=['D:\GoogleDrive\Présentations & Communications\2016-05 BCI Meeting 2016\data\bi2015\figures\']

load([inputfolder 'results.mat'],'bi2015b')
FontSize=12;
Scaling=20;
PaperFactor=1
results={bi2015b.adaptive.ACC_table};
ALLsessions=cat(3,results{:});
adaptive=mean([ALLsessions ALLsessions(:,1,:)&ALLsessions(:,2,:)],3);
ACC=mean(mean(mean(ALLsessions(:,1:2,:))));ACC=repmat(ACC,size(adaptive,1),1)
figure
plot([adaptive(:,1:3),adaptive(:,5)],'linewidth',1.5);hold on;
plot([ACC],'--k','linewidth',0.75);hold off;

% xlabel('Nb Trials','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
xlabel('Nombre de niveaux','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
% ylabel('Rate','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
ylabel('Taux','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
set(gcf, 'PaperPosition', [0 0 12.5 6]*PaperFactor,'units','normalized','outerposition',[0.4 0.1 0.5 .9])
% legend('p1','p2','failure','p1 & p2','location','northeastoutside')
legend('p1','p2','échec','p1 & p2','location','northeastoutside')
legend boxoff                   % Hides the legend's axes (legend border and background)
ylim([0 1])

print([outputfolder '2015b_adaptation_full'], '-dpng', '-r450');
print([outputfolder '2015b_adaptation_full'], '-dpdf', '-r450');

figure
plot([adaptive(:,1:2)],'linewidth',1.5);hold on;
plot([ACC],'--k','linewidth',0.75);hold off;
xlabel('Nb Trials','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
ylabel('Rate','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
set(gcf, 'PaperPosition', [0 0 12.5 6]*PaperFactor,'units','normalized','outerposition',[0.4 0.1 0.5 .9])
ylim([0 1])
legend('p1','p2','location','northeastoutside')
legend boxoff                   % Hides the legend's axes (legend border and background)
print([outputfolder 'adaptation'], '--', '-r450');print([outputfolder 'adaptation'], '-dpdf', '-r450');

%%
Step=5;
meanAdap=[];GRperf=[];close all
for indG=1:22
    for indS=1:2
GRresults={bi2015b.adaptive(indG,:).ACC_table}
GRsessions=cat(3,GRresults{:});
GRperf=mean([GRsessions(:,1:2,:) GRsessions(:,1,:)&GRsessions(:,2,:)],3);
test=GRperf(:,indS);
meanAdap(:,indS+2*(indG-1))=[mean(test(1:Step)) mean(test(1+Step:Step*2)) mean(test(Step*2+1:Step*3)) mean(test(end-Step+1:end))];
    end
end
figure
plot(meanAdap)
count(diff(meanAdap)<0.0001)
