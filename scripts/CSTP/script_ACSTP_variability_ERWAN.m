% details: use CSTPinfo
%
% *** History: 2016-02-11
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
%
% see also : ACSTP, CSTP, applyCSTP, EnsembleAverage, meanOverlap, epoch_p300, WeightsEstimation,
% CorrectLatency, ConvergenceLatencies, best_Pz, CSTPinfo
%
% This script is made to load a EEG dataset and compute the ACSTP from
% randomly selected trials in order to compute the variability. The
% variability is computed from RMSE
%
% it also compute and same the AEA and the EA ACSTP
%% (0) load early parameters

clc
close all
clear all
            Directory= ['D:\data\erwan\mat\']  ; %change path if needed

%users={'A01','A02','A03','A04','A05','A06','A07','A08'};
nbStims=[400 80]
NbStim2Take={round(nbStims*0.5) round(nbStims*0.25) round(nbStims*0.2)...
    round(nbStims*0.15) round(nbStims*0.1) round(nbStims*0.05)}

%NbStim2Take={[40 40] [20 20] [10 10]...
%   [5 5] [20 20] round(nbStims*0.01)}

NbTrials=500;

save([Directory 'parameters.mat'], 'nbStims','NbStim2Take')
%% (1) MAIN LOOP prepare data, load parameters, COMPUTE ACSTP and save RESULT

% ds
%
clear FILTER
tic
for indUser=7:24%:length(users) %select the subject
    for indSession=1%:8
        try
            close all
            clear EEG
            nameSave='_variability_';
            Directory= ['D:\data\erwan\mat\']  ; %change path if needed
            % format : ERWAN_SS1_s4_online_non-adaptive.mat
            load([Directory 'ERWAN_SS' num2str(indUser) '_s' num2str(indSession) '_training_non-adaptive.mat'])
            doFilter=1;
            % EEG is a structure with
            %              Fs: scalar (sample rate in Hz)
            %         Trigger: [nb samples x1 ] Trigger channel of '0' with '1' at the start
            %                   of each sweep. There are [nb epochs] '1'.
            %      EpochClass: [nb epochs x1] class of the sweeps (0 for Non-TARGET, 1
            %                   for TARGET).
            %        Channels: [nb samples x nb channels] preprocessed EEG recordings
            %   NoiseTrigger*: the equivalent of Trigger but for the sweep of the noise.
            %                   By default, it takes the same.
            % ElectrodesName*: {1 x nb channels} the names of the electrodes (usefull
            %                  only in case of plot, i.e. ACSTPoptions.DISPLAY=true)
            Artefacts=[];
            
            if doFilter
                %%%%%%%%%%%%%%%%%% PREPROCESSING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                decimationfactor=4; %(5) put 1 to do nothing
                f1=1; %(8) low cutoff freq  (bandpass)
                f2=20; %(9) high cutoff freq  (bandpass)
                N=4; %(10) filter order (bandpass)
                
                %%%%%%%%%%%%%%%%%%%%%% PREPROCESSING (no user input needed)%%%%%%%%%%%%%%%%
                [EEG.Channels EEG.Fs EEG.Trigger]=preprocessingEEG(double(EEG.Channels),EEG.Fs,[f1 f2 N decimationfactor],EEG.Trigger);
            end
            %%%%%%%%%%%%%%%%%% CSTP PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Window=round(EEG.Fs); %(11) the sweeps window will be 1s
            Delays=10; %(12) +/- nb shifted samples allowed for the jitter correction
            offset=0;%round(-0.500*EEG.Fs); %(13) to add an offset to all the epochs (default: 0)
            if offset~=0
                EEG.Trigger=Position2Trigger(find(EEG.Trigger)+offset,size(EEG.Channels,1));
            end
            %user parameters for the ACSTP to improve converge :
            winElec=[7,9,10,11,12,13,14,15,16]; %(13)* electrodes used for latency calculation and Pz selection
            % exemple: Cz(7),P7(9),P3(10),Pz(11),P4(12),P8(13),O1(14),Oz(15),O2(16)
            %winElec=findElectrodes(EEG.ElectrodesName,EEG.ElectrodesName);%{'Fz','Cz','Pz'})
            winTime=[floor((0.05)*EEG.Fs):ceil((0.550)*EEG.Fs)]-offset; %(14)* time window (in sample) used for latency calculation and Pz selection
            % exemple: 50ms to 550ms
            %   *optional
            
            % WARNING: FOLLOWING LINE ONLY TO COPY THE .MAT FILES TO .TXT FILES
            %FULLDIR=['D:\Mes Documents GIPSA\MATLAB\LK_TOOLBOX\data\ARNAUD\txt\'  num2str(indUser) '\' ]
            %EEG_mat2txt(FULLDIR,EEG)
            
            EEG;
            ACSTPoptions.SubspaceDim=size(EEG.Channels,2):-1:size(EEG.Channels,2)-6;
            ACSTPoptions.SubspaceDim=16:-1:12;
            ACSTPoptions.Epoch_size=Window;
            ACSTPoptions.LatencyCorr_max=3;%Delays;
            ACSTPoptions.Mask_Electrodes=winElec;
            ACSTPoptions.Mask_Time=winTime;
            ACSTPoptions.DISPLAY=0;
            ACSTPoptions.computeClassLat=[1]; %Compute only for both Class
            %ACSTPoptions.Weights=1
            %% ACSTP loop Algorithm
            EEGsave=EEG;
            for indBootstrap=1:length(NbStim2Take)
                disp(['USER' num2str(indUser) ' session' num2str(indSession) ' boot'  num2str(indBootstrap) ])
                
                clear FILTER
                outputFILE=[Directory 'results\' '2ERWAN_' nameSave '_s' num2str(indUser) '_se' num2str(indSession) '_boot'  num2str(indBootstrap) '.mat'];
                if ~exist(outputFILE) %check if 
                    
                    for indRND=1:NbTrials %do it 100 times
                        RND=[];
                        [EEG RNDseed{indRND}]=schuffle_continuousEEG(EEGsave,NbStim2Take{indBootstrap},RND);
                        [~, FILTER(indRND)]=ACSTP(EEG,ACSTPoptions);
                        
                        nbB=length(NbStim2Take);
                        inB=indBootstrap;
                        disp(['Subject ' num2str(indUser) ' Variability computation at: ' ...
                            num2str(100*(((inB-1)*NbTrials+(indRND))/(NbTrials*nbB))) ' %'])
                    end
                    INFOS.RNDseed=RNDseed;
                    INFOS.NumberStimPerClass=NbStim2Take{indBootstrap};
                    INFOS.NbTrials=NbTrials;
                    
                    save(outputFILE,'FILTER','INFOS')
                    disp(['file saved'])
                else
                    disp(['file already exists, skip...'])
                end
            end
        catch e
            disp('file not found')
            %e
        end
    end
end
toc

%% (2) Butterfly plot with the variability of each estimation
% this script plot the epochs estimations as a butterfly plot with the
% corresponding mean.

clear all
nbStims=[400 80]
pourStim=[0.25,0.10,.07,.05,.03,.01];
NbStim2Take={round(nbStims*0.5) round(nbStims*0.25) round(nbStims*0.2)...
    round(nbStims*0.15) round(nbStims*0.1) round(nbStims*0.05)};
Classes={'NT','TA'};
Fs=128;
Directory= ['D:\data\erwan\mat\']  ; %change path if needed
DirectoryFig='D:\Mes Documents GIPSA\MATLAB\figures\CSTP\ERWAN\'
load('BImapping')
Channels2Plot=[4,7,11];
FontSize=24;
Scaling=20;
PaperFactor=0.8

for indUser=3:12
    for indSession=1:1
        for indBootstrap=1:6
            try
                disp(['s' num2str(indUser) '_se' num2str(indSession) '_boot'  num2str(indBootstrap)])
                %referenceS=load(['D:\data\BNCI\ALS\results\' 'BNCI_ALS_.mat']);
                load([Directory 'results\' 'ERWAN__variability__s' num2str(indUser) '_se' num2str(indSession) '_boot'  num2str(indBootstrap) '.mat'],'FILTER','INFOS');
                
                %prepare parameters for plotting
                offset=0;
                xticklabels = (offset/Fs):.25:(Fs+offset)/Fs;
                xticks = linspace(0, (1.01)*(Fs)/Fs, numel(xticklabels));
                
                tmpcell1={FILTER(:).EA};
                tmpmat1=cat(4,tmpcell1{:});
                tmpcell2={FILTER(:).EAcstp};
                tmpmat2=cat(4,tmpcell2{:});
                
                if 0 %remove common ref
                    commonRef1=mean(tmpmat1,1)
                    commonRef1=repmat(commonRef1,[size(tmpmat1,1) 1 1 1]);
                    commonRef2=mean(tmpmat2,1)
                    commonRef2=repmat(commonRef2,[size(tmpmat2,1) 1 1 1]);
                    tmpmat1=tmpmat1-commonRef1;
                    tmpmat2=tmpmat2-commonRef2;
                end
                
                if 0 %baseline correction
                    baselineRef1=mean(tmpmat1(:,1:5,:,:),2)
                    baselineRef1=repmat(baselineRef1,[1 size(tmpmat1,2) 1 1]);
                    baselineRef2=mean(tmpmat1(:,1:5,:,:),2)
                    baselineRef2=repmat(baselineRef2,[1 size(tmpmat2,2) 1 1]);
                    tmpmat1=tmpmat1-baselineRef1;
                    tmpmat2=tmpmat2-baselineRef2;
                end
                
                for indClass=1:2 %main log for ploting
%                      clf(gcf)
                    compt=1;
                    % prepare y ticks labels
                    yticklabels = -Scaling/2:Scaling:Scaling/2;
                    yticks = linspace(-(1.01)*(Scaling), (1.01)*(Scaling), numel(yticklabels));
%                     figure
%                     set(gcf, 'Visible', 'off');
                    subplot1(length(Channels2Plot),2,'Gap',[0.01 0])
                    
                    for indCh=Channels2Plot
                        
                        
                        %plot first column
                        subplot1((compt-1)*2+1)
                        if compt==1 title('AEA','FontSize',FontSize,'fontname','times new roman','Fontangle','italic');end
                        
                        [h1 h2 h3]=plotEEGvariability(squeeze(tmpmat1(:,:,indClass,:)),'channel',indCh,'Fs',Fs);
                        axis([0 1 -Scaling Scaling]);
                        ylabel(LABELS(indCh),'fontsize',FontSize,'fontname','times new roman','Fontangle','italic');
                        set(gca,'Ytick',yticklabels,'Yticklabel',yticklabels,'FontSize',FontSize*.75,'fontname','times new roman','Fontangle','italic');
                        
                        if compt==length(Channels2Plot) % if last line, put the xlabel
                            xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
                            set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
                        end
                        
                        set(get(gca,'YLabel'),'Rotation',0) %put the ylabel horizontal
                        % Move of 15% the ylabel text on the left to avoid
                        % collision
                        h = get(gca,'ylabel');
                        pos = get(h,'position');
                        xlimits = get(gca,'xlim');
                        pos(1) = xlimits(1) -0.15 * (xlimits(2) - xlimits(1));
                        set(h,'position',pos)
                        
                        %plot second column
                        subplot1((compt-1)*2+2)
                        if compt==1 title('ACSTP','FontSize',FontSize,'fontname','times new roman','Fontangle','italic');end
                        
                        [h1 h2 h3]=plotEEGvariability(squeeze(tmpmat2(:,:,indClass,:)),'channel',indCh,'Fs',Fs);
                        axis([0 1 -Scaling Scaling])
                        if compt==length(Channels2Plot) %if last line, put the xlabel
                            xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
                            set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize*.75,'fontname','times new roman','FontAngle','italic')
                        end
                        compt=compt+1;
                    end
                    % prepare xlabel ticks
                    
                    %grid on
                    
                    %                     remining plot with the mean over electrodes
                    if 0 %plot the mean over electrodes
                        subplot1((compt-1)*2+1)
                        [h1 h2 h3]=plotEEGvariability(squeeze(tmpmat1(:,:,indClass,:)),'channel',1:16,'Fs',128);
                        axis([0 1 -Scaling Scaling])
                        ylabel('mean','fontsize',FontSize+4,'fontname','times new roman','Fontangle','italic')
                        xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
                        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
                        set(gca,'Ytick',yticks,'Yticklabel',yticklabels,'FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
                        set(get(gca,'YLabel'),'Rotation',0)
                        subplot1((compt-1)*2+2)
                        [h1 h2 h3]=plotEEGvariability(squeeze(tmpmat2(:,:,indClass,:)),'channel',1:16,'Fs',128);
                        axis([0 1 -Scaling Scaling])
                        xlabel('time (s)','FontSize',FontSize,'fontname','times new roman','Fontangle','italic')
                        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
                    end
                    set(gcf, 'PaperPosition', [0 0 20 30]*PaperFactor,'units','normalized','outerposition',[0.4 0.1 0.5 .9])
                    annotation('textbox', [0.25,0.02,0.5,.02],'units','normalized',...
                        'String',['Subj. ' cell2mat(Generate_Users_Numbers(indUser))...
                        ' - Sess. ' num2str(indSession) ' - ' Classes{indClass}...
                        ' - K='  num2str(NbStim2Take{indBootstrap}(indClass)) ''],...
                        'fontsize',(FontSize+4)*PaperFactor,'fontname','times new roman','Fontangle','italic',...
                        'EdgeColor','none');
                    print(gcf, [DirectoryFig 'butterfly\' num2str(PaperFactor*100) '_butterfly_subj' cell2mat(Generate_Users_Numbers(indUser)) '_sess' num2str(indSession) '_boot'  num2str(indBootstrap) '_' Classes{indClass}  ''],'-dtiff','-r450')
                    disp('file written')
                end 
                
            catch e
                disp(['Error: ' e.message])
                
            end
        end
    end
end
disp('Butterfly plots done')
%% (3) save the RMSE statistics (and Pz)
clear all
indSession=1
Directory= ['D:\data\erwan\mat\']  ; %change path if needed
load('BImapping')
for indUser=1:24
    indSubject=indUser;%:length(FILTER)
    for indBootstrap=1:6
        %referenceS=load(['D:\data\BNCI\ALS\results\' 'BNCI_ALS_.mat']);
        load([Directory 'results\' 'ERWAN__variability__s' num2str(indUser) '_se' num2str(indSession) '_boot'  num2str(indBootstrap) '.mat'],'FILTER','INFOS');
        
        %referenceS.FILTER(indSubject).EA;
        
        tmpcell={FILTER(:).EA};
        tmpmat=cat(4,tmpcell{:});
        Fboot(indUser,indBootstrap).EA_mean=mean(tmpmat,4);
        Fboot(indUser,indBootstrap).EA_var=std(tmpmat,[],4);
        Fboot(indUser,indBootstrap).EA_rmse=RMSE(tmpmat);
        
        
        tmpcell={FILTER(:).EAcstp};
        tmpmat=cat(4,tmpcell{:});
        Fboot(indUser,indBootstrap).EAcstp_mean=mean(tmpmat,4);
        Fboot(indUser,indBootstrap).EAcstp_var=std(tmpmat,[],4);
        Fboot(indUser,indBootstrap).EAcstp_rmse=RMSE(tmpmat);
        
        tmpcell={FILTER(:).BestPz};
        tmpmat=cat(1,tmpcell{:});
        Pz.TA(indUser,indBootstrap,:)=tmpmat(:,2);
                Pz.NT(indUser,indBootstrap,:)=tmpmat(:,1);

    end
end
save([Directory 'results\ERWAN_variability_RMSE.mat'],'Fboot','Pz')
disp('DONE')


%% (4) LOAD results + visualization RMSE + save pval in txt
close all
clear all
Directory= ['D:\data\erwan\mat\']  ; %change path if needed
DirectoryFig='D:\Mes Documents GIPSA\MATLAB\figures\CSTP\ERWAN\';
load([Directory 'results\ERWAN_variability_RMSE.mat'],'Fboot');
stdfactor=1/3;
MinMax=[0 7 0 10];
FontSize=12;
PaperFactor=0.33
load([Directory 'parameters.mat'])%figure
for indUser=1:24
    for indBootstrap=1:6
        %mean
        tmpTA1=Fboot(indUser,indBootstrap).EA_rmse(:,:,2);
        tmpNT1=Fboot(indUser,indBootstrap).EA_rmse(:,:,1);
        tmpTA2=Fboot(indUser,indBootstrap).EAcstp_rmse(:,:,2);
        tmpNT2=Fboot(indUser,indBootstrap).EAcstp_rmse(:,:,1);
        
        TA(indUser,indBootstrap)=mean(mean(tmpTA1,2),1);
        NT(indUser,indBootstrap)=mean(mean(tmpNT1,2),1);
        TAcstp(indUser,indBootstrap)=mean(mean(tmpTA2,2),1);
        NTcstp(indUser,indBootstrap)=mean(mean(tmpNT2,2),1);
        %std
        TAstd(indUser,indBootstrap)=std( reshape(tmpTA1,[],1))*stdfactor;
        NTstd(indUser,indBootstrap)=std( reshape(tmpNT1,[],1))*stdfactor;
        TAcstpstd(indUser,indBootstrap)=std( reshape(tmpTA2,[],1))*stdfactor;
        NTcstpstd(indUser,indBootstrap)=std( reshape(tmpNT2,[],1))*stdfactor;
        
        [~,TApval(indUser,indBootstrap)]=ttest(tmpTA1(:),tmpTA2(:),0.05);
        [~,NTpval(indUser,indBootstrap)]=ttest(tmpNT1(:),tmpNT2(:),0.05);
%         TApval(TApval>0.05)=nan;
%         NTpval(NTpval>0.05)=nan;
    end
end
% for indUser=1:24
dlmwrite(['D:\Mes Documents GIPSA\MATLAB\figures\CSTP\ERWAN\rmse\pval.txt'],[TApval(:)' NTpval(:)'],'delimiter','\n', 'precision', '%.24f')
% dlmwrite(['D:\Mes Documents GIPSA\MATLAB\figures\CSTP\ERWAN\rmse\pvalNTss' num2str(indUser) '.txt'],NTpval(indUser,:)','delimiter','\n', 'precision', '%.24f')
% end
%% (4b) PRINT RMSE visualization
disp('DONE')

close all
LegendHist={'AEA','ACSTP'};
XLegend={'A','B','C','D','E','F'};
XLegend1={'','40','20','16','12','8','4'}';
XLegend2={'','200','100','80','60','40','20'}';
%xticks = linspace(0,7, numel(XLegend1));

for indGroup=1:3
    close all
    hrmse=figure;
    subplot1(8,2,'Gap',[0 0.005]);
    
    
    ind=1;
    for indUser=(1:8)+(indGroup-1)*8
        clear errY
        subplot1(1+(ind-1)*2)
        y=[TA(indUser,:)' TAcstp(indUser,:)'];
        %bar([NT(1,:)' NTcstp(1,:)'])
        errY(:,:,1) = [TAstd(indUser,:)' TAcstpstd(indUser,:)'];   % 10% lower error
        remplacementVal=(MinMax(end)-errY*stdfactor)*0.5;
        y((errY*stdfactor+y)>=MinMax(end)*.8)=remplacementVal((errY*stdfactor+y)>=MinMax(end)*0.8);
        hbar=barwitherr(errY, y);    % Plot with errorbars
        axis(MinMax)
        sigstar({[0.8,1.2], [0.8,1.2]+1,[0.8,1.2]+2,[0.8,1.2]+3,[0.8,1.2]+4,[0.8,1.2]+5},TApval(indUser,:))
        set(hbar(1),'FaceColor',[1 1 1]*0.2);
        set(hbar(2),'FaceColor',[1 1 1]*0.8);
        ylabel(['s' cell2mat(Generate_Users_Numbers(indUser))],'FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
        set(get(gca,'YLabel'),'Rotation',0) %put the ylabel horizontal
        % Move of 15% the ylabel text on the left to avoid
        % collision
        h = get(gca,'ylabel');
        pos = get(h,'position');
        xlimits = get(gca,'xlim');
        ylimits = get(gca,'ylim');
        pos(1) = xlimits(1) -0.15 * (xlimits(2) - xlimits(1));
        pos(2) = ylimits(1) +0.35 * (ylimits(2) - ylimits(1));
        set(h,'position',pos)
        if ind==1
            title(['TA'],'FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
        end
        if ind==8
            set(gca,'xtick',0:(length(XLegend1)-1),'xticklabel',XLegend1,'FontSize',FontSize*.75,'fontname','times new roman','FontAngle','italic');
            set(gca,'ytick',[0 MinMax(end)],'yticklabel',[0 MinMax(end)])
             rotateXLabels(gca,45)
        else
            set(gca,'yticklabel',[])
        end
        subplot1(2+(ind-1)*2);
        y=[NT(indUser,:)' NTcstp(indUser,:)'];
        %bar([NT(1,:)' NTcstp(1,:)'])
        errY(:,:,1) = [NTstd(indUser,:)' NTcstpstd(indUser,:)'];   % 10% lower error
        y(y>=MinMax(end)*0.9)=MinMax(end)*0.01;
        
        hbar=barwitherr(errY, y);    % Plot with errorbars
        sigstar({[0.8,1.2], [0.8,1.2]+1,[0.8,1.2]+2,[0.8,1.2]+3,[0.8,1.2]+4,[0.8,1.2]+5},NTpval(indUser,:))
        axis(MinMax);
        set(hbar(1),'FaceColor',[1 1 1]*0.2);
        set(hbar(2),'FaceColor',[1 1 1]*0.8);
        if ind==1
            title(['NT'],'FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
        end
        if ind==8
            set(gca,'xtick',0:(length(XLegend2)-1),'xticklabel',XLegend2,'FontSize',FontSize*.75,'fontname','times new roman','FontAngle','italic');
             rotateXLabels(gca,45)

        end
        ind=ind+1;
    end
    set(hrmse, 'PaperPosition', [0 0 25 60]*PaperFactor,'units','normalized','outerposition',[0.45 0.1 0.5 .9])
    
    set(gcf, 'color', [1 1 1])
    % colormap('gray')
    annotation('textbox', [0.05,0.001,0.5,.05],'units','normalized',...
        'String',['RMSE - Black(AEA) White(ACSTP) - std/' num2str(1/stdfactor)],'fontsize',(FontSize+10)*PaperFactor,'fontname','arial new roman','EdgeColor','none','fontangle','italic');
    annotation('textbox', [0.05,0.001,0.5,.025],'units','normalized',...
        'String',['*p-val<0.04761138'],'fontsize',(FontSize+10)*PaperFactor,'fontname','arial new roman','EdgeColor','none','fontangle','italic');
  
    print(hrmse,[DirectoryFig '\rmse\RMSE_grey' num2str(indGroup) '_' num2str(PaperFactor*100)],'-dtiff','-r450')
end


%% (5) LOAD results + visualization Distribution Pz
close all
clear all
Directory= ['D:\data\erwan\mat\']  ; %change path if needed
DirectoryFig='D:\Mes Documents GIPSA\MATLAB\figures\CSPT\ERWAN\'
load([Directory 'results\ERWAN_variability_RMSE.mat'],'Pz')

hrmse=figure
FontSize=12
MinMax=[0 7 13 16]
LegendHist={'AEA','ACSTP'}
XLegend={'A','B','C','D','E','F'}
XLegend1={'40','20','16','12','8','4'}'
XLegend2={'200','100','80','60','40','20'}'
%xticks = linspace(0,7, numel(XLegend1));
Pzmean=mean(Pz.TA,3);
Pzst=std(Pz.TA,[],3);
ind=1

clear errY
y=mean([Pzmean(:,:)'],2)
%bar([NT(1,:)' NTcstp(1,:)'])
errY(:,:,1) = mean([Pzst(:,:)'],2);   % 10% lower error
hbar=barwitherr(errY, y);    % Plot with errorbars
axis(MinMax)
set(hbar(1),'FaceColor',[1 1 1]*0.6);

ylabel(['Subspace dimension (P_z)'],'FontSize',FontSize,'fontname','times new roman','FontAngle','italic');

title(['TA'],'FontSize',FontSize,'fontname','times new roman','FontAngle','italic');

set(gca,'xticklabel',XLegend1,'FontSize',FontSize/2,'fontname','times new roman','FontAngle','italic');


set(hrmse, 'PaperPosition', [0 0 10 8],'units','normalized','outerposition',[0.4 0.1 0.5 .9])

set(gcf, 'color', [1 1 1])
colormap('gray')
xlabel(['Nb trials'],'fontsize',FontSize,'fontname','times new roman','EdgeColor','none');
print(hrmse,[DirectoryFig '\Pz\TA'],'-dtiff','-r450')

%% (6) BUILD Table for ANOVA on RMSE versus TA/NT and correlation with Boot
close all
clear all
Directory= ['D:\data\erwan\mat\']  ; %change path if needed
DirectoryFig='D:\Mes Documents GIPSA\MATLAB\figures\CSPT\ERWAN\';
load([Directory 'results\ERWAN_variability_RMSE.mat'],'Fboot');
load([Directory 'parameters.mat'])
stdfactor=1;
MinMax=[0 7 0 14.99];
TYPE={'AEA','ACSTP'};
CLASS={'NT','TA'};
% TABLE
% ----------------
% Subject | TYPE | CLASS | Boot | RMSE
% 01      |  AEA |  TA   |  40  | val
% ...     | CSTP |  NT   |  ... | ...
% 24      |  ... | ...   |  4   | ...

%figure
for indUser=1:24
    for indBootstrap=1:6
        %mean
        tmpTA1=Fboot(indUser,indBootstrap).EA_rmse(:,:,2);
        tmpNT1=Fboot(indUser,indBootstrap).EA_rmse(:,:,1);
        tmpTA2=Fboot(indUser,indBootstrap).EAcstp_rmse(:,:,2);
        tmpNT2=Fboot(indUser,indBootstrap).EAcstp_rmse(:,:,1);
        
        TA(indUser,indBootstrap)=mean(mean(tmpTA1,2),1);
        NT(indUser,indBootstrap)=mean(mean(tmpNT1,2),1);
        TAcstp(indUser,indBootstrap)=mean(mean(tmpTA2,2),1);
        NTcstp(indUser,indBootstrap)=mean(mean(tmpNT2,2),1);
        %std
        TAstd(indUser,indBootstrap)=std( reshape(tmpTA1,[],1))*stdfactor;
        NTstd(indUser,indBootstrap)=std( reshape(tmpNT1,[],1))*stdfactor;
        TAcstpstd(indUser,indBootstrap)=std( reshape(tmpTA2,[],1))*stdfactor;
        NTcstpstd(indUser,indBootstrap)=std( reshape(tmpNT2,[],1))*stdfactor;
        
         [~,TApval(indUser,indBootstrap)]=ttest(tmpTA1(:),tmpTA2(:),0.05);
         [~,NTpval(indUser,indBootstrap)]=ttest(tmpNT1(:),tmpNT2(:),0.05);
         TApval(TApval>0.05)=nan;
         NTpval(NTpval>0.05)=nan;
    end
end
tmp1=cat(3,NT,TA);
tmp2=cat(3,NTcstp,TAcstp);
allMean=cat(4,tmp1,tmp2);
tmp1=cat(3,NTstd,TAstd);
tmp2=cat(3,NTcstpstd,TAcstpstd);
allStd=cat(4,tmp1,tmp2);
disp('computation DONE')
disp('Formatting table...')
%
indLine=1
for indUser=1:24
    for indBootstrap=1:6
        for indType=1:2
            for indClass=1:2
        %mean
        TABLE{indLine, 1}=indUser; %# subject
        TABLE{indLine, 2}=TYPE{indType}; %TYPE
        TABLE{indLine, 3}=CLASS{indClass}; %Class
        TABLE{indLine,4}=NbStim2Take{indBootstrap}(indClass);
        TABLE{indLine,5}=allMean(indUser,indBootstrap,indClass,indType);
        TABLE{indLine,6}=allStd(indUser,indBootstrap,indClass,indType);
        indLine=indLine+1;
    end
        end
    end
end
cell2csv('test.csv', TABLE,';')
disp('DONE')