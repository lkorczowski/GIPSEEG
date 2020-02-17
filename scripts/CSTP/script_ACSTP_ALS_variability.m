% details: use CSTPinfo
%
% *** History: 2015-03-19
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

users={'A01','A02','A03','A04','A05','A06','A07','A08'};
nbStims=[3500 700]
pourStim=[0.25,0.10,.07,.05,.03,.01];
NbStim2Take={round(nbStims*0.25) round(nbStims*0.10) round(nbStims*0.07)...
    round(nbStims*0.05) round(nbStims*0.03) round(nbStims*0.01)}
NbTrials=100;
%%%%%%%%%%%%%%%%%% INPUT FILE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% there are data from 4 subjects avawilable in this exemple


% Save (cell array with the following items) :
% (1) AEA
% (2) AEA after weights, latency corrections and CSTP
% (3) CSTP{iter} Bs Bt As At : the temporal and spatials filters for each
% Pz
% (4) Pzind : optimal indice of iter such as Pz is the optimal subspace
%       with respect of the mask
%
% Thus Save{2}{Save{4}) gives P300 with optimal estimation.

%% (1) MAIN LOOP prepare data, load parameters, COMPUTE ACSTP and save RESULT

% ds
%
clear FILTER
tic
for indUser=1:length(users) %select the subject
    toc
    disp(['USER' num2str(indUser)])
    close all
    clear EEG
    nameSave='ALS_variability_'
    Directory= ['D:\data\BNCI\ALS\' users{indUser} '.mat' ] ; %change path if needed
    load( [Directory])
    EEG.Fs=256;
    [EEG.Trigger EEG.EpochClass]= triggers_continuous2discrete(data.y);
    EEG.EpochClass=EEG.EpochClass-1;
    EEG.Channels=data.X;
    EEG.ElectrodesName=data.channels;
    clear data
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
        decimationfactor=2; %(5) put 1 to do nothing
        NOTCH=1; %(6) put 1 to remove 50Hz
        BANDPASS=1; %(7) put 1 to BANDPASS (see the filter below for the cutoff freq)
        f1=0.5; %(8) low cutoff freq  (bandpass)
        f2=40; %(9) high cutoff freq  (bandpass)
        N=0; %(10) filter order (bandpass)
        
        %%%%%%%%%%%%%%%%%%%%%% PREPROCESSING (no user input needed)%%%%%%%%%%%%%%%%
        [EEG.Channels EEG.Fs EEG.Trigger]=preprocessingEEG(double(EEG.Channels),EEG.Fs,[f1 f2 N decimationfactor],EEG.Trigger);
    end
    %%%%%%%%%%%%%%%%%% CSTP PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Window=round(EEG.Fs); %(11) the sweeps window will be 1s
    Delays=10; %(12) +/- nb shifted samples allowed for the jitter correction
    offset=0%round(-0.500*EEG.Fs); %(13) to add an offset to all the epochs (default: 0)
    if offset~=0
        EEG.Trigger=Position2Trigger(find(EEG.Trigger)+offset,size(EEG.Channels,1));
    end
    %user parameters for the ACSTP to improve converge :
    %winElec=[7,9,10,11,12,13,14,15,16]; %(13)* electrodes used for latency calculation and Pz selection
    % exemple: Cz(7),P7(9),P3(10),Pz(11),P4(12),P8(13),O1(14),Oz(15),O2(16)
    winElec=findElectrodes(EEG.ElectrodesName,EEG.ElectrodesName);%{'Fz','Cz','Pz'})
    winTime=[floor((0.05)*EEG.Fs):ceil((0.750)*EEG.Fs)]-offset; %(14)* time window (in sample) used for latency calculation and Pz selection
    % exemple: 50ms to 550ms
    %   *optional
    
    % WARNING: FOLLOWING LINE ONLY TO COPY THE .MAT FILES TO .TXT FILES
    %FULLDIR=['D:\Mes Documents GIPSA\MATLAB\LK_TOOLBOX\data\ARNAUD\txt\'  num2str(indUser) '\' ]
    %EEG_mat2txt(FULLDIR,EEG)
    
    EEG
    ACSTPoptions.SubspaceDim=size(EEG.Channels,2):-1:size(EEG.Channels,2)-6;
    ACSTPoptions.SubspaceDim=7;
    ACSTPoptions.Epoch_size=Window;
    ACSTPoptions.LatencyCorr_max=0;%Delays;
    ACSTPoptions.Mask_Electrodes=winElec;
    ACSTPoptions.Mask_Time=winTime;
    ACSTPoptions.DISPLAY=0;
    ACSTPoptions.computeClassLat=[0 1] %Compute only for both Class
    %ACSTPoptions.Weights=1
    %% ACSTP loop Algorithm
    EEGsave=EEG;
    for indBootstrap=1:length(NbStim2Take)
        for indRND=1:NbTrials %do it 100 times
            RND=[];
            [EEG RNDseed{indRND}]=schuffle_continuousEEG(EEGsave,NbStim2Take{indBootstrap},RND);
            [Epochs{indRND} FILTER(indRND)]=ACSTP(EEG,ACSTPoptions);
            
            nbB=length(NbStim2Take);
            inB=indBootstrap;
            disp(['Subject ' num2str(indUser) ' Variability computation at: ' ...
                num2str(100*(((inB-1)*NbTrials+(indRND))/(NbTrials*nbB))) ' %'])
        end
        INFOS.RNDseed=RNDseed;
        INFOS.NumberStimPerClass=NbStim2Take{indBootstrap};
        INFOS.NbTrials=NbTrials;
        save(['D:\data\BNCI\ALS\results\' 'BNCI_' nameSave '_s' num2str(indUser) '_boot'  num2str(indBootstrap) '.mat'],'Epochs','FILTER','INFOS')
    end
end
toc

%% (2) save the RMSE statistics
for indUser=1:8
    indSubject=indUser%:length(FILTER)
    for indBootstrap=1:6
        referenceS=load(['D:\data\BNCI\ALS\results\' 'BNCI_ALS_.mat']);
        load(['D:\data\BNCI\ALS\results\' 'BNCI_ALS_variability' '__s' num2str(indSubject) '_boot'  num2str(indBootstrap) '.mat'],'Epochs','FILTER','INFOS');
        Directory= ['D:\data\BNCI\ALS\' users{indSubject} '.mat' ] ; %change path if needed
        load(Directory)
        EEG.ElectrodesName=data.channels;
        EEG.Epoch_Size=128
        clear data
        nameSave=['Variability_boot' num2str(indBootstrap)];
        referenceS.FILTER(indSubject).EA
        
        tmpcell={FILTER(:).EA}
        tmpmat=cat(4,tmpcell{:});
        Fboot(indUser,indBootstrap).EA_mean=mean(tmpmat,4);
        Fboot(indUser,indBootstrap).EA_var=std(tmpmat,[],4)
        Fboot(indUser,indBootstrap).EA_rmse=RMSE(tmpmat)
        
        
        tmpcell={FILTER(:).EAcstp}
        tmpmat=cat(4,tmpcell{:});
        Fboot(indUser,indBootstrap).EAcstp_mean=mean(tmpmat,4);
        Fboot(indUser,indBootstrap).EAcstp_var=std(tmpmat,[],4)
        Fboot(indUser,indBootstrap).EAcstp_rmse=RMSE(tmpmat)
        
    end
end
save(['D:\data\BNCI\ALS\results\BNCI_variability_RMSE2.mat'],'Fboot')
%% (3) visualization of AEA vs CSTP for each class indepentely
for indUser=1:8
    indSubject=indUser%:length(FILTER)
    for indBootstrap=1:6
        referenceS=load(['D:\data\BNCI\ALS\results\' 'BNCI_ALS_.mat']);
        load(['D:\data\BNCI\ALS\results\' 'BNCI_ALS_variability' '__s' num2str(indSubject) '_boot'  num2str(indBootstrap) '.mat'],'Epochs','FILTER','INFOS');
        Directory= ['D:\data\BNCI\ALS\' users{indSubject} '.mat' ] ; %change path if needed
        load(Directory)
        EEG.ElectrodesName=data.channels;
        EEG.Epoch_Size=128
        clear data
        nameSave=['Variability_boot' num2str(indBootstrap)];
        referenceS.FILTER(indSubject).EA
        
        
        load(['D:\data\BNCI\ALS\results\BNCI_variability_RMSE.mat'],'Fboot') %load RMSE statistics

%%
Scale=8 %init scale

for classIND=1:2
            hfig=figure(classIND)
    for IND=1:2
        nbLINE=9
        nbCOLUMN=5
        FontSize=24;
        
        Mapping=(EEG.ElectrodesName');
        
        t = (0:(ACSTPoptions.Epoch_size-1))./EEG.Fs;
        PZs=ACSTPoptions.SubspaceDim;
        winElec2=winElec;
        winTime2=winTime;
        [Pzind]=find(FILTER(indSubject).BestPz(classIND)==PZs)
        
        switch IND
            case 1
                %PLOTX=mean(X(:,:,EpochClass==1),3);
                %PLOTX=FILTER(indSubject).EA(:,:,classIND);
                PLOTX=Fboot.EA_mean(:,:,classIND);
                PLOTXstd=Fboot.EA_var(:,:,classIND);
                info=['AEA'];
                if FILTER(indSubject).Class(classIND)==0
                    info=[info ' NT'];
                else
                    info=[info ' TA']
                end
                %Scale=FindScaling(PLOTX,winElec,winTime);
            case 2
                %PLOTX=Xcstp(:,:,EpochClass==1);
                %PLOTX=FILTER(indSubject).EAcstp(:,:,classIND);
                PLOTX=Fboot.EAcstp_mean(:,:,classIND);
                PLOTXstd=Fboot.EAcstp_var(:,:,classIND);
                info=['ACSTP(' num2str(7) ')'];
                
        end
        xticklabels = (offset/EEG.Fs):.5:(ACSTPoptions.Epoch_size+offset)/EEG.Fs;
        
        clims=[0 7]
        ClassTag=1
        set(hfig, 'PaperPosition', [0 0 25 60],'units','normalized','outerposition',[0.6 0.1 0.375 .9])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT EEG %%%%%%%%%%%%%%%
        subplot(nbLINE,nbCOLUMN,[1 2 6 7 11 12 16 17 21 22 26 27 31 32 36 37]+2*(IND-1));plotEEG(PLOTX,Scale,EEG.Fs,Mapping) % Xbar TA
        %set(gca,'yticklabel',[])
        if IND>1
            set(gca, 'color', [0.95 .95 .95])
        end
        set(gcf, 'color', [1 1 1])
        xlim([0 1-offset/EEG.Fs])
        
        title(info,'FontSize',FontSize)
        %xlabel('Time (s)')
        xticks = linspace(0, (1.01)*(ACSTPoptions.Epoch_size)/EEG.Fs, numel(xticklabels));
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
        set(gcf, 'InvertHardCopy', 'off');
        % global field power on average
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT GFP mean %%%%%%%%%%%%%%%
        subplot(nbLINE,nbCOLUMN,[21 22]+20+2*(IND-1))
        area(t,global_field_power(PLOTX'),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
        if IND==1
            %ylabel('GFP','FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
            SCALEgfp=max(global_field_power(PLOTX'));
        end
        xlim([0 1-offset/EEG.Fs])
        ylim([0 SCALEgfp]);grid on;
        if IND==1
            %ylabel('GFP','FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
        end
        %set(gca,'YtickLabel',[])
        %xticks = linspace(0, (1.01)*(ACSTPoptions.Epoch_size)/EEG.Fs, numel(xticklabels));
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
        set(gca, 'color', [0.95 .95 .95])
        
        
    end
    spaceplots
    
    saveas(hfig,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\ALS\EAvsCSTP\' nameSave '_Subject_' num2str(indSubject) 'class' num2str(classIND) '.tiff'])
    
end
    
    %% TA/NT + differences TA-NT

    close all
    indSubject=indUser%:length(FILTER)
    ClassName={'NT','TA'}
    hfigdiff=figure(1)
    nbLINE=9
    nbCOLUMN=8
    FontSize=24;
    set(hfigdiff, 'PaperPosition', [0 0 25*nbCOLUMN/5 60*nbLINE/9],'units','normalized','outerposition',[0.1 0.1 0.8 .9])
    
    for IND=1:4
        
        for classIND=1:2
            
            
            
            Mapping=(EEG.ElectrodesName');
            
            t = (0:(ACSTPoptions.Epoch_size-1))./EEG.Fs;
            PZs=ACSTPoptions.SubspaceDim;
            winElec2=winElec;
            winTime2=winTime;
            [Pzind]=find(FILTER(indSubject).BestPz(classIND)==PZs)
            
            switch IND
                case 1
                    %PLOTX=mean(X(:,:,EpochClass==1),3);
                    PLOTX=FILTER(indSubject).EA(:,:,classIND);
                    info=['AEA'];
                    if FILTER(indSubject).Class(classIND)==0
                        info=[info ''];
                    else
                        info=[info '']
                    end
                    %Scale=FindScaling(PLOTX,winElec,winTime);
                case 2
                    %PLOTX=Xcstp(:,:,EpochClass==1);
                    PLOTX=FILTER(indSubject).EAcstp(:,:,classIND);
                    info=['ACSTP(' num2str(PZs(Pzind)) ')'];
                case 3
                    PLOTX=FILTER(indSubject).EA(:,:,2)-FILTER(indSubject).EA(:,:,1);
                    info=['AEA TA-NT'];
                case 4
                    PLOTX=FILTER(indSubject).EAcstp(:,:,2)-FILTER(indSubject).EAcstp(:,:,1);
                    info=['ACSTP TA-NT'];
                    
                    
            end
            xticklabels = (offset/EEG.Fs):.5:(ACSTPoptions.Epoch_size+offset)/EEG.Fs;
            
            clims=[0 7]
            ClassTag=1
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT EEG %%%%%%%%%%%%%%%
            TemporalCells=(1:nbCOLUMN:nbCOLUMN*(nbLINE-2));TemporalCells=sort([TemporalCells,TemporalCells+1])
            if classIND==1
                LineStype='--';
                LW=1;
                Color=[0.2,0.6,0.2];
            else
                
                LineStype='-';
                LW=2;
                Color=[0.3,0.3,0.3];
                
                
                
            end
            if IND>2
                                                Color=[0.1,0.1,0.1];
            end
            hold on;
            subplot(nbLINE,nbCOLUMN,TemporalCells+2*(IND-1));h(:,IND)=plotEEG(PLOTX,Scale,EEG.Fs,Mapping,LW,Color,'none',LineStype) % Xbar TA
            hold off;
            %set(gca,'yticklabel',[])
            if mod(IND-1,2)>0
                set(gca, 'color', [0.95 .95 .95])
            end
            set(gcf, 'color', [1 1 1])
            xlim([0 1-offset/EEG.Fs])
            
            title(info,'FontSize',FontSize)
            %xlabel('Time (s)')
            xticks = linspace(0, (1.01)*(ACSTPoptions.Epoch_size)/EEG.Fs, numel(xticklabels));
            set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
            set(gcf, 'InvertHardCopy', 'off');
            % global field power on average
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT GFP mean %%%%%%%%%%%%%%%
            if IND>2
            TemporalCells=(nbCOLUMN*(nbLINE-2)+1:nbCOLUMN:nbCOLUMN*(nbLINE-1+1));
            else
               TemporalCells=(nbCOLUMN*(nbLINE-classIND)+1:nbCOLUMN:nbCOLUMN*(nbLINE-classIND+1));
            end
            TemporalCells=sort([TemporalCells,TemporalCells+1])
            hold all;
            subplot(nbLINE,nbCOLUMN,TemporalCells+2*(IND-1))
            area(t,global_field_power(PLOTX'),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
            hold off;
            
            if IND==1
                %ylabel('GFP','FontSize',FontSize,'fontname','times new roman','FontAngle','italic')
                SCALEgfp=Scale%max(global_field_power(PLOTX'));
            end
            xlim([0 1-offset/EEG.Fs])
            ylim([0 SCALEgfp]);grid on;
            if IND==1
                ylabel(ClassName(classIND),'FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
            end
            %set(gca,'YtickLabel',[])
            %xticks = linspace(0, (1.01)*(ACSTPoptions.Epoch_size)/EEG.Fs, numel(xticklabels));
            set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
            set(gca, 'color', [0.95 .95 .95])
            
            
        end
        %spaceplots
        
        annotation('textbox', [0.3,0.01,0.5,.05],'units','normalized',...
           'String',['Subject ' num2str(indSubject) ' ; NT(' num2str(NbStim2Take{indBootstrap}(1)) ') TA(' num2str(NbStim2Take{indBootstrap}(2)) ') ; Top: EEG scaling: ' ...
           num2str(Scale*2) 'µV' ' ; Bottom: GFP (µV^2)'],'fontsize',FontSize+4,'fontname','times new roman');
        
    end
    saveas(hfigdiff,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\ALS\EAvsCSTP\EAestimation\Subject_' num2str(indSubject) '_boot_' num2str(indBootstrap)  '_all.tiff'])
    
end 
end

%% (4) LOAD results + visualization RMSE
close all
        referenceS=load(['D:\data\BNCI\ALS\results\' 'BNCI_ALS_.mat']);

load(['D:\data\BNCI\ALS\results\BNCI_variability_RMSE2.mat'],'Fboot')
%figure
for indUser=1:8
    for indBootstrap=1:6
%mean
TA(indUser,indBootstrap)=mean(mean(Fboot(indUser,indBootstrap).EA_rmse(:,:,2),2),1);
NT(indUser,indBootstrap)=mean(mean(Fboot(indUser,indBootstrap).EA_rmse(:,:,1),2),1);
TAcstp(indUser,indBootstrap)=mean(mean(Fboot(indUser,indBootstrap).EAcstp_rmse(:,:,2),2),1);
NTcstp(indUser,indBootstrap)=mean(mean(Fboot(indUser,indBootstrap).EAcstp_rmse(:,:,1),2),1);
%std
TAstd(indUser,indBootstrap)=std( reshape(Fboot(indUser,indBootstrap).EA_rmse(:,:,2),[],1));
NTstd(indUser,indBootstrap)=std( reshape(Fboot(indUser,indBootstrap).EA_rmse(:,:,1),[],1));
TAcstpstd(indUser,indBootstrap)=std( reshape(Fboot(indUser,indBootstrap).EAcstp_rmse(:,:,2),[],1));
NTcstpstd(indUser,indBootstrap)=std( reshape(Fboot(indUser,indBootstrap).EAcstp_rmse(:,:,1),[],1));
    end
end
%%
close all
hrmse=figure
FontSize=24
subplot1(8,2,'Gap',[0 0])
MinMax=[0 7 0 6.99]
LegendHist={'AEA','ACSTP'}
XLegend={'A','B','C','D','E','F'}
XLegend1={'','175','70','49','35','21','7'}'
XLegend2={'','875','350','245','175','105','35'}'
%xticks = linspace(0,7, numel(XLegend1));

for indUser=1:8
    clear errY
    subplot1(1+(indUser-1)*2)
    y=[TA(indUser,:)' TAcstp(indUser,:)']
    %bar([NT(1,:)' NTcstp(1,:)'])
  errY(:,:,1) = [TAstd(indUser,:)' TAcstpstd(indUser,:)'];   % 10% lower error
  barwitherr(errY, y);    % Plot with errorbars
        axis(MinMax)

    ylabel(['s' num2str(indUser)],'FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
    if indUser==1
        title(['TA'],'FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
    end
    if indUser==8
        set(gca,'xticklabel',XLegend1,'FontSize',FontSize/2,'fontname','times new roman','FontAngle','italic');
    end
    subplot1(2+(indUser-1)*2)
    y=[NT(indUser,:)' NTcstp(indUser,:)']
    %bar([NT(1,:)' NTcstp(1,:)'])
  errY(:,:,1) = [NTstd(indUser,:)' NTcstpstd(indUser,:)'];   % 10% lower error
  barwitherr(errY, y);    % Plot with errorbars

    axis(MinMax)
    if indUser==1
        title(['NT'],'FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
    end
    if indUser==8
        set(gca,'xticklabel',XLegend2,'FontSize',FontSize/2,'fontname','times new roman','FontAngle','italic');
    end
end
    set(hrmse, 'PaperPosition', [0 0 25 60],'units','normalized','outerposition',[0.4 0.1 0.5 .9])

            set(gcf, 'color', [1 1 1])
colormap('gray')
  annotation('textbox', [0.3,0.01,0.5,.05],'units','normalized',...
           'String',['RMSE: deviation from the average of the average at a given bootstrap. Black: EA, White: EAacstp(7)'],'fontsize',FontSize-4,'fontname','times new roman');
    saveas(hrmse,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\ALS\EAvsCSTP\EAestimation\2_RMSE_grey.tiff'])

%%
subplot(211)
hist(        sum(sum(Fboot(indUser,indBootstrap).EA_rmse(:,:,2),2),1));
subplot(212)
hist(        sum(sum(Fboot(indUser,indBootstrap).EAcstp_rmse(:,:,2),2),1));

sum(sum(referenceS.FILTER(indSubject).EA(:,:,2).^2))
sum(sum(referenceS.FILTER(indSubject).EAcstp(:,:,2).^2))

figure
subplot(211)
hist(   sum(referenceS.FILTER(indSubject).EA(:,:,2)-FILTER(1).EA(:,:,2)));
subplot(212)
hist(   sum(referenceS.FILTER(indSubject).EAcstp(:,:,2)-FILTER(1).EAcstp(:,:,2)));
