% details: use CSTPinfo
%
% *** History: 2015-03-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
%
% see also : ACSTP, CSTP, applyCSTP, EnsembleAverage, meanOverlap, epoch_p300, WeightsEstimation,
% CorrectLatency, ConvergenceLatencies, best_Pz, CSTPinfo
clc
close all
clear all

users={'cba','clm','ega','fsa','gro','mba','mma','mta','pla','sce','sph','wpa'};
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

%% prepare data and load parameters
%
%indice of the subject 'data' (if several)
%%prepare data
clear FILTER
%load data
for indUser=11%:length(users) %select the subject
    close all
    nameSave='sparseGONOGO'
    Directory= ['D:\data\ARNAUD\' users{indUser} '\' users{indUser} ] ; %change path if needed
    EEG=EEG_cnt2mat( [Directory '1ff01'])
    %fid = fopen( [Directory '1ff02.dat']);
    
    % EEG is a structure with
    %              Fs: scalar (sample rate in Hz)
    %         Trigger: [nb samples x1 ] Trigger channel of '0' with '1' at the start
    %                   of each sweep. There are [nb epochs] '1'.
    %      EpochClass: [nb epochs x1] class of the sweeps (0 for Non-TARGET, 1
    %                   for TARGET).
    %        Channels: [nb samples x nb channels] preprocessed EEG recordings
    %   NoiseTrigger*: the equivalent of Trigger but for the sweep of the noise.
    %                   By default, it takes the same.
    Artefacts=[];
    
    %%%%%%%%%%%%%%%%%% PREPROCESSING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decimationfactor=4; %(5) put 1 to do nothing
    NOTCH=1; %(6) put 1 to remove 50Hz
    BANDPASS=1; %(7) put 1 to BANDPASS (see the filter below for the cutoff freq)
    f1=0.75; %(8) low cutoff freq  (bandpass)
    f2=30; %(9) high cutoff freq  (bandpass)
    N=4; %(10) filter order (bandpass)
    
    %%%%%%%%%%%%%%%%%%%%%% PREPROCESSING (no user input needed)%%%%%%%%%%%%%%%%
    [EEG.Channels EEG.Fs EEG.Trigger]=preprocessingEEG(double(EEG.Channels),EEG.Fs,[f1 f2 N decimationfactor],EEG.Trigger);
    
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
    winElec=findElectrodes(EEG.ElectrodesName,{'P3','P4','O1','O2','T3','T4','T5','T6','T5''','T6''','O1''','O2''','P3"','P4"','PZ"','OZ','I','CB1"','CB2"','CB1','CB2'})
    winTime=[floor((0.05)*EEG.Fs):ceil((0.550)*EEG.Fs)]-offset; %(14)* time window (in sample) used for latency calculation and Pz selection
    % exemple: 50ms to 550ms
    %   *optional
    
    % WARNING: FOLLOWING LINE ONLY TO COPY THE .MAT FILES TO .TXT FILES
    %FULLDIR=['D:\Mes Documents GIPSA\MATLAB\LK_TOOLBOX\data\ARNAUD\txt\'  num2str(indUser) '\' ]
    %EEG_mat2txt(FULLDIR,EEG)
    
    nameSave=[nameSave '_' num2str(EEG.Fs) 'w' num2str(f1) '-' num2str(f2) '_Lat' num2str(Delays)]
    EEG
    ACSTPoptions.SubspaceDim=size(EEG.Channels,2):-1:size(EEG.Channels,2)-6;
    ACSTPoptions.Epoch_size=Window;
    ACSTPoptions.LatencyCorr_max=Delays;
    ACSTPoptions.Mask_Electrodes=winElec;
    ACSTPoptions.Mask_Time=winTime;
    ACSTPoptions.SubspaceDim=[30];
    ACSTPoptions.computeClassLat=[0 1] %Compute only for both Class
    ACSTPoptions.Weights=1
    %% ACSTP loop Algorithm
    [Epochs{indUser} FILTER(indUser)]=ACSTP(EEG,ACSTPoptions)
    
    save(['ARNAUD_' nameSave '.mat'],'Epochs','FILTER')
    
    
    %% visualization of AEA vs CSTP
    %close all
    indSubject=indUser%1:length(FILTER)
    Scale=10 %init scale
    for classIND=1:2
        
        for IND=1:2
            nbLINE=9
            nbCOLUMN=5
            FontSize=24;
            
            Mapping=(EEG.ElectrodesName');
            
            t = (0:(ACSTPoptions.Epoch_size-1))./EEG.Fs;
            PZs=[size(EEG.Channels,2):-1:1];
            winElec2=winElec;
            winTime2=winTime;
            [Pzind]=find(FILTER(indSubject).BestPz(classIND)==PZs)
            
            switch IND
                case 1
                    %PLOTX=mean(X(:,:,EpochClass==1),3);
                    PLOTX=FILTER(indSubject).EA(:,:,classIND);
                    info=['AEA'];
                    if FILTER(indSubject).Class(classIND)==0
                        info=[info ' NT'];
                    else
                        info=[info ' TA']
                    end
                    %Scale=FindScaling(PLOTX,winElec,winTime);
                case 2
                    %PLOTX=Xcstp(:,:,EpochClass==1);
                    PLOTX=FILTER(indSubject).EAcstp(:,:,classIND);
                    info=['ACSTP(' num2str(PZs(Pzind)) ')'];
                    
            end
            xticklabels = (offset/EEG.Fs):.5:(ACSTPoptions.Epoch_size+offset)/EEG.Fs;
            
            clims=[0 7]
            ClassTag=1
            hfig=figure(classIND)
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
        
        if 1%indSubject<length(FILTER)
            saveas(hfig,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\GONOGO\EAvsCSTP\' nameSave '_Subject_' num2str(indSubject) 'class' num2str(classIND) '.tiff'])
            %saveas(hfig,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\GONOGO\EAvsCSTP\' nameSave '_Subject_' num2str(indSubject) 'class' num2str(classIND) '.fig'])
        else
            saveas(hfig,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\GONOGO\EAvsCSTP\' nameSave '_GA_class' num2str(classIND) '.tiff'])
            %saveas(hfig,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\GONOGO\EAvsCSTP\' nameSave '_Subject_GA_class' num2str(classIND) '.fig'])
        end
        
    end
    %%
    
    
    %% TOPOPLOT FOR THE COMPONENTS
    CompLABELS=Generate_Components_Numbers(1:size(EEG.Channels,2))'
    for classIND=1:length(FILTER(indSubject).Class)
        htopo=figure
        
        %[Bs]=FILTER(indSubject).fullBs{FILTER(2).BestPz(classIND)}{classIND};
        [Bs]=FILTER(indSubject).fullAs{classIND};
        Bs=normPOW(Bs');
        %[Bt]=FILTER(indSubject).fullBt{FILTER(2).BestPz(classIND)}{classIND};
        [Bt]=FILTER(indSubject).fullAt{classIND};
        Bt=(normEEG(Bt(:,1:size(Bs,1))'));
        subplot(4,5,[1 6 11 16]);plotEEG(Bt,4,EEG.Fs,CompLABELS);xlabel('Time (s)')
        text(0.35,1.05,['ss ' num2str(indSubject)], 'FontSize',FontSize*1.5,'fontname','times new roman','FontAngle','italic','units','normalized')
        for i=1:size(Bs,1)
            Count=[2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20 22 23 24 25 27 28 29 30 32 33 34 35 37 38 39 40];
            %Count=[2 4 6 8 11 13 15 17 20 22 24 26 29 31 33 35];
            subplot(9,5,[Count(i)])
            topoplot(abs(Bs(i,:)),'GONOGO_locfile.loc','plotrad',0.55);title(CompLABELS{i},'FontSize',FontSize)
            set(htopo, 'PaperPosition', [0 0 40 60])
            caxis([0 1])
        end
        hb=colorbar
        set(hb,'Units','normalized', 'position', [0.92 0.3 0.02 0.6],'FontSize',FontSize);
        colormap(gray)
        colormap(flipud(colormap))
        
        %spaceplots
        saveas(htopo,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\GONOGO\CSTP_components\' nameSave '_Subject_' num2str(indSubject) '_components_Pz' num2str((FILTER(indSubject).BestPz(classIND))) 'class' num2str(classIND) '.tiff'])
        %saveas(htopo,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\GONOGO\CSTP_components\' nameSave '_Subject_' num2str(indSubject) '_components_Pz' num2str((FILTER(indSubject).BestPz(classIND))) 'class' num2str(classIND) '.fig'])
        
        
        %close all
    end
    %% HISTOGRAMS FOR THE LATENCIES AND WEIGHTS
    
    hfighisto=figure
    
    set(hfighisto, 'PaperPosition', [0 0 25 40],'units','normalized','outerposition',[0.6 0.1 0.375 .9])
    subplot(221)
    
    hist(FILTER(indSubject).Latency(EEG.EpochClass==0),(-ACSTPoptions.LatencyCorr_max:+ACSTPoptions.LatencyCorr_max))
    text(1,1.1,['ss' num2str(indSubject)],'FontSize',FontSize,'units','normalized')
    
    xlabel('Latencies NOGO','fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
    FILTER(indSubject).Class(classIND)
    subplot(222)
    hist(FILTER(indSubject).Latency(EEG.EpochClass==1),(-ACSTPoptions.LatencyCorr_max:+ACSTPoptions.LatencyCorr_max))
    xlabel('Latencies GO','fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
    subplot(2,2,[3 4])
    hist(FILTER(indSubject).Weights)
    xlabel('Weights','fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
    saveas(hfighisto,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\GONOGO\EAvsCSTP\' nameSave '_Subject_' num2str(indSubject) '_Latencies_Weights.tiff'])
    %saveas(hfighisto,['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\GONOGO\EAvsCSTP\' nameSave '_Subject_' num2str(indSubject) '_Latencies_Weights.fig'])
end


%% compute grand average for the results
for classIND=1:2
    for indSubject=1:length(users)
        
        EA{classIND}(:,:,indSubject)=FILTER(indSubject).EA(:,:,classIND);
        EAcstp{classIND}(:,:,indSubject)=FILTER(indSubject).EAcstp(:,:,classIND);
        
    end
    EA_GA(:,:,classIND)=mean(EA{classIND},3);
    EAcstp_GA(:,:,classIND)=mean(EAcstp{classIND},3);
end
FILTER(length(users)+1)=FILTER(length(users))
FILTER(length(users)+1).EA=EA_GA
FILTER(length(users)+1).EAcstp=EAcstp_GA


