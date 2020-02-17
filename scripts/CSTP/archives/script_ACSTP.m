% details: use CSTPinfo
%
% *** History: 2015-03-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
%   related functions :
%   LatencyEstimation, weightsEstimation, ACSTP
%
%
% see also : ACSTP, applyCSTP, EnsembleAverage, meanOverlap, epoch_p300, WeightsEstimation,
% CorrectLatency, ConvergenceLatencies, best_Pz, CSTPinfo
clc
close all
clear all


%%%%%%%%%%%%%%%%%% INPUT FILE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% there are data from 4 subjects available in this exemple
Directory= '.\LK_TOOLBOX\data\ERWAN\'; %change path if needed

%subjects available in this exemple : 5, 12, 16, 24
Subjects=21 %select the subject

% Sessions available : 1 to 8 for subject 5, 1 for 12 16 24
load( [Directory 'ss' num2str(Subjects) '_ERWAN_RAW.mat'])
Session=1; %select which session for the given subject


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
 %indice of the subject in 'data' (if several)
%%prepare data

%for erwan data
Phase=find(cellfun(@(x) strcmp(x,'non-adaptive'), data.session{Session}.phasetype));
EEG=data.session{Session}.phase{Phase}.training;


% EEG is a structure with
%             Fs: scalar (sample rate in Hz)
%       StimCode*: [nb epochs x1 double] (useless for CSTP)
%     StimCodePos*: [nb epochs x1 double] (useless for CSTP) 
%           Flash: [nb samples x1 double] position of the sweeps
%          Target*: [nb samples x1 double] position of the sweeps of class
%          TARGET
%               Y: [nb epochs x1 double] class of the sweeps (0 for
%               Non-TARGET, 1 for TARGET)
%               s: [nb samples x nb channels double] raw EEG recording
%               h*: [1x1 struct] software and hardware information (useless for CSTP)
%   *optionnal
Artefacts=[];

%%%%%%%%%%%%%%%%%% INPUT EXTRACTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E=EEG.s'; %(1) full EEG data
    Fs=EEG.Fs; %(2) sample frequency
    Flash=EEG.Flash;  %(3) Sweeps' indices
    Y=EEG.Y; %(4) Class of each sweep

%%%%%%%%%%%%%%%%%% PREPROCESSING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decimationfactor=4; %(5) put 1 to do nothing
    NOTCH=1; %(6) put 1 to remove 50Hz
    BANDPASS=1; %(7) put 1 to BANDPASS (see the filter below for the cutoff freq)
    f1=1; %(8) low cutoff freq  (bandpass)
    f2=20; %(9) high cutoff freq  (bandpass)
    N=4; %(10) filter order (bandpass)

%%%%%%%%%%%%%%%%%% CSTP PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Window=1*Fs; %(11) the sweeps window will be 1s
    Delays=4; %(12) +/- nb shifted samples allowed for the jitter correction

%user parameters for the ACSTP to improve converge :
    winElec=[7,9,10,11,12,13,14,15,16]; %(13)* electrodes used for latency calculation and Pz selection
    % exemple: Cz(7),P7(9),P3(10),Pz(11),P4(12),P8(13),O1(14),Oz(15),O2(16)
    winTime=[floor((0.05)*128):ceil((0.550)*128)]; %(14)* time window (in sample) used for latency calculation and Pz selection
    % exemple: 50ms to 550ms

%   *optional

%%%%%%%%%%%%%%%%%%%%%% PREPROCESSING (no user input needed)%%%%%%%%%%%%%%%%
[Edec Fsd Flashd]=preprocessingEEG(E',Fs,[f1 f2 N decimationfactor],Flash);
Fs=Fsd;
Window=Window/decimationfactor;
clear E;
E=Edec';
Flash=Flashd;

%%%%%%% ADDITIONNALS STRUCTURES FOR ACSTP (no user input needed)%%%%%%%%%%%
WindowB=Window+2*Delays;
offset=-Delays;
[Xall]=epoch_p300(E,Flash,WindowB,offset); %prepare epoch for full window+latency
%Xall is used during jitter correction
%[Flash Y]= removeEpochArtifact(Xall,Flash,Y); %remove artifacts by simple thresholding (here 250 µV)


%% ACSTP loop Algorithm
clear Zbarz Conv FlashCorrected P LatencyCorrection Weights Xhatk As At Bs Bt iterLat cond Bs As Bt At X Save
% Save (cell array with the following items) :
% (1) AEA
% (2) AEAcstp{iter} after weights, latency corrections and CSTP
% (3) CSTP{iter} Bs Bt As At : the temporal and spatials filters
% (4) Pzind : optimal indice of iter such as Pz is the optimal subspace
%       with respect of the mask
%
% Thus Save{2}{Save{4}) gives P300 with optimal estimation.
disp(sprintf('COMPUTE ACSTP... STATE11CHAR'))
printCurrentState(sprintf('START     '))

CriteriaConvLatencies=1/length(find(Y))*Delays; %nb latencies correction allowed for convergence
maxIter_Latcor=10;
Pzfs=[size(Xall,1):-1:10];%choose the screening of the Pz
%Y=ones(size(Y));
[Xbarz Class]=meanOverlap(E,Flash,Window,Y);%estimate the arithmetic mean (0.10)
Save{1}=Xbarz; %save the AEA for future visualization

iter=1;%for each iteration, a Pz
for Pzf=Pzfs;
    %%
    
    printCurrentState(['Pzf ' num2str(Pzf) '      '])
    
    %%(A.1) #################### DATA PREPARATION ####################
    X=epoch_p300(E,Flash,Window); %get epochs
    
    %if needed, we can estimate the noise covariance from differents epochs
    %(e.g. resting state or randomly selected epochs)
    Xnoise=X;
    Ynoise=Y;
    %%(A.2) #################### LATENCY INITIALIZATION ####################
        LatencyCorrection(:,iter)=zeros(length(Y),1); %initialization latency at zero

    %%(A.3) #################### WEIGHTS INITIALIZATION ####################
    Weights=WeightsEstimation(X,Y);
    
    %%(A.4) #################### Xbarz INITIALIZATION ####################
    [Xbarz Class]=meanOverlap(E,Flash,Window,Y,Weights);
    
    
    %%(B) #################### ACSTP Initialization ####################
    printCurrentState(['Pzf ' num2str(Pzf) 'ACSTP '])
    [Bs Bt As At Class eigV]=CSTP(X,Y,Xbarz,Weights,winElec,winTime,Pzf,Xnoise,Ynoise);
    %%(3) #################### Zbarz Initialization ####################
    Zbarz{iter}=applyCSTP(Xbarz,Bs,Bt,As,At,Class);
    
    Save{2}{iter}=Zbarz{iter}; %save after CSTP (initialized)
    Save{3}{iter}={Bs Bt As At};%save CSTP (initialized)
    if iter>1 %Compute everything just for Pz<max(Pz)
        
        printCurrentState(['Pzf ' num2str(Pzf) 'WEIGHT'])
        
        %%(3) #################### WEIGHTS ESTIMATION ####################
        [Weights Xhatk]=WeightsEstimation(X,Y,Bs,Bt,As,At);
        
        %%(4) #################### LATENCY CORRECTION LOOP ####################
        classLat=Y==1;%just compute latencies for class target
        tempoLatency=[];
        STOPlat=1;
        iterLat=1;
        printCurrentState(['Pzf ' num2str(Pzf) 'LATENC'])
        
        LatencyCorrection(classLat,iter)=zeros(size(LatencyCorrection(classLat,iter-1))); %initialize latencies
        while STOPlat %converge criteria
            
            [tempoLatency Crit_Trace]=LatencyEstimation(Xall(:,:,classLat),X(:,:,classLat),Y(classLat),Weights(classLat),Bs(end),Bt(end),As(end),Bt(end),winElec,winTime);
            
            Conv(iter,iterLat)=ConvergenceLatencies(tempoLatency,LatencyCorrection(classLat,iter));
            if (iterLat>1 && Conv(iter,iterLat)<CriteriaConvLatencies) || iterLat>=maxIter_Latcor
                STOPlat=0;
            end
            LatencyCorrection(classLat,iter)=tempoLatency;
            FlashCorrected=CorrectLatency(Flash,LatencyCorrection(:,iter));
            X=epoch_p300(E,FlashCorrected,Window);
            iterLat=iterLat+1;
            
        end
        
        
        
        %%(4) #################### CORRECT LATENCIES ####################
        FlashCorrected=CorrectLatency(Flash,LatencyCorrection(:,iter));
        X=epoch_p300(E,FlashCorrected,Window);
        
        [Xbarz Class]=meanOverlap(E,FlashCorrected,Window,Y,Weights);
        % Xbarz is then the EAE averaged with weighted epochs and with
        % corrected lattencies
        
        printCurrentState(['Pzf ' num2str(Pzf) 'FINALI'])
        
        %%(3) #################### ACSTP FINAL ####################
        [Bs Bt As At Class eigV]=CSTP(X,Y,Xbarz,Weights,winElec,winTime,Pzf,Xnoise,Ynoise);
        
        %%(3) #################### FINALIZATION ####################
        % for visualization after weights, latency correction and CSTP
        Zbarz{iter}=applyCSTP(Xbarz,Bs,Bt,As,At,Class);      
        Save{2}{iter}=Zbarz{iter};
        Save{3}{iter}={Bs Bt As At};
    end
    iter=iter+1;
    
    
end
%find best Pz
    [Pzind]=best_Pz(Save{2},winElec,winTime)
    Save{4}=Pzind;
%% visualization of AEA vs CSTP

close all
Scale=1 %init scale
classIND=1
for IND=1:2
    nbLINE=9
    nbCOLUMN=5
    FontSize=24;
    Mapping={'Fp1';%1
        'Fp2';%2
        'F5';%3
        'AFz';%4
        'F6';%5
        'T7';%6
        'Cz';%7
        'T8';%8
        'P7';%9
        'P3';%10
        'Pz';%11
        'P4';%12
        'P8';%13
        'O1';%14
        'Oz';%15
        'O2'}%16
    t = (0:(128-1))./128;
    PZs=[16:-1:1];
    winElec2=[7 11 14 15 16];
    winTime2=5:60;
    [Pzind SNR]=best_Pz(Save{2},winElec,winTime,classIND);
    [Bs]=Save{3}{Pzind}{1};
    [Bt]=Save{3}{Pzind}{2};
    [As]=Save{3}{Pzind}{3};
    [At]=Save{3}{Pzind}{4};

    Yw=Y;
    
    Xcstp=applyCSTP(X,Bs,Bt,As,At,Yw);
    switch IND
        case 1
            %PLOTX=mean(X(:,:,Y==1),3);
            PLOTX=Save{1}(:,:,classIND);
            info=['ss' num2str(Subjects) ' AEA'];
            %Scale=FindScaling(PLOTX,winElec,winTime);
        case 2
            %PLOTX=Xcstp(:,:,Y==1);
            PLOTX=mean(Save{2}{Pzind}(:,:,classIND),3);
            info=['ss' num2str(Subjects) ' ACSTP(' num2str(PZs(Pzind)) ')'];
            
    end
    xticklabels = 0:.5:1;
    
    clims=[0 7]
    ClassTag=1
    hfig=figure(2)
    set(hfig, 'PaperPosition', [0 0 12.5 30],'units','normalized','outerposition',[0.6 0.1 0.375 .9])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT EEG %%%%%%%%%%%%%%%
    subplot(nbLINE,nbCOLUMN,[1 2 6 7 11 12 16 17 21 22 26 27 31 32 36 37]+2*(IND-1));plotEEG(PLOTX,Scale,128,Mapping) % Xbar TA
    %set(gca,'yticklabel',[])
    if IND>1
        set(gca, 'color', [0.95 .95 .95])
    end
    set(gcf, 'color', [1 1 1])
    
    title(info,'FontSize',FontSize)
    %xlabel('Time (s)')
    xticks = linspace(0, 1.01, numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
    set(gcf, 'InvertHardCopy', 'off');
    % global field power on average
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT GFP mean %%%%%%%%%%%%%%%
    subplot(nbLINE,nbCOLUMN,[21 22]+20+2*(IND-1))
    area(t,global_field_power(PLOTX'),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
    ylim([0 6]);grid on;
    if IND==1
        %ylabel('GFP','FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
    end
    set(gca,'YtickLabel',[])
    xticks = linspace(0, 1, numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
    set(gca, 'color', [0.95 .95 .95])
    
    
    spaceplots
    
end
