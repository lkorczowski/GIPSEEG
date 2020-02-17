function [Xhat ACSTPstruct]=ACSTP(EEG,ACSTPoptions)
% Compute the Adaptive Common Spatio-Temporal Pattern Filter (ACSTP) from the
% continuous recording EEG.
%
% The estimated ERP epochs Xhat are returned such as :
% for X the EEG epochs with corrected latencies and with weights W, the ACSTP filter, does:
% Xhat(:,:,k)=As*Bs'*W(k)*X(:,:,k)*Bt*At';
%
% With the output ACSTPstruct, you can manually apply the filter on the EEG
% such as :
% TriggerCorrected=CorrectLatency(EEG.Trigger,ACSTPstruct.Latency); % apply latency correction
% X=epoch_p300(EEG.Channels,TriggerCorrected,ACSTPstruct.Epoch_size); % extract the epochs
% Xw=applyWeights(X,ACSTPstruct.Weights); %apply the weights
% Xhat=applyCSTP(Xw,ACSTPstruct.Bs,ACSTPstruct.Bt,ACSTPstruct.As,ACSTPstruct.At,EEG.EpochClass); %apply the ACSTP
%
%
% INPUTS :
% ------
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
%
% ACSTPoptions is a structure with
%      Epoch_size: scalar, the length of the epoch window (in samples)
% LatencyCorr_max: scalar, the maximum of samples for the latency
%                   correction. Set 0 to disable the Latency correction.
% Mask_Electrodes: vector of the selectionned electrodes. Usefull for the
%                  automatic subspace selection (BestPz) and latency
%                  correction (Latency).
%       Mask_Time: vector of the selectionned sample. Usefull for the
%                  automatic subspace selection (BestPz) and latency
%                  correction (Latency).
% MaxIterLatency*: scalar, the maximum iteration allowed to compute the
%                   latency correction. Default: 10.
%    SubspaceDim*: vector, containing all the subspace dimension (<nb electrodes)
%                  to test in descent order.
%                   By default, it is equal to (nb_channels:-1:(nb_channels/2))
%computeClassLat*: vector, containing all the class tag in which you want
%                   to compute the latency correction. By default, it
%                   computes it for all classes but it could be long (for
%                   instance you can skip the non-target).
%        Weights*: Default: true.
%                  option1(given) [nb epochs x1] vector, containing the weights for each
%                   epoch if it is computed from an external function.
%                  option2 (true/false) boolean. If true (default) the ACSTP compute the
%                   weight for each epoch. If false, the weights are
%                   desactivated (i.e. set to 1).
%        DISPLAY*: Boolean, should the comparative result between the arithmetic ensemble
%                   average and the ACSTP should be display at the end. Default: true.
%
%   *optional
%
% OUTPUT :
% ------
% Xhat with corrected latencies, weighted and filtered
% ACSTPstruct is a structure with
%              EA: the ensemble average before ACSTP
%          EAcstp: the ensemble average corrected with latencies, weighted
%                   and filtered + with the effect of overlapping
%                   correction
%     As Bs Bt At: such as Xhat(:,:,k)=As{indClass}*Bs{indClass}'*W(k)*X(:,:,k)*Bt{indClass}*At{indClass}'
%                   Each filter is a cell containing a matrix filter for
%                   each class
%           Class: The tag and the order in which the filter are sorted
%          BestPz: Orders of the best subspace for the ACSTP for the given
%                  Class
%         Weights: [nb epochs x1] the weights of each epoch
%         Latency: [nb epochs x1] the corrected offset of each epoch
%      Epoch_size: scalar, the length of the epoch window (in samples)
%
%
% *** History: 16-Apr-2015
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
% see also : CSTP, applyCSTP, EnsembleAverage, meanOverlap, epoch_p300, WeightsEstimation,
% CorrectLatency, ConvergenceLatencies, best_Pz, CSTPinfo
Xhat=[];
ACSTPstruct=struct;
if ~isfield(ACSTPoptions,'DISPLAY')
    DISPLAY=1; %put '1' to plot the results at the end of the computation or '0' to not.
else
    DISPLAY=ACSTPoptions.DISPLAY;
end

DISPLAYtime=0;%set on for debug (optimization)
if( DISPLAYtime) tic; end
%%%%%%%%%%%%%%%%%% INPUT EXTRACTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEGsignals=EEG.Channels'; %(1) full EEG data
Fs=EEG.Fs; %(2) sample frequency
Trigger=EEG.Trigger;  %(3) Sweeps' indices
EpochClass=EEG.EpochClass; %(4) Class of each sweep

%Check if the EEG INPUTS are correct
if size(EEGsignals,2)<size(EEGsignals,1)
    warning('ACSTP WARNING1: EEG INPUT ERROR: the temporal dimension of EEG.Channels is smaller than the spatial dimension.')
end

if ~(size(EEGsignals,2)==length(Trigger))
    warning('ACSTP WARNING2: EEG INPUT ERROR. Size of EEG.Trigger')
    if length(Trigger)==length(EpochClass)
        %we assume that Trigger gives the position of the EpochClass instead of
        %being a trigger channel at the sample rate Fs
        tmp=zeros(size(EEGsignals,2),1);
        tmp(Trigger)=1;
        Trigger=tmp;
        clear tmp
        warning('EEG.Trigger has been assumed to be the EEG.EpochClass locations instead of a trigger channel at the sample rate Fs')
    else
        error('ACSTP ERROR3: EEG INPUT ERROR.')
        
    end
end


%%%%%%%%%%%%%%%%%% CSTP PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<2
    Window=1*Fs; %(11) the sweeps window will be 1s
    Delays=4; %(12) +/- nb shifted samples allowed for the jitter correction
    winElec=[7,9,10,11,12,13,14,15,16]; %(13)* electrodes used for latency calculation and Pz selection
    % exemple: Cz(7),P7(9),P3(10),Pz(11),P4(12),P8(13),O1(14),Oz(15),O2(16)
    winTime=[floor((0.05)*128):ceil((0.550)*128)]; %(14)* time window (in sample) used for latency calculation and Pz selection
    % exemple: 50ms to 550ms
else
    Window=ACSTPoptions.Epoch_size;
    Delays=ACSTPoptions.LatencyCorr_max;
    if ~isfield(ACSTPoptions,'Mask_Electrodes')
        winElec=1:size(EEGsignals,1);%put all electrodes
    elseif isempty(ACSTPoptions.Mask_Electrodes)
        winElec=1:size(EEGsignals,1);%put all electrodes
    else
        winElec=ACSTPoptions.Mask_Electrodes;
    end
    
    if ~isfield(ACSTPoptions,'Mask_Time')
        winElec=1:Window;%put all samples
    elseif isempty(ACSTPoptions.Mask_Time)
        winElec=1:Window;%put all samples
        
    else
        winTime=ACSTPoptions.Mask_Time;
    end
end
%Check if the ACSTP INPUTS are correct
if find(Trigger,inf,'last')+Window>size(EEGsignals,2)
    error('ACSTP ERROR4: ACSTP INPUT ERROR. The size of the epochs choosen in ACSTPoptions.Epoch_size is too big for the last event in EEG.Trigger.');
    
end
if (find(Trigger,inf,'last')+Window+Delays>size(EEGsignals,2))
    error('ACSTP ERROR5: ACSTP INPUT ERROR. The delays choosen in ACSTPoptions.LatencyCorr_max is too big for the last event in EEG.Trigger.');
    
end
if (find(Trigger,1)-Delays<1)
    error('ACSTP ERROR6: ACSTP INPUT ERROR. The delays choosen in ACSTPoptions.LatencyCorr_max is too big for the first event in EEG.Trigger.');
    
end
if (max(winElec)>size(EEGsignals,1))
    error('ACSTP ERROR7: ACSTP INPUT ERROR. ACSTPoptions.Mask_Electrodes incorrect');
    
end
if (max(winTime)>Window)
    error('ACSTP ERROR8: ACSTP INPUT ERROR. ACSTPoptions.Mask_Time incorrect');
    
end
%%%%%%% ADDITIONNALS STRUCTURES FOR ACSTP (no user input needed)%%%%%%%%%%%
WindowB=Window+2*Delays;
offset=-Delays;
[Xall]=epoch_p300(EEGsignals,Trigger,WindowB,offset); %prepare epoch for full window+latency
if ~isfield(ACSTPoptions,'MaxIterLatency')
    maxIter_Latcor=10;
else
    maxIter_Latcor=ACSTPoptions.MaxIterLatency
end
if ~isfield(ACSTPoptions,'SubspaceDim')
    Pzfs=[size(Xall,1):-1:size(Xall,1)/2];%choose the screening of the Pz
else
    if max(ACSTPoptions.SubspaceDim)>size(EEGsignals,1)
        error('ACSTP ERROR9: ACSTP INPUT ERROR. The subspace dimension(s) choosen in ACSTPoptions.SubspaceDim is/are greater than the number of electrodes.');
        return
    elseif min(ACSTPoptions.SubspaceDim)<1
        error('ACSTP ERROR9: ACSTP INPUT ERROR. The subspace dimension(s) choosen in ACSTPoptions.SubspaceDim is/are smaller than 1.');
        return
    else
        Pzfs=ACSTPoptions.SubspaceDim;
    end
end

%{
%the constant offset has been desactivated. If you'll like to add an
offset, please consider modifying Trigger as well as Epoch_size and Mask_Time
accordingly in ACSTPoptions.

if ~isfield(ACSTPoptions,'EpochsOffset')
    EpochsOffset=0;%offset has been set to zero
else
    EpochsOffset=ACSTPoptions.EpochsOffset;
    disp(['ACSTP INFO: The epochs will be with an offset of: ' num2str(EpochsOffset) ' samples.']);
    if (find(Trigger,1)-Delays+EpochsOffset<1)
        disp('ACSTP ERROR9: ACSTP INPUT ERROR. ACSTPoptions.EpochsOffset is too big for the first event in EEG.Trigger.');
    end
end
%}

if ~isfield(ACSTPoptions,'computeClassLat')
    ClassLat=unique(EpochClass);%all the classes
    disp(['ACSTP INFO: Latencies will be computed for class(es): ' num2str(ClassLat) '.']);
    disp('ACSTP INFO: Please consider only classes that need Latency correction to speed up the procedure.');
    disp( 'ACSTP INFO: See setting ACSTPoptions.computeClassLat');
    
else
    ClassLat=ACSTPoptions.computeClassLat;
    disp(['ACSTP INFO: Latencies will be computed only for class(es): ' num2str(ClassLat)]);
    
end

% CHECK the Weights' estimation (true/false/given).
if ~isfield(ACSTPoptions,'Weights') % true (default)
    ComputeWeights=1; %ACSTP will estimate the Weights
    disp(['ACSTP INFO: Weights will be computed: ' num2str(ComputeWeights) ' (1=true/0=false)']);
else
    if(length(ACSTPoptions.Weights)==length(EpochClass)) %given
        ComputeWeights=0; %ACSTP will NOT estimate the Weights
        Weights=ACSTPoptions.Weights;
        warning(['ACSTP WARNING10: ACSTP INPUT WARNING. ACSTPoptions.Weights has been given and won''t be computed.'...
            ' Please check carefully your weights to avoid a scaling issue in output of the CSTP.']);
    elseif(length(ACSTPoptions.Weights)==1) %set true/false
        ComputeWeights=ACSTPoptions.Weights;%ACSTP does(true)/doesn't(false) estimate the Weights
        disp(['ACSTP INFO: Weights will be computed: ' num2str(ComputeWeights) ' (1=true/0=false)']);
    else %wrong size
        error('ACSTP ERROR11: ACSTP INPUT ERROR. ACSTPoptions.Weights size incorrect');
    end
    
end
if (~exist('Weights','var') && ~ComputeWeights) %if Weights estimation false
    Weights=ones(length(EpochClass),1);
end
%% ACSTP loop Algorithm

disp(['ACSTP INFO: Latency estimation: ' num2str(~(Delays==0)) ' (1=activated/0=disabled)']);

% Heuristic criterion for the Latency convergence:
CriteriaConvLatencies=1/length(find(EpochClass))*Delays; %nb latencies correction allowed for convergence.


[ACSTPstruct.EA]=EnsembleAverage(EEGsignals,EpochClass,Trigger,Window);

%EpochClass=ones(size(EpochClass));
[Xbarz Class]=meanOverlap(EEGsignals,Trigger,Window,EpochClass);%estimate the arithmetic mean (0.10)
Output{1}=Xbarz; %save the AEA for future visualization
if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
disp((['COMPUTE ACSTP... STATE20CHAR         ' dtime]))
printCurrentState(['START               ' dtime])
iter=1;%for each iteration, a Pz
for Pzf=Pzfs;
    %%
    
    
    %%(A.1) #################### DATA PREPARATION ####################
    X=epoch_p300(EEGsignals,Trigger,Window); %get epochs
    
    %if needed, we can estimate the noise covariance from differents epochs
    %(e.g. resting state or randomly selected epochs)
    Xnoise=X;
    Ynoise=EpochClass;
    %%(A.2) #################### LATENCY INITIALIZATION ####################
    
    LatencyCorrection(:,iter)=zeros(length(EpochClass),1); %initialization latency at zero
    
    %%(A.3) #################### WEIGHTS INITIALIZATION ####################
    if (ComputeWeights)
        Weights=WeightsEstimation(X,EpochClass);
    end
    
    %%(A.4) #################### Xbarz INITIALIZATION ####################
    
    [Xbarz Class]=meanOverlap(EEGsignals,Trigger,Window,EpochClass,Weights);
    
    if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
    %%(B) #################### ACSTP Initialization ####################
    printCurrentState(['Pzf ' num2str(Pzf) 'ACSTP          ' dtime])
    [Bs Bt As At Class eigV]=CSTP(X,EpochClass,Xbarz,Weights,winElec,winTime,Pzf,Xnoise,Ynoise);
    %%(3) #################### Zbarz Initialization ####################
    Zbarz{iter}=applyCSTP(Xbarz,Bs,Bt,As,At,Class);
    
    Output{2}{iter}=Zbarz{iter}; %save after CSTP (initialized)
    Output{3}{iter}={Bs Bt As At};%save CSTP (initialized)
    
    if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
    printCurrentState(['Pzf ' num2str(Pzf) 'WEIGHT         ' dtime])
    
    %%(3) #################### WEIGHTS ESTIMATION ####################
    if (ComputeWeights)
        Weights=WeightsEstimation(X,EpochClass,Bs,Bt,As,At);
    end
    
    %%(4) #################### LATENCY CORRECTION LOOP ####################
    for indCl=1:length(ClassLat)
        class_for_lat=ClassLat(indCl);
        classLat_indices=(EpochClass==class_for_lat);%just compute latencies for class target
        tempoLatency=[];
        STOPlat=1;
        iterLat=1;
        if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
        printCurrentState(['Pzf ' num2str(Pzf) 'LATEN' num2str(class_for_lat) '         ' dtime])
        
        LatencyCorrection(classLat_indices,iter)=zeros(size(LatencyCorrection(classLat_indices,iter))); %initialize latencies
        if Delays>0
            while STOPlat %converge criteria
                
                [tempoLatency Crit_Trace]=LatencyEstimation(Xall(:,:,classLat_indices),X(:,:,classLat_indices),EpochClass(classLat_indices),Weights(classLat_indices),Bs(end),Bt(end),As(end),Bt(end),winElec,winTime);
                
                Conv(iter,iterLat)=ConvergenceLatencies(tempoLatency,LatencyCorrection(classLat_indices,iter));
                if (iterLat>1 && Conv(iter,iterLat)<CriteriaConvLatencies) || iterLat>=maxIter_Latcor
                    STOPlat=0;
                end
                LatencyCorrection(classLat_indices,iter)=tempoLatency;
                TriggerCorrected=CorrectLatency(Trigger,LatencyCorrection(:,iter));
                X=epoch_p300(EEGsignals,TriggerCorrected,Window);
                iterLat=iterLat+1;
                
            end
        else
            if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
            
            printCurrentState(['NOLAT_esti         ' dtime])
        end
        
    end
    
    
    %%(4) #################### CORRECT LATENCIES ####################
    TriggerCorrected=CorrectLatency(Trigger,LatencyCorrection(:,iter));
    X=epoch_p300(EEGsignals,TriggerCorrected,Window);
    
    [Xbarz Class]=meanOverlap(EEGsignals,TriggerCorrected,Window,EpochClass,Weights);
    
    % Xbarz is then the EAE averaged with weighted epochs and with
    % corrected lattencies
    
    if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
    printCurrentState(['Pzf ' num2str(Pzf) 'FINALI         ' dtime])
    
    %%(3) #################### ACSTP FINAL ####################
    
    [Bs Bt As At Class eigV fullBs{iter} fullBt{iter} fullAs{iter} fullAt{iter}]=CSTP(X,EpochClass,Xbarz,Weights,winElec,winTime,Pzf,Xnoise,Ynoise);
    
    %%(3) #################### FINALIZATION ####################
    % for visualization after weights, latency correction and CSTP
    Zbarz{iter}=applyCSTP(Xbarz,Bs,Bt,As,At,Class);
    Output{2}{iter}=Zbarz{iter};
    Output{3}{iter}={Bs Bt As At};
    Output{4}{iter}=Weights;
    %end
    iter=iter+1;
    
    
end
% find the best Pz for each class
nbClass=size(Output{2}{1},3) ;
ACSTPstruct.EAoverlap=Output{1};%standard arithmetic ensemble average
for indClass=1:size(Output{2}{1},3)
    
    [Pzind SNR]=best_Pz(Output{2},winElec,winTime,indClass);
    
    ACSTPstruct.Bs{indClass}=Output{3}{Pzind}{1}{indClass};
    ACSTPstruct.Bt{indClass}=Output{3}{Pzind}{2}{indClass};
    ACSTPstruct.As{indClass}=Output{3}{Pzind}{3}{indClass};
    ACSTPstruct.At{indClass}=Output{3}{Pzind}{4}{indClass};
    ACSTPstruct.EAcstp(:,:,indClass)=mean(Output{2}{Pzind}(:,:,indClass),3);
    ACSTPstruct.BestPz(indClass)=Pzfs(Pzind);
end
ACSTPstruct.Latency=LatencyCorrection(:,Pzind);
ACSTPstruct.Weights=Output{4}{Pzind};
ACSTPstruct.Class=Class;
ACSTPstruct.fullBs=fullBs{Pzind};
ACSTPstruct.fullBt=fullBt{Pzind};
ACSTPstruct.fullAs=fullAs{Pzind};
ACSTPstruct.fullAt=fullAt{Pzind};


% do the final filtering from initial data from the optimal
% SubspaceDimension
TriggerCorrected=CorrectLatency(Trigger,ACSTPstruct.Latency);
X=epoch_p300(EEGsignals,TriggerCorrected,Window); %apply latency correction
Xw=applyWeights(X,Output{4}{Pzind}); %apply the weights
Xhat=applyCSTP(Xw,ACSTPstruct.Bs,ACSTPstruct.Bt,ACSTPstruct.As,ACSTPstruct.At,EpochClass); % apply the ACSTP

if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
printCurrentState(['ACSTP has finished  ' dtime])

% ACSTPstruct is a structure with
% EnsembleAverage: the ensemble average corrected with latencies, weighted
%                   and filtered + with the effect of overlapping
%                   correction
%     As Bs Bt At: such as Xhat(:,:,k)=As{indClass}*Bs{indClass}'*W(k)*X(:,:,k)*Bt{indClass}*At{indClass}'
%                   Each filter is a cell containing a matrix filter for
%                   each class
%           Class: The tag and the order in which the filter are sorted
%         Weights: [nb epochs x1] the weights of each epoch
%         Latency: [nb epochs x1] the corrected offset of each epoch



%% find best Pz and PLOT
if DISPLAY
    
    for i=1:size(Output{2}{1},3) %plot each class results
        
        %% visualization of AEA vs CSTP
        %    close all
        Scale=10 %init scale
        for classIND=1:2
            
            for IND=1:2
                nbLINE=9
                nbCOLUMN=5
                FontSize=24;
                if size(EEGsignals,1)==16
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
                elseif size(EEGsignals,1)==32
                    Mapping={'Fp1';%1
                        'Fp2';%2
                        'AFz';%3
                        'F7';%4
                        'F3';%5
                        'F4';%6
                        'F8';%7
                        'FC5';%8
                        'FC1';%9
                        'FC2';%10
                        'FC6';%11
                        'T7';%12
                        'C3';%13
                        'Cz';%14
                        'C4';%15
                        'T8';%16
                        'Cp5';%17
                        'Cp1';%18
                        'Cp2';%19
                        'Cp6';%20
                        'P7';%21
                        'P3';%22
                        'Pz';%23
                        'P4';%24
                        'P8';%25
                        'PO7';%26
                        'O1';%27
                        'Oz';%28
                        'O2';%29
                        'PO8';%30
                        'PO9';%31
                        'PO10' };%32
                    
                elseif isfield(EEG,'ElectrodesName')
                    Mapping=(EEG.ElectrodesName');
                else
                    Mapping=[]
                    
                end
                t = (0:(Window-1))./Window;
                PZs=(size(EEGsignals,1):-1:1);
                [Pzind]=best_Pz(Output{2},winElec,winTime,classIND);
                
                switch IND
                    case 1
                        %PLOTX=mean(X(:,:,EpochClass==1),3);
                        PLOTX=Output{1}(:,:,classIND);
                        info=['AEA'];
                        if Class(classIND)==0
                            info=[info ' NT'];
                        else
                            info=[info ' TA']
                        end
                        %Scale=FindScaling(PLOTX,winElec,winTime);
                    case 2
                        %PLOTX=Xcstp(:,:,EpochClass==1);
                        PLOTX=mean(Output{2}{Pzind}(:,:,classIND),3);
                        info=['ACSTP(' num2str(PZs(Pzind)) ')'];
                        
                end
                xticklabels = 0:.5:1;
                
                clims=[0 7]
                ClassTag=1
                hfig=figure(classIND*1000)
                set(hfig, 'PaperPosition', [0 0 12.5 30],'units','normalized','outerposition',[0.6 0.1 0.375 .9])
                %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT EEG %%%%%%%%%%%%%%%
                subplot(nbLINE,nbCOLUMN,[1 2 6 7 11 12 16 17 21 22 26 27 31 32 36 37]+2*(IND-1));plotEEG(PLOTX,Scale,EEG.Fs,Mapping) % Xbar TA
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
                
                
                %spaceplots
                
            end
        end
    end
end