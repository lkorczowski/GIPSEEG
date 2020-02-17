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
for indUser=1%:length(users) %select the subject
    close all
    nameSave='GONOGO'
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
    Window=round(2*EEG.Fs); %(11) the sweeps window will be 1s
    Delays=10; %(12) +/- nb shifted samples allowed for the jitter correction
    %user parameters for the ACSTP to improve converge :
    %winElec=[7,9,10,11,12,13,14,15,16]; %(13)* electrodes used for latency calculation and Pz selection
    % exemple: Cz(7),P7(9),P3(10),Pz(11),P4(12),P8(13),O1(14),Oz(15),O2(16)
    winElec=findElectrodes(EEG.ElectrodesName,{'P3','P4','O1','O2','T3','T4','T5','T6','T5''','T6''','O1''','O2''','P3"','P4"','PZ"','OZ','I','CB1"','CB2"','CB1','CB2'})
    winTime=[floor((0.05)*EEG.Fs):ceil((0.550)*EEG.Fs)]; %(14)* time window (in sample) used for latency calculation and Pz selection
    % exemple: 50ms to 550ms
       EEG
    ACSTPoptions.SubspaceDim=size(EEG.Channels,2):-1:size(EEG.Channels,2)-12;
    ACSTPoptions.Epoch_size=Window;
    ACSTPoptions.LatencyCorr_max=Delays;
    ACSTPoptions.Mask_Electrodes=winElec;
    ACSTPoptions.Mask_Time=winTime;
    ACSTPoptions.computeClassLat=[0 1] %Compute only for both Class
    ACSTPoptions.Weights=1
    %% ACSTP loop Algorithm
    [Epochs{indUser} FILTER(indUser)]=ACSTP(EEG,ACSTPoptions)
end