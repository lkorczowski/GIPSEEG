% details: use CSTPinfo
%
% *** History: 16-Apr-2015
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
%
% see also : ACSTP, CSTP, applyCSTP, EnsembleAverage, meanOverlap, epoch_p300, WeightsEstimation,
% CorrectLatency, ConvergenceLatencies, best_Pz, CSTPinfo
clc
close all
clear all


%%%%%%%%%%%%%%%%%% INPUT FILE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% there are data from 4 subjects avawilable in this exemple
Directory= '.\SVN\ACSTP\'; %change path if needed

%subjects available in this exemple : 5, 12, 16, 24
Subjects=21 %select the subject

% Sessions available : 1 to 8 for subject 5, 1 for 12 16 24
load( [Directory 'RAW_screening_data_S' num2str(Subjects) '.mat'])
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
%Phase=find(cellfun(@(x) strcmp(x,'non-adaptive'), data.session{Session}.phasetype));
EEG=subject.session(Session);


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
    winElec=[14,21,22,23,24,25,27,28,29]; %(13)* electrodes used for latency calculation and Pz selection
    % exemple: Cz(7),P7(9),P3(10),Pz(11),P4(12),P8(13),O1(14),Oz(15),O2(16)
    winTime=[floor((0.050)*128):ceil((0.550)*128)]; %(14)* time window (in sample) used for latency calculation and Pz selection
    % exemple: 50ms to 550ms

%   *optional


%%%%%%%%%%%%%%%%%%%%%% PREPROCESSING (no user input needed)%%%%%%%%%%%%%%%%
[Edec Fsd Flashd]=preprocessingEEG(E',Fs,[f1 f2 N decimationfactor],Flash);
EEG.Fs=Fsd;
Window=Window/decimationfactor;
clear E;
EEG.s=Edec;
EEG.Flash=Flashd;



    ACSTPoptions.Epoch_size=Window;
    ACSTPoptions.LatencyCorr_max=Delays;
    ACSTPoptions.Mask_Electrodes=winElec;
    ACSTPoptions.Mask_Time=winTime;
    ACSTPoptions.computeClassLat=1
%% ACSTP loop Algorithm
[Epochs FILTER]=ACSTP(EEG,ACSTPoptions);