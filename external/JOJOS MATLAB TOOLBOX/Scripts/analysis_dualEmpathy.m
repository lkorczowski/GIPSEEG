% -------------------------------------------------------------------------
% Analysis of dual-EEG dual-empathy experiment data (Gispa-Lab 2012)
% -------------------------------------------------------------------------
%
% This script performs analysis on data from MULTIPLE couples.
%
% Here we do following analysis procedure:
% 0) Data importation from all subject couples
%
% *** C/ Analysis of behavioral data
% C1) Harmony in couple
% C2) Emotion intensity and perceived coupling
% C3) Empathy questionnaire
% C4) Interactions between items
%
% *** B/ Analysis of physiological signals
% B1) Pulse Rate Variability (PRV)
% B2) Breath-to-Breath Variability (BBV) 
% B3) Electrodermal Activity (EDA) 
%
% *** A/ Group BSS and analysis of component spectra
% A1) Estimate average ICA components for all subjects -> B_global
% A2) Visualize components information and put IC of interest in first positions
% A3) Compute t-test on components spectra estimated in different conditions
% A4) [TODO] Refine BSS: estimate individual ICA components for each subjects,
% using or not using B_global
% A5) [TODO] Visualize components information (spectra, scalp maps) and select
% individual ones that correspond to global ones using time-correlations
% A6) [TODO] Compute statistics on these components? 
% 
% History:
% --- 2013-04-23
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


clear all
close all


addpath(genpath([fileparts(pwd) '\JOJOS MATLAB TOOLBOX'])); % add all subfolders below main functions folder
cfgANA.DataPath = ''; % <-- PUT YOUR PATH HERE (where data is located on your disk)


%% -------------------------- Parametrization
% -- Analysis parameters
cfgANA.couple_set     	= 'all';   % can be 'all', 'best', 'verybest', 'eyesopen', 'eyesclosed', 'worst' or a vector with couples IDs
cfgANA.dataType         = 'EXT_only';   % 'all' or 'EEG_only' or 'EXT_only' : choice of data on which to perform analysis 
cfgANA.procedureA_groupBSS  = 0;        % perform procedure A/ (see above)
cfgANA.procedureB_physio    = 1;    	% perform procedure B/ (see above)
cfgANA.procedureC_behav  	= 0;        % perform procedure C/ (see above)

% -- Analysis conditions to contrast
cfgANA.EFFECT_ANALYSIS  = {'TOUCH','EMOTION','CONGRUENCE'};      % effect to analyse, either 'TOUCH', 'EMOTION', 'CONGRUENCE', 'VALENCE', or 'ROLE', or 'TOUCH_NoControl'
cfgANA.role_contrast1   = 'all';        % role to consider within contrast 1 ('all' or 'empathizer' or 'transmitter')
cfgANA.role_contrast2   = 'all';        % role to consider within contrast 2 ('all' or 'empathizer' or 'transmitter')

% -- Physiological signal analysis
cfgANA.physio_sig           = 'EDA';    % within procedure B/ -> perform analysis of 'PRV', 'RVT' or 'EDA'       
cfgANA.coupling_measure     = 'xcorr';  % within procedure B/ -> choice of coupling measure: 'xcorr' or 'PLV'
cfgANA.perform_jointPhysioQuest = 1;  	% within procedure B/ -> perform joint analysis physio signals and questionnaires data
cfgANA.do_PLV               = 0;        % apply Phase Locking Value on physio measure
cfgANA.PLV_winSize          = 4096;     % PLV: size of the windows in number of samples (@Fs=128Hz -> 1024=8s, 4096=32s)
cfgANA.PRV_logTransform     = 0;        % apply (cross-)correlation on log(-PP_interval) instead of raw values
cfgANA.do_xcorr_hann        = 0;        % apply a hanning weighting function before processing xcorr
cfgANA.derived_measures     = 1;        % analyse (time and/or freq) measures derived from physio signals
cfgANA.do_xcorr_interCouple = 1;        % control xcorr intra couple with xcorr computed across couples 
cfgANA.nBoot                = 10000;    % number of bootstrap sampling   
cfgANA.trialLim             = [30 120]; % put 'all' or trial [start end] in seconds, useful for removing trial begin and/or end.
cfgANA.xcorr_lag_type       = 'max';    % 'max' or 'mean', type of lag to be considered in cross-correlation. 'max': take only max xcorr value. 'mean': compute mean in following interval
cfgANA.xcorr_lag_interval   = [-5 5];   % cross-correlation lag to consider for 'mean' xcorr, either '0' or interval in sec (e.g. [-5 5])

% -- Analysis t-test parameters
cfgANA.nbPerms        	= 1024;         % number of permutation for t-max test
cfgANA.tTest_alpha      = 0.05;        	% t-test significance level
cfgANA.tTest_doAlphaCorr = 1;        	% perform correction (Bonferroni) of t-test significance level according to number of frequencies and IC tested
cfgANA.tTest_freqRange  = [2 30];       % Freq range for testing frequency    
cfgANA.tTest_freqRes    = 1;            % frequency resolution (in Hz) used for t-test
cfgANA.tTest_powType    = 'absolute';     % spectrum normalization: either 'absolute' or 'normal' or 'relative' 

% -- Data importation parameters
cfgANA.useSubjectiveReports = 1;        % import online self report on emotion intensity
cfgANA.keepEmotionTrial = [];           % empty or scalar (1 to 8), use online self report on emotion intensity to keep only X (on 8) more intense trials for each subject
cfgANA.manualArtifRej   = 0;            % '1' to consider manual artifact rejection from cartool marker file (i.e. remove these periods from .set files)
cfgANA.nbSmoothSamples  = 8;            % specifies nb of samples to smooth at begin and end of artifact-free periods
cfgEEG.timeOffset       = 30;           % Crop time data before this startpoint (in seconds), '0' for nothing
cfgANA.cropRestPeriods  = 1;            % Crop inter-trial rest periods to 30s (remove data in their middle)
cfgANA.standardizeEEG   = 'none';       % EEG standardization, either 'subject', 'frequency' or 'none'.
cfgANA.BPfilter         = []; %[.5 45];     	% Frequency band [fmin fmax] for band-pass filtering ([] for nothing) 

% -- EEG data parameters
cfgEEG.Fs               = 128 ;         % sampling rate (assumed identical for each subject)
cfgEEG.NbChannels       = 46;           % number of channels per subject (assumed identical for each subject)
cfgEEG.lengthTrial      = 120;          % Length of trial (condition) periods in seconds
cfgEEG.lengthRest       = 30;          % Length of inter-trial (rest) periods in seconds
cfgEEG.chan_str         = { 'FP1','FP2','AFZ','F7','F3','FZ','F4','F8','T7','C3','CZ','C4',...
                            'T8','P7','P3','PZ','P4','P8','O1','O2','Resp','Pulse','GSR',...
                            'FP1','FP2','AFZ','F7','F3','FZ','F4','F8','T7','C3','CZ','C4',...
                            'T8','P7','P3','PZ','P4','P8','O1','O2','Resp','Pulse','GSR'};
cfgEEG.EEG_channels     = [1:20,24:43]; % indexes of EEG channels
cfgEEG.EXT_channels     = [21:23,44:46];% indexes of external channels (Resp, Pulse, GSR)

% -- Configuration for cospectra processing
cfgEEG.FreqRes          = .5;           % Cospectra frequency resolution in Hz
cfgEEG.WinOverlap       = .5;           % windows overlap (possible values: 0->.99)
cfgEEG.WinType          = 'hanning' ;   % 'hanning' 'welch' ...
cfgEEG.FreqRangeCosp    = [0 48];       % Freq range for cospectra processing
cfgEEG.FreqRangeNorm    = [5 48];       % Freq range for spectrum normalization (mainly for visualization purpose)

% -- Whitening and BSS configuration
cfgBSS.ApplyWhitening   = 0;        % Perform pre-whitening of cospectra for a steady initialisation of JBSS  
cfgBSS.ApplyNonDiagFct 	= 1;        % Weight cospectra with non-diagonality function   
cfgBSS.PercentVariance  = 0.01;   	% Specify minimum variance explained to reduce data after whitening (0 for no data reduction)
cfgBSS.FreqRangeBSS  	= [5 28];   % Freq range for whitening and (j)BSS (usually [5 28] Hz)      
cfgBSS.jDiagContrast    = 1;      	% '1' for BSS to simultaneously diagonalize cospectra from both contrast conditions 
%                                           (concatenate them), '0' to take only first contrast
 
%% --- Various initializations
% RawDataPath contains path with couples' EEG raw .mat files, (root folder
% from where STUDY files are stored)
% cfgANA.RawDataPath = [ fileparts(cfgANA.DataPath(1:end-1)) '\']; 
cfgANA.RawDataPath = cfgANA.DataPath;

% Here we initialize couples ID based on experimental priors
if strcmp(cfgANA.couple_set ,'all')    
    cfgANA.coupleID = [2:5,7:16];              % keep all couple data (couple n°6 has no data)
elseif strcmp(cfgANA.couple_set ,'best')
    cfgANA.coupleID = [7 9 10 11 13 14 15 16]; % choice of best couples is based on subjective experimental assessment
elseif strcmp(cfgANA.couple_set ,'verybest')
    cfgANA.coupleID = [9 10 14 15 16];         % choice of best couples is based on subjective experimental assessment
elseif strcmp(cfgANA.couple_set ,'eyesopen')
    cfgANA.coupleID = [2 3 4 5 7];             % subjects with eyes open during the experiment
elseif strcmp(cfgANA.couple_set ,'eyesclosed')
    cfgANA.coupleID = [9 10 11 13 14 15 16];   % subjects with eyes open during the experiment
elseif strcmp(cfgANA.couple_set ,'worst')    
    cfgANA.coupleID = [2 3 4 5 8];        
end

% Data segmentation
cfgEEG.CondTrigVal      = [8:12,16:20]; % Trigger indexes used for all conditions in the experiment
cfgEEG.nbCondTot        = length(cfgEEG.CondTrigVal);
cfgEEG.subRoles(1,:)    = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32]; % this are subject index
cfgEEG.subRoles(2,:)  	= [0 1 0 1 1 0 1 0 1 0  0  1  0  1  1  0  0  1  0  1  1  0  0  1  0  1  1  0  0  1  0  1];  % subject roles: '0' for Empathizer, '1' for Transmitter 
couple_select           = sort([2*cfgANA.coupleID-1,2*cfgANA.coupleID]);
cfgEEG.subRoles         = cfgEEG.subRoles(:,couple_select);
cfgEEG.effect_labels    = { 'Touch_Congruent' ; 'Touch_noCongruent' ; 'Touch_Control' ; ...
                            'noTouch_Congruent'; 'noTouch_noCongruent'; 'noTouch_Control';...
                            'Touch_ValencePos' ; 'Touch_ValenceNeg' ; 'noTouch_ValencePos' ; 'noTouch_ValenceNeg'};                       
cfgEEG.effect_cond_ix = {[17 20] ; [18 19] ; 16; [9 12] ; [10 11] ; 8; [17 18]; [19 20] ; [9 10]; [11 12]};   % these run indexes correspond to 'effect_labels'
cfgEEG.cond_labels    = {   'noTouch-Neutral' ; 'noTouch-Congruent(+)' ; 'noTouch-noCongruent(+)' ; ...
                            'noTouch-noCongruent(-)'; 'noTouch-Congruent(-)'; 'Touch-Neutral' ; ...
                            'Touch-Congruent(+)'; 'Touch-noCongruent(+)'; 'Touch-noCongruent(-)'; 'Touch-Congruent(-)'};

if strcmp(cfgANA.trialLim,'all')
       cfgANA.trialLim = [0 120]; 
end
switch cfgANA.EFFECT_ANALYSIS{1}
    % 1. contrast TOUCH / no TOUCH
    case 'TOUCH'
        cfgANA.effect_contrast1     = {'Touch_Congruent';'Touch_noCongruent';'Touch_Control'};   
        cfgANA.effect_contrast2     = {'noTouch_Congruent';'noTouch_noCongruent';'noTouch_Control'};  
        cfgANA.effect_contrast1_str = 'touch';
        cfgANA.effect_contrast2_str = 'no touch';
    % 2. contrast EMOTION / CONTROL
    case 'EMOTION'
        cfgANA.effect_contrast1    = {'Touch_Congruent';'noTouch_Congruent';'Touch_noCongruent';'noTouch_noCongruent'};     
        cfgANA.effect_contrast2    = {'Touch_Control';'noTouch_Control'}; 
      	cfgANA.effect_contrast1_str = 'emotion';
        cfgANA.effect_contrast2_str = 'control';
    % 3. contrast CONGRUENT / no CONGRUENT
    case 'CONGRUENCE'
        cfgANA.effect_contrast1   = {'Touch_Congruent';'noTouch_Congruent'};      
        cfgANA.effect_contrast2   = {'Touch_noCongruent';'noTouch_noCongruent'}; 
      	cfgANA.effect_contrast1_str = 'congruent';
        cfgANA.effect_contrast2_str = 'no congruent';
    % 4. contrast EMOTION POS / EMOTION NEG
    case 'VALENCE'
        cfgANA.effect_contrast1   = {'Touch_ValencePos' ; 'noTouch_ValencePos'};      
        cfgANA.effect_contrast2   = {'Touch_ValenceNeg' ; 'noTouch_ValenceNeg'};  
       	cfgANA.effect_contrast1_str = 'emotion pos';
        cfgANA.effect_contrast2_str = 'emotion neg';
    % 5. contrast EMPATHIZER / TRANSMITTER (all conditions except CONTROL)
    case 'ROLE' 
        cfgANA.effect_contrast1   = {'Touch_Congruent';'Touch_noCongruent'; 'noTouch_Congruent';'noTouch_noCongruent'};     
        cfgANA.effect_contrast2   = {'Touch_Congruent';'Touch_noCongruent';'noTouch_Congruent';'noTouch_noCongruent'};
        cfgANA.effect_contrast1_str = 'empathizer';
        cfgANA.effect_contrast2_str = 'transmitter';
   % 6. contrast TOUCH / no TOUCH (WITHOUT CONTROL!!!)
    case 'TOUCH_NoControl'
        cfgANA.effect_contrast1     = {'Touch_Congruent';'Touch_noCongruent'};   
        cfgANA.effect_contrast2     = {'noTouch_Congruent';'noTouch_noCongruent'};  
        cfgANA.effect_contrast1_str = 'touch (no control)';
        cfgANA.effect_contrast2_str = 'no touch (no control)'; 
    % 7. contrast CONTROL TOUCH / CONTROL no TOUCH 
    case 'CONTROL'
        cfgANA.effect_contrast1     = {'Touch_Control'};   
        cfgANA.effect_contrast2     = {'noTouch_Control'};  
        cfgANA.effect_contrast1_str = 'control touch';
        cfgANA.effect_contrast2_str = 'control no touch'; 
    otherwise
        error('unrecognized effect type')
end
if ((cfgANA.derived_measures==1) && (strcmp(cfgANA.physio_sig,'PRV') == 1) && (cfgANA.trialLim(2)-cfgANA.trialLim(1))<60)
    error('You must select a bigger run length to perform analysis of PRV frequency measures.')
%     cfgANA.derived_measures = 0;
end
if (strcmp(cfgANA.effect_contrast1,'all') & cfgBSS.jDiagContrast==0)
     cfgANA.effect_contrast1 = cfgEEG.effect_labels;
end
cfgANA.cond_contrast1 = sort(cell2mat(reshape(cfgEEG.effect_cond_ix(ismember(cfgEEG.effect_labels,cfgANA.effect_contrast1)),1,[])));
cfgANA.cond_contrast2 = sort(cell2mat(reshape(cfgEEG.effect_cond_ix(ismember(cfgEEG.effect_labels,cfgANA.effect_contrast2)),1,[])));
cfgANA.nb_cond_contrast1 = length(cfgANA.cond_contrast1);
cfgANA.nb_cond_contrast2 = length(cfgANA.cond_contrast2);
if strcmp(cfgANA.role_contrast1,'all')
    cfgANA.sub_contrast1 = cfgEEG.subRoles(1,:);
elseif strcmp(cfgANA.role_contrast1,'empathizer')
    cfgANA.sub_contrast1 = cfgEEG.subRoles(1,~cfgEEG.subRoles(2,:));
elseif strcmp(cfgANA.role_contrast1,'transmitter')
    cfgANA.sub_contrast1 = cfgEEG.subRoles(1,~~cfgEEG.subRoles(2,:));   
end
if strcmp(cfgANA.role_contrast2,'all')
    cfgANA.sub_contrast2 = cfgEEG.subRoles(1,:);
elseif strcmp(cfgANA.role_contrast2,'empathizer')
    cfgANA.sub_contrast2 = cfgEEG.subRoles(1,~cfgEEG.subRoles(2,:));
elseif strcmp(cfgANA.role_contrast2,'transmitter')
    cfgANA.sub_contrast2 = cfgEEG.subRoles(1,~~cfgEEG.subRoles(2,:));   
end    

% Channel selection
if strcmp(cfgANA.dataType,'EEG_only')
    cfgEEG.chan_str = cfgEEG.chan_str(cfgEEG.EEG_channels);
    cfgEEG.chanSel = cfgEEG.EEG_channels;
elseif strcmp(cfgANA.dataType,'EXT_only')
    cfgEEG.chan_str = cfgEEG.chan_str(cfgEEG.EXT_channels);
    cfgEEG.chanSel = cfgEEG.EXT_channels;
end
cfgEEG.nbChan   = length(cfgEEG.chanSel);    	% total number of channels
cfgEEG.nbChanSub = cfgEEG.nbChan/2;          	% number of channels per subject
cfgEEG.S        = [cfgEEG.nbChanSub  ; cfgEEG.nbChanSub];   % size of each subject channel sets
cfgEEG.chanSet  = [1:cfgEEG.nbChanSub ; cfgEEG.nbChanSub+1:cfgEEG.nbChan];  % subject1 = 1st row, subject2 = 2nd row
cfgEEG.nbSub    = length(cfgEEG.S);             % number of subject per group (here couple so obviously = 2)
cfgEEG.nbCouple = length(cfgANA.coupleID);      % number of couples
cfgEEG.nbSubTot = cfgEEG.nbSub*cfgEEG.nbCouple; % total number of subjects


% read electrode coordinates
try
    cfgEEG.eloc = readlocs('EEGLAB_electrodes_Dual_Empathy.loc','filetype','loc'); 
catch
    error('You forgot to launch eeglab first! (again...)');  % (do not forget to launch EEGlab first)
end

% BP filter design
if ~isempty(cfgANA.BPfilter)
    [cfgEEG.H_bp, cfgEEG.b_bp, cfgEEG.a_bp]  = designBandPassFilter_FIR_strongCut2(cfgANA.BPfilter(1), cfgANA.BPfilter(2), cfgEEG.Fs); 
%     fvtool(cfgEEG.H_bp,'Fs',cfgEEG.Fs,'Analysis','freq','Legend','on') % show filter characterists (magnitude and phase plot)
end

% various init
if (cfgANA.xcorr_lag_interval==0)
	cfgANA.xcorr_lag_interval = [-1/cfgEEG.Fs 1/cfgEEG.Fs]; % setting minimum lag around 0 for instantaneous cross-correlation
end
fprintf('\n \n');


%% --- 0) Importing and segmenting raw EEG for each subject / condition

% here we use data that was imported first using 'Import_gTEC_data.m' 
for couple_ix = cfgANA.coupleID
    couple_ix_total = find(ismember(cfgANA.coupleID,couple_ix));
    % --- load data
    dataFile = ['couple' int2str(couple_ix) '_raw.mat'];
    dataPath = [cfgANA.RawDataPath 'Couple_' int2str(couple_ix) '\' ];
%     fprintf('\n');
%     disp(['Loading current datafile:    ' dataFile '...']);
    try
        load([dataPath dataFile]); 
    catch
        error(['INEXISTANT DATA FILE: ' dataFile '. Have you set the right path?...']); 
    end

    % --- keep only selected data type (EEG, EXT or both)
    EEG.data = EEG.data(cfgEEG.chanSel,:);

    % --- BP-filter data
    if ~isempty(cfgANA.BPfilter)
        disp('Filtering EEG data...');  
        EEG.data_filt = filtfilt(cfgEEG.b_bp, cfgEEG.a_bp, EEG.data')';  % zero-phase filtering
        disp('EEG data filtered');  
    end
    
    % --- load markers for manual artifact rejection
    if cfgANA.manualArtifRej == 1
        mrkFile = [dataPath 'couple' int2str(couple_ix) '_filtered.ep.MANUAL-ARTIF-REJ.mrk'];
        if ~exist(mrkFile,'file')
            error('Marker file for manual artifact rejection does not exist!');
        end
        EEG.mask_mrk  = import_cartool_mrk(mrkFile, size(EEG.data,2));
%         disp('Marker file loaded for manual artifact rejection.');
    else
        EEG.mask_mrk  = true(1, size(EEG.data,2));    % select everything
%         disp('No manual artifact rejection (everything selected).')
    end
    
    % --- format trigger to specified conditions and trial begin/end
    % - EEG.trigVal(:,1) -> trigger name (condition id)
    % - EEG.trigVal(:,2) -> start sample 
    % - EEG.trigVal(:,3) -> end sample
    % - EEG.trigVal(:,4) -> condition order during the experiment
    % - EEG.trigIndex    -> cell array for indexing 'EEG.data' to condition
    %                       periods, {nb conditions} ( [1,trialLength] , ID )
    % - EEG.maskedTrigIndex -> idem with mask set using manual artifact rej
    % - EEG.mask_mrk_inCond -> mask with artifact free in condition periods
    EEG.trigVal = EEG.trigVal(ismember(EEG.trigVal(:,1),[cfgEEG.CondTrigVal 7]),:); % remove sometime spurious trigger at begin
    EEG.trigVal = EEG.trigVal(~ismember(EEG.trigVal(:,2),1),:);
    
    % -- Crop inter-trial rest periods to 30s (remove data in their middle)
    if cfgANA.cropRestPeriods == 1  
        begin_rest_samples  = EEG.trigVal(ismember(EEG.trigVal(:,1),7),2);
        end_rest_samples    = EEG.trigVal(ismember(EEG.trigVal(:,1),cfgEEG.CondTrigVal),2);
        end_rest_samples    = [ end_rest_samples(2:end) ; size(EEG.data,2) ]; 
        middle_rest_samples = begin_rest_samples + 15*cfgEEG.Fs;    % crop from 15s after begin of rest period
        uper_crop_samples   = middle_rest_samples + (end_rest_samples-begin_rest_samples-30*cfgEEG.Fs); % keep 30s only
        mask_crop_rest      = true(1,size(EEG.data,2));
        for rest_ix = 1:length(begin_rest_samples)
            mask_crop_rest(middle_rest_samples(rest_ix):uper_crop_samples(rest_ix)) = false;
        end
        % smoothing data at edges removed periods
        EEG.data = EEG.data .* (ones(cfgEEG.nbChan,1) * smoothMask(mask_crop_rest,cfgANA.nbSmoothSamples,2)); 
        % crop data and adjust triggers
        EEG.data = EEG.data(:,mask_crop_rest);
        cumulated_delays = sort(repmat(cumsum(uper_crop_samples-middle_rest_samples),[2,1]));
        EEG.trigVal(3:end,2) = EEG.trigVal(3:end,2) - cumulated_delays(1:end-2); % only adjust from 2nd condition period to end
        EEG.mask_mrk = EEG.mask_mrk(mask_crop_rest);
        clear begin_rest_samples end_rest_samples middle_rest_samples uper_crop_samples mask_crop_rest cumulated_delays rest_ix
    end
    
    % -- Crop data at begin if necessary (keep only 'timeOffset' seconds before first condition onset)
    if cfgEEG.timeOffset > 0
        crop_beg_sample = EEG.trigVal(1,2) - cfgEEG.Fs*cfgEEG.timeOffset; 
        length_experiment = cfgEEG.Fs*cfgEEG.nbCondTot*(cfgEEG.lengthTrial+cfgEEG.lengthRest)+cfgEEG.Fs*10; % keep 10s after last period
        crop_end_sample = crop_beg_sample + length_experiment -1; % remove some samples at the end, to obtain exact same data length for all subjects
        if(crop_beg_sample <= 0)
           error('Parameter ''cfgEEG.timeOffset'' is too high (not enough data at begin)'); 
        end
        if (crop_end_sample > size(EEG.data,2))
           error('Impossible to keep 10s at experiment end: not enough data'); 
        end
        EEG.data        = EEG.data(:,crop_beg_sample:crop_end_sample);
        % adjust masks and triggers after cropping data
        EEG.mask_mrk    = EEG.mask_mrk(crop_beg_sample:crop_end_sample);
        EEG.trigVal(:,2)= EEG.trigVal(:,2) - crop_beg_sample;
        clear crop_beg_sample crop_end_sample length_experiment
    end
    
    % -- create trigger related structures
    startOffset         = cfgANA.trialLim(1) * cfgEEG.Fs ;
    endOffset           = cfgANA.trialLim(2) * cfgEEG.Fs ;
    cfgANA.trial_sample = cfgEEG.Fs*(cfgANA.trialLim(2)-cfgANA.trialLim(1)) ;  % data length per trial in samples
    EEG.trigVal         = EEG.trigVal(ismember(EEG.trigVal(:,1), cfgEEG.CondTrigVal),:);  % this rejects inter trial rest periods
    EEG.trigVal(:,2)    = EEG.trigVal(:,2) + startOffset;
    EEG.trigVal(:,3)    = EEG.trigVal(:,2) + endOffset - startOffset -1;
    [temp condOrder]    = sort(EEG.trigVal(:,1),'ascend');
    EEG.trigVal         = EEG.trigVal(condOrder,:);
    EEG.trigVal(:,4)    = condOrder;
    EEG.trigIndex     	= zeros(cfgEEG.nbCondTot, cfgANA.trial_sample); 
    EEG.trigIndex    	= num2cell(EEG.trigIndex,2);
    EEG.maskedTrigIndex = zeros(cfgEEG.nbCondTot, cfgANA.trial_sample); 
    EEG.maskedTrigIndex = num2cell(EEG.trigIndex,2);
    EEG.mask_mrk_inCond = false(1,size(EEG.data,2));
    for cond_ix = 1:cfgEEG.nbCondTot
        EEG.trigIndex{cond_ix}      = EEG.trigVal(cond_ix,2):EEG.trigVal(cond_ix,3); % EEG.trigVal(cond_ix,2) + (EEG.trigVal(cond_ix,2):EEG.trigVal(cond_ix,3));
        EEG.trigIndex{cond_ix,2}    = EEG.trigVal(cond_ix,1) ;      % condition ID
        EEG.maskedTrigIndex{cond_ix}= EEG.trigVal(cond_ix,2) + ...
                                            find( EEG.mask_mrk(EEG.trigVal(cond_ix,2):EEG.trigVal(cond_ix,3)) );
        EEG.maskedTrigIndex{cond_ix,2} = EEG.trigVal(cond_ix,1) ;   % condition ID     
        EEG.mask_mrk_inCond(EEG.maskedTrigIndex{cond_ix}) = true;
    end
    % weight EEG with smoothing windows to remove time discontinuities when selecting artifact periods
    if (cfgANA.manualArtifRej == 1) 
        EEG.mask_mrk_smooth = smoothMask(EEG.mask_mrk,cfgANA.nbSmoothSamples,2);	
        EEG.data = EEG.data .* (ones(cfgEEG.nbChan,1) * EEG.mask_mrk_smooth); % smoothing EEG at edges of artifact periods
    end
    
    % --- create a cell array for all trigger related structures
    for sub_ix = 1:cfgEEG.nbSub  
      	sub_ix_total = cfgEEG.nbSub*(couple_ix_total-1)+sub_ix; 
        TRIG.trigSub(sub_ix_total).trigVal = EEG.trigVal;
        TRIG.trigSub(sub_ix_total).dispTrig = cell(2*cfgEEG.nbCondTot,2);
        TRIG.trigSub(sub_ix_total).dispTrig2 = cell(cfgEEG.nbCondTot,1);
        TRIG.trigSub(sub_ix_total).restVal  = zeros(cfgEEG.nbCondTot+1,2);
        for cond_ix = 1:cfgEEG.nbCondTot
            TRIG.trigSub(sub_ix_total).dispTrig{2*cond_ix-1,1}  = EEG.trigVal(cond_ix,2);
            TRIG.trigSub(sub_ix_total).dispTrig{2*cond_ix-1,2}  = int2str(EEG.trigVal(cond_ix,1)); 
            TRIG.trigSub(sub_ix_total).dispTrig{2*cond_ix,1}    = EEG.trigVal(cond_ix,3);
            TRIG.trigSub(sub_ix_total).dispTrig{2*cond_ix,2}    = int2str(EEG.trigVal(cond_ix,1));
            TRIG.trigSub(sub_ix_total).dispTrig2{cond_ix}       = cfgEEG.cond_labels{EEG.trigVal(cond_ix,4)};  
        end
        [temp sort_ix] = sort(EEG.trigVal(:,4));
        rest_beg = [0 ; EEG.trigVal(sort_ix,3)+(cfgEEG.lengthTrial-cfgANA.trialLim(2))*cfgEEG.Fs];
        rest_end = [EEG.trigVal(sort_ix,2)-cfgANA.trialLim(1)*cfgEEG.Fs; size(EEG.data,2)];
        TRIG.trigSub(sub_ix_total).restVal = [rest_beg rest_end]; % this contains begin and end samples of rest periods
    end
    clear temp sort_ix rest_beg rest_end
    
    % --- load subject answer to various questionnaires
    if (cfgANA.procedureC_behav==1) || (cfgANA.useSubjectiveReports==1)
        for sub_ix = 1:cfgEEG.nbSub  
            sub_ix_total    = cfgEEG.nbSub*(couple_ix_total-1)+sub_ix; 
            sub_ix_noEmpty  = cfgEEG.nbSub*(couple_ix-1)+sub_ix; 
            % load subject answer empathy questionnaire
            dataFile_SR = [cfgANA.RawDataPath 'QUESTIONNAIRES\Subject_' int2str(sub_ix_total) '_Answers_Quest_Empathy.txt'];
            rawEmpathyScores = readSelfReport(dataFile_SR,1); 
            QUEST.score_empathy{sub_ix_total} = cmpEmpathyFavre(rawEmpathyScores);
            % load subject answer couple's harmony
            dataFile_SR = [cfgANA.RawDataPath 'QUESTIONNAIRES\Subject_' int2str(sub_ix_noEmpty) '_Answers_Quest_Harmony.txt'];
            QUEST.score_harmony{sub_ix_total} = readSelfReport(dataFile_SR,1); % colums: questions, rows: trials
            % load subject answer to online self reports on emotion intensity and perceived coupling
            dataFile_SR = [cfgANA.RawDataPath 'QUESTIONNAIRES\Subject_' int2str(sub_ix_noEmpty) '_Answers_Quest_Self_Report.txt'];
            QUEST.score_onlineEmo{sub_ix_total} = readSelfReport(dataFile_SR,2); % colums: questions, rows: trials
        end
        clear sub_ix_noEmpty sub_ix rawEmpathyScores dataFile_SR
    end
    
    % -- keep only X (over 8 no-control trials) more intense trials for each subject 
    if cfgANA.useSubjectiveReports == 1
        % -- ... and keep preferably the same trials (that's the heavy part)
        if ~isempty(cfgANA.keepEmotionTrial)
            [sortValSub1 ixSub1] = sort(QUEST.score_onlineEmo{1,sub_ix_total-1}(:,1),'descend');
            [sortValSub2 ixSub2] = sort(QUEST.score_onlineEmo{1,sub_ix_total}(:,1),'descend');
            % remove control condition (not pertinent for emotion intensity)
            ixSub1 = ixSub1(~ismember(ixSub1,EEG.trigVal([1 6],4))); sortValSub1 = sortValSub1(~ismember(ixSub1,EEG.trigVal([1 6],4)));
            ixSub2 = ixSub2(~ismember(ixSub2,EEG.trigVal([1 6],4))); sortValSub2 = sortValSub2(~ismember(ixSub2,EEG.trigVal([1 6],4)));
            % keep index of (X+?) conditions with highest scores 
            ixSub1 = ixSub1(sortValSub1>=sortValSub1(cfgANA.keepEmotionTrial));
            ixSub2 = ixSub2(sortValSub2>=sortValSub2(cfgANA.keepEmotionTrial));
            % when possible place the same trials first
            lastSameVal1 = ixSub1(sortValSub1==sortValSub1(cfgANA.keepEmotionTrial));
            lastSameVal2 = ixSub2(sortValSub2==sortValSub2(cfgANA.keepEmotionTrial));
            ixSub1 = ixSub1(sortValSub1>sortValSub1(cfgANA.keepEmotionTrial));
            ixSub2 = ixSub2(sortValSub2>sortValSub2(cfgANA.keepEmotionTrial));
            lastSameVal1 = [lastSameVal1(ismember(lastSameVal1,lastSameVal2)) ; lastSameVal1(~ismember(lastSameVal1,lastSameVal2))];
            lastSameVal2 = [lastSameVal2(ismember(lastSameVal2,lastSameVal1)) ; lastSameVal2(~ismember(lastSameVal2,lastSameVal1))];
            % CONSTRUCT FINAL VECTOR OF X CONDITIONS TO KEEP  
            ixSub1 = [ixSub1 ; lastSameVal1(1:cfgANA.keepEmotionTrial-length(ixSub1))];
            ixSub2 = [ixSub2 ; lastSameVal2(1:cfgANA.keepEmotionTrial-length(ixSub2))];
            cfgEEG.SR_select_sub{sub_ix_total-1}= sort([cfgEEG.CondTrigVal(ixSub1) 8 16]); % indexes of most intense conditions + control conditions
            cfgEEG.SR_select_sub{sub_ix_total}  = sort([cfgEEG.CondTrigVal(ixSub2) 8 16]);
            clear sortValSub1 sortValSub2 ixSub1 ixSub2 lastSameVal1 lastSameVal2
        else
            cfgEEG.SR_select_sub{sub_ix_total-1}= cfgEEG.CondTrigVal;
            cfgEEG.SR_select_sub{sub_ix_total}  = cfgEEG.CondTrigVal;
        end
    end
        
    % --- data standardization INSIDE COUPLE using each subject total variance from all chans
    if strcmp(cfgANA.standardizeEEG,'subject')
        % here we standardize only with EEG inside experimental conditions (less noise with high amplitude)
        EEG.chanSTD = std(double(EEG.data(:,EEG.mask_mrk_inCond)),0,2);     % double: less round up error when (single) data is very long
        meanPow(couple_ix_total) = mean(EEG.chanSTD);                       % mean total power for both subject
        meanPow1   = mean(EEG.chanSTD(cfgEEG.chanSet(1,:)));    % mean variance for subject 1
        meanPow2   = mean(EEG.chanSTD(cfgEEG.chanSet(2,:)));    % mean variance for subject 2
        EEG.data(cfgEEG.chanSet(1,:),:) = EEG.data(cfgEEG.chanSet(1,:),:) ./ ((meanPow1/meanPow(couple_ix_total))*ones(cfgEEG.S(1),size(EEG.data,2))); 
        EEG.data(cfgEEG.chanSet(2,:),:) = EEG.data(cfgEEG.chanSet(2,:),:) ./ ((meanPow2/meanPow(couple_ix_total))*ones(cfgEEG.S(2),size(EEG.data,2)));
        disp('EEG standardized with each subject average power at all electrodes.');
        disp(['mean power before std: sub1= ' int2str(meanPow1) 'µV², sub2=' int2str(meanPow2) 'µV²']);
    end
    
    % --- perform EEG data segmentation to condition periods
    for sub_ix = 1:cfgEEG.nbSub  
        sub_ix_total = cfgEEG.nbSub*(couple_ix_total-1)+sub_ix; % this index increments for each new subject   
        for cond_ix = 1:cfgEEG.nbCondTot
          	mask_cond = EEG.maskedTrigIndex{cond_ix,1};
            EEG_sub(sub_ix_total).(['cond_' int2str(EEG.maskedTrigIndex{cond_ix,2})]) = EEG.data(cfgEEG.chanSet(sub_ix,:),mask_cond);	% store EEG for this specific subject and condition
        end   
        EEG_sub(sub_ix_total).('all') = EEG.data(cfgEEG.chanSet(sub_ix,:),:); % if not used, comment this to save memory 
        TRIG.trigSub(sub_ix_total).trigIndex = EEG.trigIndex; % this is to save trigIndex for each subject
    end   
    cond_SRreport_index = ismember(cell2mat(EEG.maskedTrigIndex(:,2)),cfgEEG.SR_select_sub{sub_ix_total}); % conditions with more intense trials
    EEG.lengthOK = length(cell2mat(EEG.maskedTrigIndex(cond_SRreport_index,1)')); % in samples
    cfgEEG.lengthExpe = size(EEG.data,2); % total data kept after segmentation (in samples)
    cfgEEG.t_axis_allExpe = (0:cfgEEG.lengthExpe-1) ./ cfgEEG.Fs; % time axis in seconds
    cfgEEG.lengthCond = length(mask_cond); % total data kept after segmentation (in samples)
%     disp([  'Total data length kept for analysis: ' int2str(EEG.lengthOK/cfgEEG.Fs) ...
%             ' seconds (' int2str(100*EEG.lengthOK/(cfgEEG.nbCondTot*cfgANA.trial_sample)) '%)']);
    disp(['Loaded:    ' dataFile ' Data kept for analysis: ' int2str(EEG.lengthOK/cfgEEG.Fs) 's (' int2str(100*EEG.lengthOK/(cfgEEG.nbCondTot*cfgANA.trial_sample)) '%)']);
end
% create a wheighting window with smooth hanning function at extremities only.
if (cfgANA.do_xcorr_hann == 1) 
    hann_length = 4096; % length of data to smooth at begin and end of time window before xcorr
	hann_win = hann(hann_length); 
    cfgANA.xcorr_wheighting_win = ones(1,cfgANA.trial_sample); 
    cfgANA.xcorr_wheighting_win([1:hann_length/2,end-hann_length/2:end]) = hann_win([1:hann_length/2,end-hann_length/2:end]);
end
% concatenate trigger indexes for each couple
TRIG.trigAll_couples    = zeros(cfgEEG.nbCouple, cfgANA.trial_sample, cfgEEG.nbCondTot);   % conditions are ordered such as: {8;9;10;11;12;16;17;18;19;20}
TRIG.trigAll_sub        = zeros(cfgEEG.nbSubTot, cfgANA.trial_sample, cfgEEG.nbCondTot);    % conditions are ordered such as: {8;9;10;11;12;16;17;18;19;20}
for couple_ix = 1:cfgEEG.nbCouple
    TRIG.trigAll_couples(couple_ix,:,:) = cell2mat(TRIG.trigSub(couple_ix*2-1).trigIndex(:,1))';
    TRIG.trigAll_sub(couple_ix*2-1,:,:) = TRIG.trigAll_couples(couple_ix,:,:); 
    TRIG.trigAll_sub(couple_ix*2,:,:) = TRIG.trigAll_couples(couple_ix,:,:);
end
% TRIG.trigAll_sub = reshape(repmat(TRIG.trigAll_couples',2,1,1),cfgEEG.nbCondTot,cfgEEG.nbSubTot)'; % format: (nbSubTot,nbTotalCond)
TRIG.run_begin_samples_couple   = squeeze(TRIG.trigAll_couples(:,1,:))-cfgANA.trialLim(1)*cfgEEG.Fs; % format: (nbCouples,nbTotalCond)
TRIG.run_begin_samples_sub      = reshape(repmat(TRIG.run_begin_samples_couple',2,1),cfgEEG.nbCondTot,cfgEEG.nbSubTot)'; % format: (nbSubTot,nbTotalCond)

% clean up data structures
clear EEG dataFile couple_ix startOffset endOffset setName subName condOrder temp cond_effect_index cond_SRreport_index final_cond_index
clear mrkFile tempTrigVal cond_ix couple_ix dataFile dataPath effect_ix mask_effect couple_ix_total couple_select mask_cond
clear sub_ix outFileName outFilePath cfgEEG.effect_labels effect_cond_ix sub_ix_total meanPow1 meanPow2 hann_win hann_length
EEG.data_sub    = EEG_sub;
if strcmp(cfgANA.standardizeEEG,'subject')
    EEG.meanPow = meanPow;
end
clear EEG_sub EEG_couple meanPow


%% % *** C/ Analysis of questionnaires   
if cfgANA.procedureC_behav == 1
    disp('Analysis of behavorial data.'),
    
    disp_harmony            = 1;
    disp_selfReports        = 1;
    disp_empathy            = 1;
    disp_corrSR             = 1;
    disp_corrEmpathySR      = 1;
    cfgANA.questAnalysis_effect = {'TOUCH','EMOTION','CONGRUENCE','VALENCE'};%,'ROLE'}
    %{'TOUCH','EMOTION','CONGRUENCE','CONGRUENCE-TOUCH','CONGRUENCE-NoTOUCH'};% {'TOUCH','EMOTION','CONGRUENCE','VALENCE'};%,'ROLE'};    
    cfgANA.questAnalysis_effect_str =  {'NO TOUCH  ','TOUCH  ','NEUTRAL  ','EMOTION  ','INCONGRUENT  ','CONGRUENT  ','NEGATIVE  ','POSITIVE  '};
    %{ 'NO TOUCH  ','TOUCH  ','NEUTRAL  ','EMOTION  ','INCONGRUENT  ','CONGRUENT  ',...
    % 'INCONGRUENT-TOUCH  ','CONGRUENT-TOUCH  ','INCONGRUENT-NOTOUCH  ','CONGRUENT-NOTOUCH  '};% {'NO TOUCH  ','TOUCH  ','NEUTRAL  ','EMOTION  ','INCONGRUENT  ','CONGRUENT  ','NEGATIVE  ','POSITIVE  '};
    
    % ** extract empathy scores from each subject
    for sub_ix = 1:cfgEEG.nbSubTot
        QUEST.scoreEmpathy.grouped_effect(1,sub_ix)	= QUEST.score_empathy{sub_ix}.empathy;
        QUEST.scoreEmpathy.grouped_effect(2,sub_ix)	= QUEST.score_empathy{sub_ix}.contagion;
        QUEST.scoreEmpathy.grouped_effect(3,sub_ix)	= QUEST.score_empathy{sub_ix}.emoCut;
        QUEST.scoreEmpathy.grouped_effect(4,sub_ix) = QUEST.score_empathy{sub_ix}.global;
    end
    
	% ** extract self-report scores from each effect / subject
    QUEST.scoreEmotion.grouped_effect   = zeros(length(cfgANA.questAnalysis_effect),2,cfgEEG.nbSubTot); % (nb effects, 2 contrasts, nb subjects) 
    QUEST.scoreCoupling.grouped_effect  = zeros(length(cfgANA.questAnalysis_effect),2,cfgEEG.nbSubTot); 
    QUEST.scoreEmotion.all_cond         = zeros(10,cfgEEG.nbSubTot);        
    QUEST.scoreCoupling.all_cond       	= zeros(10,cfgEEG.nbSubTot);   
    for sub_ix = 1:cfgEEG.nbSubTot
        QUEST.scoreEmotion.all_cond(:,sub_ix)  	= QUEST.score_onlineEmo{1,sub_ix}(:,1);
      	QUEST.scoreCoupling.all_cond(:,sub_ix)	= QUEST.score_onlineEmo{1,sub_ix}(:,2);
    end
    for effect_ix = 1:length(cfgANA.questAnalysis_effect) 
        effect = cfgANA.questAnalysis_effect{effect_ix};
        switch effect
        % 1. contrast TOUCH / no TOUCH
        case 'TOUCH'
            effect_contrast1     = {'noTouch_Congruent';'noTouch_noCongruent';'noTouch_Control'};
            effect_contrast2     = {'Touch_Congruent';'Touch_noCongruent';'Touch_Control'};   
            effect_contrast2_str = 'touch'; effect_contrast1_str = 'noTouch';
        % 2. contrast EMOTION / CONTROL
        case 'EMOTION'
            effect_contrast1    = {'Touch_Control';'noTouch_Control'}; 
            effect_contrast2    = {'Touch_Congruent';'noTouch_Congruent';'Touch_noCongruent';'noTouch_noCongruent'};     
            effect_contrast2_str = 'emotion'; effect_contrast1_str = 'control';
        % 3. contrast CONGRUENT / no CONGRUENT
        case 'CONGRUENCE'
            effect_contrast1   = {'Touch_noCongruent';'noTouch_noCongruent'}; 
            effect_contrast2   = {'Touch_Congruent';'noTouch_Congruent'};      
            effect_contrast2_str = 'congruent'; effect_contrast1_str = 'noCongruent';
        % 4. contrast EMOTION POS / EMOTION NEG
        case 'VALENCE'
            effect_contrast1   = {'Touch_ValenceNeg' ; 'noTouch_ValenceNeg'};  
            effect_contrast2   = {'Touch_ValencePos' ; 'noTouch_ValencePos'};      
            effect_contrast2_str = 'emotionPos'; effect_contrast1_str = 'emotionNeg';
        % 5. contrast EMPATHIZER / TRANSMITTER (all conditions)
        case 'ROLE' 
            effect_contrast1   = {'Touch_Control';'Touch_Congruent';'Touch_noCongruent';'noTouch_Control';'noTouch_Congruent';'noTouch_noCongruent'};     
            effect_contrast2   = {'Touch_Control';'Touch_Congruent';'Touch_noCongruent';'noTouch_Control';'noTouch_Congruent';'noTouch_noCongruent'}; 
            effect_contrast1_str = 'empathizer'; effect_contrast2_str = 'transmitter';
       % 6. contrast TOUCH / no TOUCH (WITHOUT CONTROL!!!)
        case 'TOUCH_NoControl'
            effect_contrast1     = {'noTouch_Congruent';'noTouch_noCongruent'};  
            effect_contrast2     = {'Touch_Congruent';'Touch_noCongruent'};   
            effect_contrast2_str = 'touch_noControl'; effect_contrast1_str = 'noTouch_noControl'; 
        % 7. contrast CONTROL TOUCH / CONTROL no TOUCH 
        case 'CONTROL'
            effect_contrast1     = {'noTouch_Control'};  
            effect_contrast2     = {'Touch_Control'};   
            effect_contrast2_str = 'controlTouch'; effect_contrast1_str = 'controlNoTouch';     
        % 8. contrast (CONGRUENT+TOUCH) / (NoCONGRUENT+TOUCH)
        case 'CONGRUENCE-TOUCH'
            effect_contrast1   = {'Touch_noCongruent'}; 
            effect_contrast2   = {'Touch_Congruent'};      
            effect_contrast2_str = 'congruentTouch'; effect_contrast1_str = 'noCongruentTouch';
    	% 9. contrast (CONGRUENT+NoTOUCH) / (NoCONGRUENT+NoTOUCH)
        case 'CONGRUENCE-NoTOUCH'
            effect_contrast1   = {'noTouch_noCongruent'}; 
            effect_contrast2   = {'noTouch_Congruent'};      
            effect_contrast2_str = 'congruentNoTouch'; effect_contrast1_str = 'noCongruentNoTouch';
        otherwise
            error('unrecognized effect type')
        end
        cond_contrast1 = sort(cell2mat(reshape(cfgEEG.effect_cond_ix(ismember(cfgEEG.effect_labels,effect_contrast1)),1,[])));
        cond_contrast2 = sort(cell2mat(reshape(cfgEEG.effect_cond_ix(ismember(cfgEEG.effect_labels,effect_contrast2)),1,[])));
        nb_cond_contrast1 = length(cond_contrast1);
        nb_cond_contrast2 = length(cond_contrast2);         
        QUEST.scoreEmotion.(effect_contrast1_str)   = zeros(nb_cond_contrast1,cfgEEG.nbSubTot);
        QUEST.scoreEmotion.(effect_contrast2_str)	= zeros(nb_cond_contrast2,cfgEEG.nbSubTot);
        QUEST.scoreCoupling.(effect_contrast1_str)	= zeros(nb_cond_contrast1,cfgEEG.nbSubTot);
        QUEST.scoreCoupling.(effect_contrast2_str)	= zeros(nb_cond_contrast1,cfgEEG.nbSubTot);
        empathizer_ix = 1;  transmitter_ix = 1;
        for sub_ix = 1:cfgEEG.nbSubTot
            if strcmp(effect,'ROLE')  % handle the special case where we test subject role: all 10 conditions but only nbSubTot/2 values
               if cfgEEG.subRoles(2,sub_ix) == 0
                   QUEST.scoreEmotion.(effect_contrast1_str)(:,[empathizer_ix,2*empathizer_ix])  = [QUEST.score_onlineEmo{1,sub_ix}(:,1), QUEST.score_onlineEmo{1,sub_ix}(:,1)]; % trick here: we double the data so it has the same size as in other conditions
                   QUEST.scoreCoupling.(effect_contrast1_str)(:,[empathizer_ix,2*empathizer_ix]) = [QUEST.score_onlineEmo{1,sub_ix}(:,2), QUEST.score_onlineEmo{1,sub_ix}(:,2)]; % (this does not change subsequent correlations)
                   empathizer_ix = empathizer_ix+1;
               else
                   QUEST.scoreEmotion.(effect_contrast2_str)(:,[transmitter_ix,2*transmitter_ix])  = [QUEST.score_onlineEmo{1,sub_ix}(:,1),QUEST.score_onlineEmo{1,sub_ix}(:,1)];
                   QUEST.scoreCoupling.(effect_contrast2_str)(:,[transmitter_ix,2*transmitter_ix]) = [QUEST.score_onlineEmo{1,sub_ix}(:,2),QUEST.score_onlineEmo{1,sub_ix}(:,2)];
                   transmitter_ix = transmitter_ix+1;
               end
            else    % for all other cases
                for cond_ix = 1:nb_cond_contrast1
                    curr_cond = TRIG.trigSub(1,sub_ix).trigVal(TRIG.trigSub(1,sub_ix).trigVal(:,1)==cond_contrast1(cond_ix),4);
                    QUEST.scoreEmotion.(effect_contrast1_str)(cond_ix,sub_ix) = QUEST.score_onlineEmo{1,sub_ix}(curr_cond,1);
                    QUEST.scoreCoupling.(effect_contrast1_str)(cond_ix,sub_ix) = QUEST.score_onlineEmo{1,sub_ix}(curr_cond,2);
                end
                for cond_ix = 1:nb_cond_contrast2
                    curr_cond = TRIG.trigSub(1,sub_ix).trigVal(TRIG.trigSub(1,sub_ix).trigVal(:,1)==cond_contrast2(cond_ix),4);
                    QUEST.scoreEmotion.(effect_contrast2_str)(cond_ix,sub_ix) = QUEST.score_onlineEmo{1,sub_ix}(curr_cond,1);
                    QUEST.scoreCoupling.(effect_contrast2_str)(cond_ix,sub_ix) = QUEST.score_onlineEmo{1,sub_ix}(curr_cond,2);
                end
            end
        end
        QUEST.scoreEmotion.(effect_contrast1_str)       = mean(QUEST.scoreEmotion.(effect_contrast1_str),1);     % averaging over conditions within effect
        QUEST.scoreEmotion.(effect_contrast2_str)       = mean(QUEST.scoreEmotion.(effect_contrast2_str),1);     % averaging over conditions within effect
        QUEST.scoreCoupling.(effect_contrast1_str)      = mean(QUEST.scoreCoupling.(effect_contrast1_str),1);    % averaging over conditions within effect
        QUEST.scoreCoupling.(effect_contrast2_str)      = mean(QUEST.scoreCoupling.(effect_contrast2_str),1);    % averaging over conditions within effect
        QUEST.scoreEmotion.grouped_effect(effect_ix,1,:) 	= QUEST.scoreEmotion.(effect_contrast1_str);
        QUEST.scoreEmotion.grouped_effect(effect_ix,2,:) 	= QUEST.scoreEmotion.(effect_contrast2_str);
        QUEST.scoreCoupling.grouped_effect(effect_ix,1,:) 	= QUEST.scoreCoupling.(effect_contrast1_str);
        QUEST.scoreCoupling.grouped_effect(effect_ix,2,:) 	= QUEST.scoreCoupling.(effect_contrast2_str);
    end
	clear cond_contrast1 cond_contrast2 cond_ix curr_cond curr_sub_ix effect_contrast1 effect_contrast1_str sub_ix
    clear effect_contrast2 effect_contrast2_str effect_ix mask_cond nb_cond_contrast1 nb_cond_contrast2 effect empathizer_ix transmitter_ix

    
    % *** Compute significance levels from online self reports in each condition
    QUEST.scoreEmotion.paired_test_p = zeros(1,length(cfgANA.questAnalysis_effect));
    QUEST.scoreCoupling.paired_test_p = zeros(1,length(cfgANA.questAnalysis_effect));
    for effect_ix = 1:length(cfgANA.questAnalysis_effect) 
        vec1 = squeeze(QUEST.scoreEmotion.grouped_effect(effect_ix,1,:));
        vec2 = squeeze(QUEST.scoreEmotion.grouped_effect(effect_ix,2,:));
        QUEST.scoreEmotion.paired_test_p(effect_ix) = signrank(vec1,vec2); % Wilcoxon signed rank test for paired difference test on ordinal data
        vec1 = squeeze(QUEST.scoreCoupling.grouped_effect(effect_ix,1,:));
        vec2 = squeeze(QUEST.scoreCoupling.grouped_effect(effect_ix,2,:));
        QUEST.scoreCoupling.paired_test_p(effect_ix) = signrank(vec1,vec2);
    end
    clear effect_ix vec1 vec2
    
    
  	% *** Display score from online self reports on emotion intensity and perceived coupling    
    if (disp_selfReports == 1)
        xtick = sort([(1:length(cfgANA.questAnalysis_effect))-.2,(1:length(cfgANA.questAnalysis_effect))+.2]);
        % emotion intensity
        emotion_mean_display  	= mean(QUEST.scoreEmotion.grouped_effect,3);
        emotion_sd_display      = std(QUEST.scoreEmotion.grouped_effect,0,3);
        figure('Name',['Emotion intensity scores, couples: ' cfgANA.couple_set],'WindowStyle','docked','NumberTitle','off'),
        barPlot(emotion_mean_display,QUEST.scoreEmotion.paired_test_p,[.05 .01 .001],9), hold on, 
        errorb(emotion_mean_display,emotion_sd_display,'top');
        ylim([0 12]), title('Emotion intensity'),
        set(gca,'Xtick',xtick,'XTickLabel',cfgANA.questAnalysis_effect_str,'fontsize',8), 
        rotateticklabel(gca,50);  set(gca,'OuterPosition',[0 0.3 1 .6]), 
        ylabel('Emotion intensity during previous run'), colormap summer,
        % perceived coupling
        coupling_mean_display  	= mean(QUEST.scoreCoupling.grouped_effect,3);
        coupling_sd_display     = std(QUEST.scoreCoupling.grouped_effect,0,3);
        figure('Name',['Perceived coupling, couples: ' cfgANA.couple_set],'WindowStyle','docked','NumberTitle','off'),
        barPlot(coupling_mean_display,QUEST.scoreCoupling.paired_test_p,[.05 .01 .001],9), hold on, 
        errorb(coupling_mean_display,coupling_sd_display,'top');
        set(gca,'Xtick',xtick,'XTickLabel',cfgANA.questAnalysis_effect_str,'fontsize',8), 
        rotateticklabel(gca,50);  set(gca,'OuterPosition',[0 0.3 1 .6]), ylim([0 10]),
        title('Perceived coupling with romantic partner'), ylim([0 12]), 
        ylabel('Perceived coupling during previous run'), colormap summer,
        clear emotion_mean_display emotion_sd_display coupling_mean_display coupling_sd_display xtick
    end  
    
    
    % *** Compute and display correlation between emotion intensity and perceived coupling
    if(disp_corrSR == 1)
        % organize data
        scoreEmotion            = mean(QUEST.scoreEmotion.all_cond);
        [scoreCoupling sort_ix] = sort(mean(QUEST.scoreCoupling.all_cond)); % this sort is for drawing a straight regression line 
        scoreEmotion            = scoreEmotion(sort_ix);
        % process spearman correlation and confidence interval
        [QUEST.corrEmotionCoupling.rho,QUEST.corrEmotionCoupling.p] = corr(scoreCoupling',scoreEmotion','type','Spearman');
        QUEST.corrEmotionCoupling.upper = tanh(atanh(QUEST.corrEmotionCoupling.rho)+(1.96/sqrt(length(scoreEmotion)-3))); % 95% confidence interval for spearman rank correlation
        QUEST.corrEmotionCoupling.lower = tanh(atanh(QUEST.corrEmotionCoupling.rho)-(1.96/sqrt(length(scoreEmotion)-3)));
        % display correlation scattered plot along with regression line
        rho_str = ['r = ' num2str(QUEST.corrEmotionCoupling.rho,'%0.3f')];	p_str = ['p = ' num2str(QUEST.corrEmotionCoupling.p,'%0.3f')];
        figure('Name','Correlation emotion / coupling','WindowStyle','docked','NumberTitle','off'),
        plot(scoreCoupling,scoreEmotion,'o'), hold on,
        p = polyfit(scoreCoupling,scoreEmotion,1); f = polyval(p,scoreCoupling); plot(scoreCoupling,f,'color','r','linewidth',1),
        xlabel('perceived coupling'), ylabel('emotion intensity'), text(.8,.2,{rho_str;p_str},'Units','Normalized','fontsize',9),
        clear scoreEmotion scoreCoupling sort_ix rho_str p_str p f 
    end
    
    
  	% *** Compute correlation between empathy scores and effects obtained on self reports
    effect_emotion  = squeeze(QUEST.scoreEmotion.grouped_effect(:,2,:)) - squeeze(QUEST.scoreEmotion.grouped_effect(:,1,:));
   	effect_coupling = squeeze(QUEST.scoreCoupling.grouped_effect(:,2,:)) - squeeze(QUEST.scoreCoupling.grouped_effect(:,1,:));
    n_samples       = length(effect_emotion);
    % compute (partial) correlations and significance for each type of empathy scores and for each effect
    all_variables_emotion           = [effect_emotion ; QUEST.scoreEmpathy.grouped_effect]; % matrix with all effects on emotion intensity along with empathy scores 
    [QUEST.scoreEmotion.partial_corr_rho, QUEST.scoreEmotion.corr_rho, QUEST.scoreEmotion.corr_p] = procPartialCorr(all_variables_emotion','Spearman');
    QUEST.scoreEmotion.corr_upper	= tanh(atanh(QUEST.scoreEmotion.corr_rho)+(1.96/sqrt(n_samples-3))); % 95% confidence interval for spearman rank correlation
  	QUEST.scoreEmotion.corr_lower   = tanh(atanh(QUEST.scoreEmotion.corr_rho)-(1.96/sqrt(n_samples-3)));
    all_variables_coupling          = [effect_coupling ; QUEST.scoreEmpathy.grouped_effect]; % matrix with all effects on perceived coupling along with empathy scores 
    [QUEST.scoreCoupling.partial_corr_rho, QUEST.scoreCoupling.corr_rho, QUEST.scoreCoupling.corr_p] = procPartialCorr(all_variables_coupling','Spearman');
    QUEST.scoreCoupling.corr_upper 	= tanh(atanh(QUEST.scoreCoupling.corr_rho)+(1.96/sqrt(n_samples-3))); % 95% confidence interval for spearman rank correlation
  	QUEST.scoreCoupling.corr_lower 	= tanh(atanh(QUEST.scoreCoupling.corr_rho)-(1.96/sqrt(n_samples-3)));
    clear effect_emotion effect_coupling n_samples 
    
    
    % *** Display correlation between empathy scores and effects obtained on self reports
    if (disp_corrEmpathySR==1)

        % first display the (partial) correlation matrices (bottom: correlations, top: partial correlations)
        empathy_str     = {'empathy','contagion','emoCut','global'};
        label_all       = [cfgANA.questAnalysis_effect empathy_str];
        corr_emotion_display = triu(QUEST.scoreEmotion.corr_rho) + tril(QUEST.scoreEmotion.partial_corr_rho); 
        corr_emotion_display(QUEST.scoreEmotion.corr_p>=.05)=0;
    	corr_coupling_display = triu(QUEST.scoreCoupling.corr_rho) + tril(QUEST.scoreCoupling.partial_corr_rho);
        corr_coupling_display(QUEST.scoreCoupling.corr_p>=.05)=0;
        figure('Name','(partial) correlations emotion intensity / empathy scores','WindowStyle','docked','NumberTitle','off'),
        imagescMatrixSelection(corr_emotion_display,[length(cfgANA.questAnalysis_effect),length(empathy_str)]), hold on,
        set(gca,'YTick',1:length(corr_emotion_display)), set(gca,'XTick',[]), caxis([-1 1]),
        set(gca,'YTickLabel',label_all), title('(partial) correlations emotion intensity / empathy scores'),
        figure('Name','(partial) correlations perceived coupling / empathy scores','WindowStyle','docked','NumberTitle','off'),
        imagescMatrixSelection(corr_coupling_display,[length(cfgANA.questAnalysis_effect),length(empathy_str)]), hold on,
        set(gca,'YTick',1:length(corr_coupling_display)), set(gca,'XTick',[]),caxis([-1 1]),
        set(gca,'YTickLabel',label_all), title('(partial) correlations perceived coupling / empathy scores'),
                
        % then display scattered plot for all significant (p<0.05) correlations
        % do it for emotion intensity
        significant_corr = QUEST.scoreEmotion.corr_p + tril(ones(size(QUEST.scoreEmotion.corr_p,1),size(QUEST.scoreEmotion.corr_p,1))); % mask botton left triangular part
        [sig_row,sig_col] = find(significant_corr(1:length(cfgANA.questAnalysis_effect),:) < .05);   % keep only values that are < 0.05 and do not consider correlations within empathy scores 
        for corr_ix = 1:length(sig_col)
            [var1 sort_ix]  = sort(all_variables_emotion(sig_col(corr_ix),:)); % this sort is for drawing a straight regression line 
            var2            = all_variables_emotion(sig_row(corr_ix),sort_ix);
            rho_str = ['r = ' num2str(QUEST.scoreEmotion.corr_rho(sig_col(corr_ix),sig_row(corr_ix)),'%0.3f')];
            p_str   = ['p = ' num2str(QUEST.scoreEmotion.corr_p(sig_col(corr_ix),sig_row(corr_ix)),'%0.3f')];
            figure('Name','Correlation scattered plot','WindowStyle','docked','NumberTitle','off'),
            plot(var1,var2,'o'), hold on, line(get(gca,'Xlim'),[0 0],'color','k','linewidth',1,'linestyle','--');
            p = polyfit(var1,var2,1); f = polyval(p,var1); plot(var1,f,'color','r','linewidth',1),
          	text(.8,.2,{rho_str;p_str},'Units','Normalized','fontsize',9),
            xlabel(label_all{sig_col(corr_ix)}), ylabel(label_all{sig_row(corr_ix)}), title('emotion intensity'), 
        end
        % do it for perceived coupling
        significant_corr = QUEST.scoreCoupling.corr_p + tril(ones(size(QUEST.scoreCoupling.corr_p,1),size(QUEST.scoreCoupling.corr_p,1))); % mask botton left triangular part
        [sig_row,sig_col] = find(significant_corr(1:length(cfgANA.questAnalysis_effect),:) < .05);   % keep only values that are < 0.05 and do not consider correlations within empathy scores 
        for corr_ix = 1:length(sig_col)
            [var1 sort_ix]  = sort(all_variables_coupling(sig_col(corr_ix),:)); % this sort is for drawing a straight regression line 
            var2            = all_variables_coupling(sig_row(corr_ix),sort_ix);
            rho_str = ['r = ' num2str(QUEST.scoreCoupling.corr_rho(sig_col(corr_ix),sig_row(corr_ix)),'%0.3f')];
            p_str   = ['p = ' num2str(QUEST.scoreCoupling.corr_p(sig_col(corr_ix),sig_row(corr_ix)),'%0.3f')];
            figure('Name','Correlation scattered plot','WindowStyle','docked','NumberTitle','off'),
            plot(var1,var2,'o'), hold on, line(get(gca,'Xlim'),[0 0],'color','k','linewidth',1,'linestyle','--');
            p = polyfit(var1,var2,1); f = polyval(p,var1); plot(var1,f,'color','r','linewidth',1),
            text(.8,.2,{rho_str;p_str},'Units','Normalized','fontsize',9),
            xlabel(label_all{sig_col(corr_ix)}), ylabel(label_all{sig_row(corr_ix)}), title('perceived coupling'), 
        end
        clear empathy_str corr_emotion_display corr_coupling_display p p_str rho_str sig_col 
        clear sig_row significant_corr sort_ix var1 var2 label_all f corr_ix all_variables_coupling all_variables_emotion
    end

    
	% *** Display empathy scores 
    if (disp_empathy == 1)
        empathy_mean_display	= mean(QUEST.scoreEmpathy.grouped_effect,2);
        empathy_sd_display      = std(QUEST.scoreEmpathy.grouped_effect,0,2);
        figure('Name',['Empathy scores, couples: ' cfgANA.couple_set],'WindowStyle','docked','NumberTitle','off'),
        bar(empathy_mean_display), hold on, errorb(empathy_mean_display,empathy_sd_display);
        ylim([0 10]), set(gca,'XTickLabel',{'empathy','contagion','cut','global'}), title('Empathy scores to Favre''s questionnaire'),
        ylabel('score'), colormap summer,
        clear empathy_mean_display empathy_sd_display
    end   
    
    
    % *** Display results on harmony in the romantic couple
    if (disp_harmony == 1)
        score_harmony_display = reshape(cell2mat(QUEST.score_harmony),[length(cfgANA.coupleID) 2]);
        couple_labels = [repmat('#',[length(cfgANA.coupleID) 1]) , int2str(cfgANA.coupleID')];
        figure('Name',['Harmony scores, couples: ' cfgANA.couple_set],'WindowStyle','docked','NumberTitle','off'),
        bar(score_harmony_display,'group','BarWidth',1), ylim([1 5]), legend('sub1','sub2'), colormap summer,
        xlabel('Couple index'), set(gca,'XTickLabel',couple_labels), ylabel('Harmony score [1-5]'), 
        clear score_harmony_display couple_labels
    end
end




%% % *** B/ Analysis of physiological signals   
if cfgANA.procedureB_physio == 1

    % --- B1) B2) B3) Analysis of Pulse Rate Variability (PRV), Breath-to-Breath Variability (BBV), or Electrodermal Activity (EDA) 
    cfgDISP.physio_allExpe          = 0;   	% display physio signal in time for each subject for all experiment
    cfgDISP.physioTrend_allExpe     = 0;    % display physiological trend for each subject during all experiment
    cfgDISP.xcorr_sub_allExpe       = 0;  	% display for each subject physio auto- and cross-correlation for all experiment
    cfgDISP.xcorr_couple_cond       = 0;    % display for each couple physio cross-correlation averaged over conditions of interest
    cfgDISP.xcorr_lags              = 0;    % display distribution of lags of max xcorr values in time
    cfgDISP.xcorr_couple_averages   = 1;    % display average physio xcorr for each condition in INTRA and INTER couple
    cfgDISP.couple_correlations   	= 0;    % display correlations and partial correlation across all subjects in contrasted conditions
    cfgDISP.PRV_PSD_sub             = 0;    % display PSD of PRV on all experiment for each subject

    
    % *** various init
	if (~strcmp(cfgANA.dataType,'all') && ~strcmp(cfgANA.dataType,'EXT_only')) 
        error('Need external channels for analysis of physiological signals!'),
    end
    if (cfgANA.manualArtifRej == 1)
        error('Manual artifact rejection must be set to 0 for analysis of physiological signals.'),
    end        
    switch cfgANA.physio_sig
        % 1. Pulse Rate Variability
        case 'PRV'
            cfgANA.physio_str  = 'Pulse Rate Variability (PRV)';
            cfgANA.physio_str2 = '\Delta{P-P} (s)';         % this is a figure label
            physio_sig = 'PRV';
            physio_chan = ismember(cfgEEG.chan_str(cfgEEG.chanSet(1,:)),'Pulse');
        % 2. Respiration Volume per Time
        case 'RVT'
         	cfgANA.physio_str  = 'Respiration Volume per Time (RVT)';
          	cfgANA.physio_str2 = '\Delta{V} / \Delta{B-B}'; % this is a figure label
            physio_sig = 'RVT';
            physio_chan = ismember(cfgEEG.chan_str(cfgEEG.chanSet(1,:)),'Resp');
        % 3. Electrodermal Activity
        case 'EDA'
         	cfgANA.physio_str  = 'Electrodermal Activity (EDA)';
         	cfgANA.physio_str2 = 'normalized EDA';          % this is a figure label
            physio_sig = 'EDA';
            physio_chan = ismember(cfgEEG.chan_str(cfgEEG.chanSet(1,:)),'GSR');
        otherwise
            error('unrecognized signal type')
    end
    % concatenate data from all subjects
    EXT.([physio_sig '_allCouples']).all_raw = zeros(cfgEEG.nbSubTot,cfgEEG.lengthExpe);   % struct with raw physio data for all experiment
    for sub_ix = 1 : cfgEEG.nbSubTot
         EXT.([physio_sig '_allCouples']).all_raw(sub_ix,:) = EEG.data_sub(sub_ix).all(physio_chan,:); % format: (nbSubs, lengthExpe)
    end
    % TODO: create generic structure for all measures derived from physio data 
    % cfgANA.derived_measures = 0;        % process and analyse (time and/or freq) measures derived from physio signals 

 
    % *** process physio data for all experiment
 	disp(' '); disp(['Calculation of ' cfgANA.physio_str '...']);
    switch physio_sig
        case 'PRV'
            [EXT.([physio_sig '_allCouples']).all,  EXT.([physio_sig '_allCouples']).all_trend] = procPRV(EXT.([physio_sig '_allCouples']).all_raw,cfgEEG.Fs,'derivative2'); % process PRV
            EXT.([physio_sig '_allCouples']).PRV_measures_all = procPRV_measures(EXT.([physio_sig '_allCouples']).all, cfgEEG.Fs, cfgDISP.PRV_PSD_sub);
        case 'RVT'
         	[EXT.([physio_sig '_allCouples']).all,  EXT.([physio_sig '_allCouples']).all_trend] = procRVT(EXT.([physio_sig '_allCouples']).all_raw,cfgEEG.Fs); % process RVT
        case 'EDA'
         	[EXT.([physio_sig '_allCouples']).all,  EXT.([physio_sig '_allCouples']).all_trend] = procEDA(EXT.([physio_sig '_allCouples']).all_raw,cfgEEG.Fs,0); % preprocess EDA
         	EXT.([physio_sig '_allCouples']).AMP = procEDA_measures(EXT.([physio_sig '_allCouples']).all,TRIG.run_begin_samples_sub,cfgEEG.Fs,0);
        otherwise
            error('unrecognized signal type')
    end
    
    % *** compute derived PRV variables for each contrasted condition
    if (cfgANA.derived_measures == 1) && (strcmp(cfgANA.physio_sig,'PRV'))
        for cond_ix = 1:cfgEEG.nbCondTot
            curr_data_cond = indexMat(EXT.([physio_sig '_allCouples']).all, 1:cfgEEG.nbSubTot, TRIG.trigAll_sub(:,:,cond_ix));
            curr_PRV_measures = procPRV_measures(curr_data_cond, cfgEEG.Fs, cfgDISP.PRV_PSD_sub);
            EXT.([physio_sig '_allCouples']).PRV_measures_cond.MEAN(:,cond_ix)  = curr_PRV_measures.MEAN;
            EXT.([physio_sig '_allCouples']).PRV_measures_cond.SD(:,cond_ix)    = curr_PRV_measures.SD;
            EXT.([physio_sig '_allCouples']).PRV_measures_cond.PSD(:,:,cond_ix) = curr_PRV_measures.PSD;
            EXT.([physio_sig '_allCouples']).PRV_measures_cond.VLF(:,cond_ix)   = curr_PRV_measures.VLF;
            EXT.([physio_sig '_allCouples']).PRV_measures_cond.LF(:,cond_ix)    = curr_PRV_measures.LF;
            EXT.([physio_sig '_allCouples']).PRV_measures_cond.HF(:,cond_ix)    = curr_PRV_measures.HF;    
        end
        EXT.([physio_sig '_allCouples']).PRV_measures_cond.PSD_faxis = curr_PRV_measures.PSD_faxis;
    end
    
    % *** process coupling on physio data for all experiment, and for all conditions separately
    switch cfgANA.coupling_measure
      	% 1. cross-correlation
        case 'xcorr' 
            disp('Computing corr, partial corr, and xcorr across all couples...'),
            % ** various init
            EXT.coupling_length_all     = 2*cfgEEG.lengthExpe-1;    % xcorr length for all experiment
            EXT.coupling_length_cond    = 2*cfgEEG.lengthCond-1;    % xcorr length for each condition
            EXT.coupling_timeAxis_all   = ((0:EXT.coupling_length_all-1) - EXT.coupling_length_all/2) ./ cfgEEG.Fs;     % time axis for all experiment
            EXT.coupling_timeAxis_cond  = ((0:EXT.coupling_length_cond-1) - EXT.coupling_length_cond/2) ./ cfgEEG.Fs;  	% time axis for each condition
            cfgANA.nb_coupling_inter     = .5*cfgEEG.nbSubTot*(cfgEEG.nbSubTot-2); % N(N-2)/2 , with N = nsubs
            
            % ** compute corr and partial corr coefs
            [EXT.([physio_sig '_allCouples']).partial_corr_all, EXT.([physio_sig '_allCouples']).corr_all] = ...
                procPartialCorr(EXT.([physio_sig '_allCouples']).all','Pearson'); 
            
            % ** compute autocorr for whole experiment
            EXT.([physio_sig '_allCouples']).autocorr_all = zeros(cfgEEG.nbSubTot,EXT.coupling_length_all);
            for sub_ix = 1:cfgEEG.nbSubTot
                EXT.([physio_sig '_allCouples']).autocorr_all(sub_ix,:) = xcorr(EXT.([physio_sig '_allCouples']).all(sub_ix,:),'coeff'); % autocorrelation. format: (2xlengthExpe-1,nbSub)
            end
            
            % ** compute for whole experiment xcorr within each couple (far too heavy to compute across all couples!)
            EXT.([physio_sig '_allCouples']).xcorr_all_intra = zeros(cfgEEG.nbCouple,EXT.coupling_length_all);
            for couple_ix = 1:cfgEEG.nbCouple
                EXT.([physio_sig '_allCouples']).xcorr_all_intra(couple_ix,:) = xcorr(EXT.([physio_sig '_allCouples']).all(2*couple_ix-1,:),EXT.([physio_sig '_allCouples']).all(2*couple_ix,:),'coeff'); % xcorr between all subjects. format: (2xlengthExpe-1,nbCouples)
            end
            
            % ** compute xcorr across all couples for each condition separately
            cfgANA.index_autocorr = sub2ind([cfgEEG.nbSubTot cfgEEG.nbSubTot],1:cfgEEG.nbSubTot,1:cfgEEG.nbSubTot); % mask vector to keep only N autocorr 
            cfgANA.index_xcorr_intra = sub2ind([cfgEEG.nbSubTot cfgEEG.nbSubTot], 1:2:cfgEEG.nbSubTot, 2:2:cfgEEG.nbSubTot); % mask vector to keep only N xcorr corresponding to subjects from same couples 
            cfgANA.index_xcorr_inter = tril(ones(cfgEEG.nbSubTot,cfgEEG.nbSubTot),-1);
            for sub_ix = 1:2:cfgEEG.nbSubTot
                cfgANA.index_xcorr_inter(sub_ix+1,sub_ix) = 0;
            end
            cfgANA.index_xcorr_inter = logical(cfgANA.index_xcorr_inter(:)); % mask vector to keep only N(N-2)/2 xcorr corresponding to subjects from different couples 
            xcorr_cond = zeros(cfgEEG.nbSubTot^2,EXT.coupling_length_cond,cfgEEG.nbCondTot);
            for cond_ix = 1:cfgEEG.nbCondTot
                curr_data = indexMat(EXT.([physio_sig '_allCouples']).all, 1:cfgEEG.nbSubTot, TRIG.trigAll_sub(:,:,cond_ix));
                xcorr_cond(:,:,cond_ix) = xcorr(curr_data','coeff')';
            end
            EXT.([physio_sig '_allCouples']).xcorr_cond_intra =  xcorr_cond(cfgANA.index_xcorr_intra,:,:); % format: (nbCouples , 2xlengthExpe-1, nbCond)
            EXT.([physio_sig '_allCouples']).xcorr_cond_inter =  xcorr_cond(cfgANA.index_xcorr_inter,:,:); % format: (nbCouples(nbCouples-1) , 2xlengthExpe-1, nbCond)
            disp('--> xcorr computed INTRA and INTER couples.'),
            clear xcorr_all_intra autocorr_all couple_ix sub_ix row_ix col_ix curr_cond_ix cond_ix xcorr_cond curr_data

      	% 2. phase locking value
        case 'PLV'
            % ** various init
            EXT.coupling_length_all     = cfgEEG.lengthExpe;    % xcorr length for all experiment
            EXT.coupling_length_cond    = cfgEEG.lengthCond;    % xcorr length for each condition
            EXT.coupling_timeAxis_all   = (0:EXT.coupling_length_all-1) ./ cfgEEG.Fs;     % time axis for all experiment
            EXT.coupling_timeAxis_cond  = (0:EXT.coupling_length_cond-1)./ cfgEEG.Fs;      % time axis for each condition
            
            % ** INTRA: compute PLV between subject from same couple
            disp('Computing PLV INTRA couples...'),
            dataset1 = EXT.([physio_sig '_allCouples']).all((1:2:cfgEEG.nbSubTot),:); % contains physio data for first half of all subjects
            dataset2 = EXT.([physio_sig '_allCouples']).all((2:2:cfgEEG.nbSubTot),:); % contains physio data for second half of all subjects
            EXT.([physio_sig '_allCouples']).PLV_all_intra = procPLV_mat(dataset1, dataset2, cfgANA.PLV_winSize, .75, 0); 

            % ** INTER: compute PLV between subject from different couple (can be only a random part of all possible combinations)
            part_inter_to_compute = 10; % part of coupling inter couple to compute (in percents, from 1 to 100)
            inter_ix = nchoosek(1:cfgEEG.nbSubTot,2); intra_ix = find([1;diff(inter_ix(:,1))]); 
            intra_ix = intra_ix(1:2:end); inter_ix(intra_ix,:)=[]; % inter_ix: indexes of subjects from different couples, format: (N(N-2)/2, 2)
            cfgANA.nb_coupling_inter = fix(length(inter_ix)*(part_inter_to_compute/100)); 
            rand_part = randperm(length(inter_ix),cfgANA.nb_coupling_inter);
            % permuting condition periods before calculation of PLV so we have identical conditions across subjects from different couples. MARCHE POOOOOO!!!!!!!
%             trigAllConcatenate = reshape(TRIG.trigAll_couples,cfgEEG.nbCouple,EXT.coupling_length_cond*cfgEEG.nbCondTot);
%          	  row_ix = repmat(1:cfgANA.nb_coupling_inter,cfgEEG.nbCondTot*cfgEEG.lengthCond,1)';
%             col_ix = trigAllConcatenate(round(inter_ix(rand_part,1)/2),:); 
%             cond_permute_ix = sub2ind([cfgANA.nb_coupling_inter cfgEEG.lengthExpe],row_ix(:),col_ix(:)); 
%             dataset1 = EXT.([physio_sig '_allCouples']).all(round(inter_ix(rand_part,1)/2),:);
%             dataset1 = reshape(dataset1(cond_permute_ix),cfgEEG.nbCondTot*cfgEEG.lengthCond,cfgANA.nb_coupling_inter)'; 
%             col_ix = trigAllConcatenate(round(inter_ix(rand_part,2)/2),:); 
%             cond_permute_ix = sub2ind([cfgANA.nb_coupling_inter cfgEEG.lengthExpe],row_ix(:),col_ix(:)); 
%             dataset2 = EXT.([physio_sig '_allCouples']).all(round(inter_ix(rand_part,2)/2),:);
%             dataset2 = reshape(dataset2(cond_permute_ix),cfgEEG.nbCondTot*cfgEEG.lengthCond,cfgANA.nb_coupling_inter)'; 
            dataset1 = EXT.([physio_sig '_allCouples']).all(inter_ix(rand_part,1),:);  % same without permutation
            dataset2 = EXT.([physio_sig '_allCouples']).all(inter_ix(rand_part,2),:);  % same without permutation
            disp(['Computing PLV INTER couples... (' int2str(part_inter_to_compute) '% of ' int2str(length(inter_ix)) ' total combinations)']),
            EXT.([physio_sig '_allCouples']).PLV_all_inter = procPLV_mat(dataset1, dataset2, cfgANA.PLV_winSize, .75, 0); 
%             EXT.([physio_sig '_allCouples']).PLV_cond_inter = reshape(EXT.([physio_sig '_allCouples']).PLV_all_inter',cfgANA.nb_coupling_inter,EXT.coupling_length_cond,cfgEEG.nbCondTot);
            
            % ** TEMPO INTER: extract PLV for each condition separately % same without permutation
            EXT.([physio_sig '_allCouples']).PLV_cond_inter = zeros(cfgANA.nb_coupling_inter,EXT.coupling_length_cond,cfgEEG.nbCondTot);
            for cond_ix = 1:cfgEEG.nbCondTot
                col_ix = repmat(TRIG.trigAll_couples(:,:,cond_ix),27,1); col_ix=col_ix(1:cfgANA.nb_coupling_inter,:); % très sale!
                EXT.([physio_sig '_allCouples']).PLV_cond_inter(:,:,cond_ix) = indexMat(EXT.([physio_sig '_allCouples']).PLV_all_inter, 1:cfgANA.nb_coupling_inter, col_ix);
            end

            % ** INTRA: extract PLV for each condition separately
            EXT.([physio_sig '_allCouples']).PLV_cond_intra = zeros(cfgEEG.nbCouple,EXT.coupling_length_cond,cfgEEG.nbCondTot);
            for cond_ix = 1:cfgEEG.nbCondTot
                EXT.([physio_sig '_allCouples']).PLV_cond_intra(:,:,cond_ix) = indexMat(EXT.([physio_sig '_allCouples']).PLV_all_intra, 1:cfgEEG.nbCouple, TRIG.trigAll_couples(:,:,cond_ix));
            end
            % example:  row=repmat((1:3),4,1), col=([1 2 3 4; 2 3 4 5; 3 4 5 1])', a=[1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15], b=a(sub2ind([3 5],row(:),col(:))), c=reshape(b,4,3)'
            clear dataset1 dataset2 cond_ix col_ix curr_cond_ix curr_data part_inter_to_compute inter_ix intra_ix nb_inter_part rand_part
            disp('--> PLV computed INTRA couples.'),
            
        % 3. other coupling measure (mutual information, etc.)
        % see section on PLV to add a new measure
        % 
            
        otherwise
            error('unrecognized coupling measure')
    end
    
    
   % *** perform analysis across contrasted conditions
   for effect_ix = 1:length(cfgANA.EFFECT_ANALYSIS) 
        
        % ** various init
        effect = cfgANA.EFFECT_ANALYSIS{effect_ix};
        switch effect
        % 1. contrast TOUCH / no TOUCH
        case 'TOUCH'
            cfgANA.effect_contrast1     = {'noTouch_Congruent';'noTouch_noCongruent';'noTouch_Control'};
            cfgANA.effect_contrast2     = {'Touch_Congruent';'Touch_noCongruent';'Touch_Control'};   
            cfgANA.effect_contrast2_str = 'touch'; cfgANA.effect_contrast1_str = 'no touch'; 
            
        % 2. contrast EMOTION / CONTROL
        case 'EMOTION'
            cfgANA.effect_contrast1    = {'Touch_Control';'noTouch_Control'}; 
            cfgANA.effect_contrast2    = {'Touch_Congruent';'noTouch_Congruent';'Touch_noCongruent';'noTouch_noCongruent'};  
            cfgANA.effect_contrast2_str = 'emotion'; cfgANA.effect_contrast1_str = 'neutral';
        % 3. contrast CONGRUENT / no CONGRUENT
        case 'CONGRUENCE'
            cfgANA.effect_contrast1   = {'Touch_noCongruent';'noTouch_noCongruent'}; 
            cfgANA.effect_contrast2   = {'Touch_Congruent';'noTouch_Congruent'};      
            cfgANA.effect_contrast2_str = 'congruent'; cfgANA.effect_contrast1_str = 'no congruent';
        % 4. contrast EMOTION POS / EMOTION NEG
        case 'VALENCE'
            cfgANA.effect_contrast1   = {'Touch_ValenceNeg' ; 'noTouch_ValenceNeg'};  
            cfgANA.effect_contrast2   = {'Touch_ValencePos' ; 'noTouch_ValencePos'};      
            cfgANA.effect_contrast2_str = 'emotionPos'; cfgANA.effect_contrast1_str = 'emotionNeg';
        otherwise
            error('unrecognized effect type')
        end
        cfgANA.cond_contrast1 = ismember(cfgEEG.CondTrigVal,cell2mat(reshape(cfgEEG.effect_cond_ix(ismember(cfgEEG.effect_labels,cfgANA.effect_contrast1)),1,[]))); % logical, e.g. [1 1 1 1 1 0 0 0 0 0]
        cfgANA.cond_contrast2 = ismember(cfgEEG.CondTrigVal,cell2mat(reshape(cfgEEG.effect_cond_ix(ismember(cfgEEG.effect_labels,cfgANA.effect_contrast2)),1,[])));
        cfgANA.nb_cond_contrast1 = length(find(cfgANA.cond_contrast1)); % scalar
        cfgANA.nb_cond_contrast2 = length(find(cfgANA.cond_contrast2)); % scalar
        cfgANA.all_effect_str{2*effect_ix-1} = [cfgANA.effect_contrast1_str '   ']; cfgANA.all_effect_str{2*effect_ix} = [cfgANA.effect_contrast2_str '   '];
        disp(['*** Testing ' effect '... ***']);

        
    	% *** SUBJECT LEVEL: extract derived EDA variables for each contrasted conditions
        if (cfgANA.derived_measures == 1) && (strcmp(cfgANA.physio_sig,'EDA'))
            EXT.([physio_sig '_allCouples']).AMP_cond(effect_ix,1,:) = mean(EXT.([physio_sig '_allCouples']).AMP(:,cfgANA.cond_contrast1),2); % EDA AMP mean for all subjects in both contrasts (nbEffects,2,nbSubTot)
         	EXT.([physio_sig '_allCouples']).AMP_cond(effect_ix,2,:) = mean(EXT.([physio_sig '_allCouples']).AMP(:,cfgANA.cond_contrast2),2); % EDA AMP mean for all subjects in both contrasts (nbEffects,2,nbSubTot)
        end
        
        % *** SUBJECT LEVEL: compute (geometrical) mean of PRV measures within each group of contrast conditions 
        if (cfgANA.derived_measures == 1) && (strcmp(cfgANA.physio_sig,'PRV'))
            EXT.([physio_sig '_allCouples']).PRV_measures_effect.MEAN(effect_ix,1,:) = mean(EXT.([physio_sig '_allCouples']).PRV_measures_cond.MEAN(:,cfgANA.cond_contrast1),2); % PRV MEAN average for all subjects in both contrasts (nbEffects,2,nbSubTot)
            EXT.([physio_sig '_allCouples']).PRV_measures_effect.MEAN(effect_ix,2,:) = mean(EXT.([physio_sig '_allCouples']).PRV_measures_cond.MEAN(:,cfgANA.cond_contrast2),2);
            EXT.([physio_sig '_allCouples']).PRV_measures_effect.SD(effect_ix,1,:) = geomean(EXT.([physio_sig '_allCouples']).PRV_measures_cond.SD(:,cfgANA.cond_contrast1),2); % PRV SD geometrical mean for all subjects in both contrasts (nbEffects,2,nbSubTot)
            EXT.([physio_sig '_allCouples']).PRV_measures_effect.SD(effect_ix,2,:) = geomean(EXT.([physio_sig '_allCouples']).PRV_measures_cond.SD(:,cfgANA.cond_contrast2),2);
            EXT.([physio_sig '_allCouples']).PRV_measures_effect.LF(effect_ix,1,:) = geomean(EXT.([physio_sig '_allCouples']).PRV_measures_cond.LF(:,cfgANA.cond_contrast1),2); % PRV LF geometrical mean for all subjects in both contrasts (nbEffects,2,nbSubTot)
            EXT.([physio_sig '_allCouples']).PRV_measures_effect.LF(effect_ix,2,:) = geomean(EXT.([physio_sig '_allCouples']).PRV_measures_cond.LF(:,cfgANA.cond_contrast2),2);
            EXT.([physio_sig '_allCouples']).PRV_measures_effect.HF(effect_ix,1,:) = geomean(EXT.([physio_sig '_allCouples']).PRV_measures_cond.HF(:,cfgANA.cond_contrast1),2); % PRV HF geometrical mean for all subjects in both contrasts (nbEffects,2,nbSubTot)
            EXT.([physio_sig '_allCouples']).PRV_measures_effect.HF(effect_ix,2,:) = geomean(EXT.([physio_sig '_allCouples']).PRV_measures_cond.HF(:,cfgANA.cond_contrast2),2);
        end
        
        % *** COUPLE LEVEL: perform analysis across contrasted conditions depending on chosen coupling measure
        switch cfgANA.coupling_measure
            % 1. cross-correlation
            case 'xcorr' 
                coupling_intra_contrast1 = EXT.([physio_sig '_allCouples']).xcorr_cond_intra(:,:,cfgANA.cond_contrast1);    % format: (nbCouples, lengthCond, nbCondContrast1)
                coupling_intra_contrast2 = EXT.([physio_sig '_allCouples']).xcorr_cond_intra(:,:,cfgANA.cond_contrast2);    % format: (nbCouples, lengthCond, nbCondContrast2)
                coupling_inter_contrast1 = EXT.([physio_sig '_allCouples']).xcorr_cond_inter(:,:,cfgANA.cond_contrast1);    % format: (nbCouples(nbCouples-1), lengthCond, nbCondContrast1)
                coupling_inter_contrast2 = EXT.([physio_sig '_allCouples']).xcorr_cond_inter(:,:,cfgANA.cond_contrast2);     % format: (nbCouples(nbCouples-1), lengthCond, nbCondContrast2)
                
                if strcmp(cfgANA.xcorr_lag_type,'mean')	 % only testing MEAN xcorr within specified interval around lag0
                    cfgANA.xcorr_lag_index = find((EXT.coupling_timeAxis_cond>=cfgANA.xcorr_lag_interval(1))&(EXT.coupling_timeAxis_cond<=cfgANA.xcorr_lag_interval(2)));
                    coupling_intra_contrast1_mean = squeeze(zMean(coupling_intra_contrast1(:,cfgANA.xcorr_lag_index,:),2));	% format: (nbCouples, nbCondContrast1)
                    coupling_intra_contrast2_mean = squeeze(zMean(coupling_intra_contrast2(:,cfgANA.xcorr_lag_index,:),2));	% format: (nbCouples, nbCondContrast2)
                    coupling_intra_allCond_mean = squeeze(zMean(EXT.([physio_sig '_allCouples']).xcorr_cond_intra(:,cfgANA.xcorr_lag_index,:),2)); % format: (nbCouples, nbCondTot)
                    coupling_inter_contrast1_mean = squeeze(zMean(coupling_inter_contrast1(:,cfgANA.xcorr_lag_index,:),2));	% format: (nbCouples(nbCouples-1), nbCondContrast1)
                    coupling_inter_contrast2_mean = squeeze(zMean(coupling_inter_contrast2(:,cfgANA.xcorr_lag_index,:),2));	% format: (nbCouples(nbCouples-1), nbCondContrast2)
                    coupling_inter_allCond_mean = squeeze(zMean(EXT.([physio_sig '_allCouples']).xcorr_cond_inter(:,cfgANA.xcorr_lag_index,:),2)); % format: (nbCouples(nbCouples-1), nbCondTot) 
                elseif strcmp(cfgANA.xcorr_lag_type,'max')  % only testing MAX xcorr    
                    [max_intra_contrast1,max_intra_contrast1_ix] = max(coupling_intra_contrast1,[],2); % getting max BEFORE averaging over conditions. 
                    [max_intra_contrast2,max_intra_contrast2_ix] = max(coupling_intra_contrast2,[],2); 
                    [max_intra,max_intra_ix] = max(EXT.([physio_sig '_allCouples']).xcorr_cond_intra,[],2); 
                    [max_inter_contrast1,max_inter_contrast1_ix] = max(coupling_inter_contrast1,[],2); % getting max BEFORE averaging over conditions. 
                    [max_inter_contrast2,max_inter_contrast2_ix] = max(coupling_inter_contrast2,[],2); 
                    [max_inter,max_inter_ix] = max(EXT.([physio_sig '_allCouples']).xcorr_cond_inter,[],2); 
                    coupling_intra_contrast1_mean = squeeze(max_intra_contrast1); coupling_intra_contrast2_mean = squeeze(max_intra_contrast2); coupling_intra_allCond_mean = squeeze(max_intra); 
                    coupling_inter_contrast1_mean = squeeze(max_inter_contrast1); coupling_inter_contrast2_mean = squeeze(max_inter_contrast2); coupling_inter_allCond_mean = squeeze(max_inter); 
                    max_intra_contrast1_ix = squeeze(max_intra_contrast1_ix); max_intra_contrast2_ix = squeeze(max_intra_contrast2_ix); max_intra_ix = squeeze(max_intra_ix);
                    max_inter_contrast1_ix = squeeze(max_inter_contrast1_ix); max_inter_contrast2_ix = squeeze(max_inter_contrast2_ix); max_inter_ix = squeeze(max_inter_ix);
                end
                coupling_str = ['xcorr_' cfgANA.xcorr_lag_type];
                meanFct = @zMean;     % necessary to z-transform correlation coefficients before suming/averaging
          	% 2. phase locking value    
            case 'PLV'
                coupling_intra_contrast1 = EXT.([physio_sig '_allCouples']).PLV_cond_intra(:,:,cfgANA.cond_contrast1);    % format: (nbCouples, lengthCond, nbCondContrast1)
                coupling_intra_contrast2 = EXT.([physio_sig '_allCouples']).PLV_cond_intra(:,:,cfgANA.cond_contrast2);    % format: (nbCouples, lengthCond, nbCondContrast2)
                coupling_inter_contrast1 = EXT.([physio_sig '_allCouples']).PLV_cond_inter(:,:,cfgANA.cond_contrast1);    % format: (nbCouples(nbCouples-1), lengthCond, nbCondContrast1)
                coupling_inter_contrast2 = EXT.([physio_sig '_allCouples']).PLV_cond_inter(:,:,cfgANA.cond_contrast2);    % format: (nbCouples(nbCouples-1), lengthCond, nbCondContrast2)
                coupling_intra_contrast1_mean   = squeeze(mean(coupling_intra_contrast1,2));	% format: (nbCouples, nbCondContrast1)
                coupling_intra_contrast2_mean   = squeeze(mean(coupling_intra_contrast2,2));	% format: (nbCouples, nbCondContrast2)
                coupling_intra_allCond_mean     = squeeze(mean(EXT.([physio_sig '_allCouples']).PLV_cond_intra,2)); % format: (nbCouples, nbCondTot)
                coupling_inter_contrast1_mean   = squeeze(mean(coupling_inter_contrast1,2));	% format: (nbCouples(nbCouples-1), nbCondContrast1)
                coupling_inter_contrast2_mean   = squeeze(mean(coupling_inter_contrast2,2));	% format: (nbCouples(nbCouples-1), nbCondContrast2)
                coupling_inter_allCond_mean     = squeeze(mean(EXT.([physio_sig '_allCouples']).PLV_cond_inter,2)); % format: (nbCouples(nbCouples-1), nbCondTot) 
                coupling_str = 'PLV';
                meanFct = @mean;
            otherwise
                error('unrecognized coupling measure')
        end
            
        % *** INTRA: is there a difference in intersubject coupling between conditions? -> compute significance using paired t-test
        result_str = [physio_sig '_' cfgANA.coupling_measure '_effectSize_Intra'];
        RESULTS.(result_str).(effect).diff_obs = meanFct(meanFct(coupling_intra_contrast1_mean,2)) - meanFct(meanFct(coupling_intra_contrast2_mean,2));
        coupling_couple_contrast1 = meanFct(coupling_intra_contrast1_mean,2)';
        coupling_couple_contrast2 = meanFct(coupling_intra_contrast2_mean,2)';
        if (strcmp(cfgANA.couple_set,'all')) % only feasible if all couples considered so we have enough paired permutations (alpha min = 1/(2^nCouples))
            [RESULTS.(result_str).(effect).h, RESULTS.(result_str).(effect).p] = pairedPermTest(coupling_couple_contrast1, coupling_couple_contrast2, cfgANA.tTest_alpha,cfgANA.nbPerms);      
            disp(['INTRA couples, effect of ' effect ' on ' coupling_str  ': p=' num2str(RESULTS.(result_str).(effect).p)]);
        end 

        % *** INTER: is there a difference in average intersubject coupling (all conditions) computed in intra vs inter? 
        % -> compute significance with bootstrap analysis (to heavy to compute paired permutation tests on stats INTER)
     	grandMean_intra = meanFct(meanFct(coupling_intra_allCond_mean,2));
        result_str = [physio_sig '_' cfgANA.coupling_measure '_average_IntraInter'];
        cfgANA.bootMat_mean = bootSampling(cfgANA.nb_coupling_inter, cfgEEG.nbCouple, cfgANA.nBoot);
        RESULTS.(result_str).bootRes = zeros(cfgANA.nBoot,1); % contains results of tests for all bootstrap samples
        for boot_ix = 1:cfgANA.nBoot
            bootSample = coupling_inter_allCond_mean(cfgANA.bootMat_mean(boot_ix,:));
            RESULTS.(result_str).bootRes(boot_ix) = meanFct(bootSample);
        end
        RESULTS.(result_str).p = (1/cfgANA.nBoot)*(1+length(find(RESULTS.(result_str).bootRes > grandMean_intra)));   
        RESULTS.(result_str).h = RESULTS.(result_str).p < cfgANA.tTest_alpha;  
        disp(['INTRA vs INTER couples, difference in average ' coupling_str ': p=' num2str(RESULTS.(result_str).p)]);

        % *** INTER: is there a difference in intersubject coupling between conditions? 
        % -> compute significance with bootstrap analysis (to heavy to compute paired permutation tests on stats INTER)
        result_str = [physio_sig '_' cfgANA.coupling_measure '_effectSize_Inter'];
        RESULTS.(result_str).(effect).mean_effect = meanFct(meanFct(coupling_inter_contrast1_mean,2)) - meanFct(meanFct(coupling_inter_contrast2_mean,2));
        cfgANA.bootMat_contrast1 = bootSampling(cfgANA.nb_coupling_inter, cfgEEG.nbCouple, cfgANA.nBoot);
        RESULTS.(result_str).(effect).bootRes = zeros(cfgANA.nBoot,1); % contains results of xcorr diff at specified lags for all bootstrap samples    
        for boot_ix = 1:cfgANA.nBoot
            bootSample1 = coupling_inter_contrast1_mean(cfgANA.bootMat_mean(boot_ix,:));
            bootSample2 = coupling_inter_contrast2_mean(cfgANA.bootMat_mean(boot_ix,:));
         	RESULTS.(result_str).(effect).bootRes(boot_ix) = real(meanFct(bootSample1-bootSample2));
        end                
        RESULTS.(result_str).(effect).p = (1/cfgANA.nBoot)*(1+length(find(RESULTS.(result_str).(effect).bootRes > RESULTS.(result_str).(effect).mean_effect)));   
        RESULTS.(result_str).(effect).h = RESULTS.(result_str).(effect).p < cfgANA.tTest_alpha;         
        disp(['INTER couples, effect of ' effect ' on ' coupling_str ': p=' num2str(RESULTS.(result_str).(effect).p)]);

        % *** INTRA: is there a difference in intersubject coupling between conditions after mean coupling INTER was removed? 
        % -> compute significance using paired t-test
        result_str = [physio_sig '_' cfgANA.coupling_measure '_effectSize_IntraNoInter'];
        RESULTS.(result_str).(effect).diff_obs = meanFct(meanFct(coupling_intra_contrast1_mean,2)) - meanFct(meanFct(coupling_intra_contrast2_mean,2));
        if (strcmp(cfgANA.couple_set,'all')) % only feasible if all couples considered so we have enough paired permutations (alpha min = 1/(2^nCouples))
            coupling_mean_intraNoInter_contrast1 = meanFct(coupling_intra_contrast1_mean,2)' - meanFct(meanFct(coupling_inter_contrast1_mean,2));
            coupling_mean_intraNoInter_contrast2 = meanFct(coupling_intra_contrast2_mean,2)' - meanFct(meanFct(coupling_inter_contrast2_mean,2)); 
            [RESULTS.(result_str).(effect).h, RESULTS.(result_str).(effect).p] = ...
                pairedPermTest(coupling_mean_intraNoInter_contrast1, coupling_mean_intraNoInter_contrast2, cfgANA.tTest_alpha,cfgANA.nbPerms);      
            disp(['INTRA couples with INTER removed, effect of ' effect ' on ' coupling_str ': p=' num2str(RESULTS.(result_str).(effect).p)]);
        end 

        % *** VARIOUS DISPLAY RELATED TO TESTED EFFECT
        % ** INTRA (at couple level): display physio cross-correlation averaged over conditions of interest
        if (cfgDISP.xcorr_couple_cond == 1) && (strcmp(cfgANA.coupling_measure,'xcorr'))
            for couple_ix = 1:cfgEEG.nbCouple
            	coupling_couple_mean            = zMean(squeeze(EXT.([physio_sig '_allCouples']).xcorr_cond_intra(couple_ix,:,:)),2);
                coupling_couple_contrast1_mean  = zMean(squeeze(coupling_intra_contrast1(couple_ix,:,:)),2);
                coupling_couple_contrast2_mean  = zMean(squeeze(coupling_intra_contrast2(couple_ix,:,:)),2);
                coupling_couple_contrast1_prctl = [zPrctile(squeeze(coupling_intra_contrast1(couple_ix,:,:))',5); zPrctile(squeeze(coupling_intra_contrast1(couple_ix,:,:))',95)]; 
                coupling_couple_contrast2_prctl = [zPrctile(squeeze(coupling_intra_contrast2(couple_ix,:,:))',5); zPrctile(squeeze(coupling_intra_contrast2(couple_ix,:,:))',95)]; 
                figure('Name',['Couple n°' int2str(couple_ix) ' ' physio_sig ' xcorr within cond'],'WindowStyle','docked','NumberTitle','off'),
                plot(EXT.coupling_timeAxis_cond, coupling_couple_mean, 'color','b','linewidth',2),
                hold on,
                plot(EXT.coupling_timeAxis_cond, coupling_couple_contrast2_mean,'color','r','linewidth',2),
                plot(EXT.coupling_timeAxis_cond, coupling_couple_contrast1_mean,'color','g','linewidth',2),
                legend('mean',cfgANA.effect_contrast2_str,cfgANA.effect_contrast1_str),
                jbfill(EXT.coupling_timeAxis_cond,coupling_couple_contrast2_prctl(2,:),coupling_couple_contrast2_prctl(1,:),[255 150 150]./255,'r');
                jbfill(EXT.coupling_timeAxis_cond,coupling_couple_contrast1_prctl(2,:),coupling_couple_contrast1_prctl(1,:),[150 255 150]./255,'g');
                grid on, ylim([-1 1]), line([0 0],get(gca,'Ylim'),'color','k','linewidth',1)
                title([physio_sig ' cross-correlation averaged on condition periods']), xlabel('lag (s)'),ylabel('correlation coefficient'),
                clear coupling_couple_mean coupling_couple_contrast1_mean coupling_couple_contrast2_mean coupling_couple_contrast1_prctl coupling_couple_contrast2_prctl
            end
        end
        
        % ** INTRA & INTER: visualize distribution of lags of max xcorr values in time
        if strcmp(cfgANA.coupling_measure,'xcorr') && strcmp(cfgANA.xcorr_lag_type,'max') && (cfgDISP.xcorr_lags == 1)
            fontSize = 11;
            lags_max_intra_contrast1 = EXT.coupling_timeAxis_cond(max_intra_contrast1_ix); % this gets time values from xcorr indexes
            lags_max_intra_contrast2 = EXT.coupling_timeAxis_cond(max_intra_contrast2_ix);
            lags_max_inter_contrast1 = EXT.coupling_timeAxis_cond(max_inter_contrast1_ix); % this gets time values from xcorr indexes
            lags_max_inter_contrast2 = EXT.coupling_timeAxis_cond(max_inter_contrast2_ix);
            bins = 1.2*min(EXT.coupling_timeAxis_cond) : 2 : 1.2*max(EXT.coupling_timeAxis_cond);
            figure('Name',['INTRA couples, maxima of ' physio_sig ' xcorr.'],'WindowStyle','docked','NumberTitle','off'),
            n_elem_contrast1 =  hist(lags_max_intra_contrast1(:),bins); n_elem_contrast2 =  hist(lags_max_intra_contrast2(:),bins);
            max_y = 1.2*max([n_elem_contrast1 n_elem_contrast2]);
            subplot(2,1,1), hist(lags_max_intra_contrast1(:),bins), ylim([0 max_y]), set(findobj(gca,'Type','patch'),'EdgeColor','g','FaceColor','g');
            set(gca,'layer','top','FontSize',fontSize), legend(['INTRA ' cfgANA.effect_contrast1_str]),%title(['histogram of maximum correlation lags, ' cfgANA.effect_contrast1_str]), xlabel('time(s)'),
            subplot(2,1,2), hist(lags_max_intra_contrast2(:),bins), ylim([0 max_y]), set(findobj(gca,'Type','patch'),'EdgeColor','r','FaceColor','r'), 
            set(gca,'layer','top','FontSize',fontSize), legend(['INTRA ' cfgANA.effect_contrast2_str]),
            ylim([0 max_y]), xlabel('lag (s)','FontSize',fontSize), ylabel('xcorr max occurrences','FontSize',fontSize), xlim([-100 100]), % title(['histogram of maximum correlation lags, ' cfgANA.effect_contrast2_str]), 
            figure('Name',['INTER couples, maxima of ' physio_sig ' xcorr.'],'WindowStyle','docked','NumberTitle','off'),
            n_elem_contrast1 =  hist(lags_max_inter_contrast1(:),bins); n_elem_contrast2 =  hist(lags_max_inter_contrast2(:),bins);
            max_y = 1.2*max([n_elem_contrast1 n_elem_contrast2]);
            subplot(2,1,1), hist(lags_max_inter_contrast1(:),bins), ylim([0 max_y]), set(findobj(gca,'Type','patch'),'EdgeColor','g','FaceColor','g');
            set(gca,'layer','top','FontSize',fontSize), legend(['INTER ' cfgANA.effect_contrast1_str]),%title(['histogram of maximum correlation lags, ' cfgANA.effect_contrast1_str]), xlabel('time(s)'),
            subplot(2,1,2), hist(lags_max_inter_contrast2(:),bins), ylim([0 max_y]), set(findobj(gca,'Type','patch'),'EdgeColor','r','FaceColor','r'), 
            set(gca,'layer','top','FontSize',fontSize), legend(['INTER ' cfgANA.effect_contrast2_str]),
            ylim([0 max_y]), xlabel('lag (s)','FontSize',fontSize), ylabel('xcorr max occurrences','FontSize',fontSize), xlim([-100 100]), % title(['histogram of maximum correlation lags, ' cfgANA.effect_contrast2_str]), 
            clear fontSize lags_max_intra_contrast1 lags_max_intra_contrast2 lags_max_inter_contrast1 lags_max_inter_contrast2 bins n_elem_contrast1 n_elem_contrast2 max_y
        end

        % ** INTRA, effect size on coupling: display bar plot with effect size and deviation
        fontSize = 9;
        result_str = [physio_sig '_' cfgANA.coupling_measure '_effectSize_Intra'];
        coupling_display_mean   = [zMean(coupling_intra_contrast2_mean(:)), zMean(coupling_intra_contrast1_mean(:))];
        coupling_display_std    = [zStd(coupling_intra_contrast2_mean(:)), zStd(coupling_intra_contrast1_mean(:))];
        figure('Name',['INTRA, effect size on ' coupling_str],'WindowStyle','docked','NumberTitle','off'),
        barPlot(coupling_display_mean,RESULTS.(result_str).(effect).p,[.05 .01 .001],1), hold on, 
        errorb(coupling_display_mean,coupling_display_std,'top');
        set(gca,'XTickLabel',{[cfgANA.effect_contrast2_str '  '], [cfgANA.effect_contrast1_str '  ']},'fontsize',10), 
        rotateticklabel(gca,50);  set(gca,'OuterPosition',[0 0.1 1 .9],'FontSize',fontSize);
        ylim([0 1.2]), ylabel(['average ' coupling_str],'FontSize',fontSize), colormap summer,
        clear result_str coupling_display_mean coupling_display_std

        % ** INTRA & INTER, effect size on coupling: display physio cross-correlation averaged over conditions of interest
        if strcmp(cfgANA.coupling_measure,'xcorr') && (cfgDISP.xcorr_couple_averages==1)
            coupling_intra_contrast1_allTrials = reshape(permute(coupling_intra_contrast1,[2 1 3]),size(coupling_intra_contrast1,2),[])';  % concatenate 3rd dimension (trials) on the rows -> new format:  (nbCouples x nbCondContrast1, lengthCouplingCond)
            coupling_intra_contrast2_allTrials = reshape(permute(coupling_intra_contrast2,[2 1 3]),size(coupling_intra_contrast2,2),[])';  % concatenate 3rd dimension (trials) on the rows
            coupling_inter_contrast1_allTrials = reshape(permute(coupling_inter_contrast1,[2 1 3]),size(coupling_inter_contrast1,2),[])';  % concatenate 3rd dimension (trials) on the rows
            coupling_inter_contrast2_allTrials = reshape(permute(coupling_inter_contrast2,[2 1 3]),size(coupling_inter_contrast2,2),[])';  % concatenate 3rd dimension (trials) on the rows
            coupling_intra_contrast1_display_mean = zMean(coupling_intra_contrast1_allTrials); coupling_intra_contrast2_display_mean = zMean(coupling_intra_contrast2_allTrials);
            coupling_inter_contrast1_display_mean = zMean(coupling_inter_contrast1_allTrials); coupling_inter_contrast2_display_mean = zMean(coupling_inter_contrast2_allTrials);
            coupling_intra_display_mean = zMean(zMean(EXT.([physio_sig '_allCouples']).xcorr_cond_intra,3),1);                                          % mean over all subjects and trials
            coupling_intra_contrast1_display_prctl = [zPrctile(coupling_intra_contrast1_allTrials,5); zPrctile(coupling_intra_contrast1_allTrials,95)]; % 5th and 95th percentiles
            coupling_intra_contrast2_display_prctl = [zPrctile(coupling_intra_contrast2_allTrials,5); zPrctile(coupling_intra_contrast2_allTrials,95)];
            coupling_inter_contrast1_display_prctl = [zPrctile(coupling_inter_contrast1_allTrials,5); zPrctile(coupling_inter_contrast1_allTrials,95)];
            coupling_inter_contrast2_display_prctl = [zPrctile(coupling_inter_contrast2_allTrials,5); zPrctile(coupling_inter_contrast2_allTrials,95)];
            figure('Name',['INTRA, effect size on ' coupling_str],'WindowStyle','docked','NumberTitle','off'),
            plot(EXT.coupling_timeAxis_cond, coupling_intra_display_mean, 'color','b','linewidth',2),
            hold on,
            plot(EXT.coupling_timeAxis_cond, coupling_intra_contrast2_display_mean,'color','r','linewidth',2),
            plot(EXT.coupling_timeAxis_cond, coupling_intra_contrast1_display_mean,'color','g','linewidth',2),
            legend('mean',cfgANA.effect_contrast2_str,cfgANA.effect_contrast1_str),
            jbfill(EXT.coupling_timeAxis_cond,coupling_intra_contrast2_display_prctl(2,:),coupling_intra_contrast2_display_prctl(1,:),[255 150 150]./255,'r');
            jbfill(EXT.coupling_timeAxis_cond,coupling_intra_contrast1_display_prctl(2,:),coupling_intra_contrast1_display_prctl(1,:),[150 255 150]./255,'g');
            grid on, ylim([-1 1]), line([0 0],get(gca,'Ylim'),'color','k','linewidth',1)
            title(['INTRA couples: ' physio_sig ' cross-correlation averaged on condition periods']), xlabel('lag (s)'),ylabel('correlation coefficient'),
            figure('Name',['INTER, effect size on ' coupling_str],'WindowStyle','docked','NumberTitle','off'),
            plot(EXT.coupling_timeAxis_cond, coupling_intra_display_mean, 'color','b','linewidth',2),
            hold on,
            plot(EXT.coupling_timeAxis_cond, coupling_inter_contrast2_display_mean,'color','r','linewidth',2),
            plot(EXT.coupling_timeAxis_cond, coupling_inter_contrast1_display_mean,'color','g','linewidth',2),
            legend('mean intra',cfgANA.effect_contrast2_str,cfgANA.effect_contrast1_str),
            jbfill(EXT.coupling_timeAxis_cond,coupling_inter_contrast2_display_prctl(2,:),coupling_inter_contrast2_display_prctl(1,:),[255 150 150]./255,'r');
            jbfill(EXT.coupling_timeAxis_cond,coupling_inter_contrast1_display_prctl(2,:),coupling_inter_contrast1_display_prctl(1,:),[150 255 150]./255,'g');
            grid on, ylim([-1 1]), line([0 0],get(gca,'Ylim'),'color','k','linewidth',1)
            title(['INTER couples: ' physio_sig ' cross-correlation averaged on condition periods']), xlabel('lag (s)'),ylabel('correlation coefficient'),
            clear coupling_intra_contrast1_allTrials coupling_intra_contrast2_allTrials coupling_inter_contrast1_allTrials coupling_inter_contrast2_allTrials
            clear coupling_intra_contrast1_display_mean coupling_intra_contrast2_display_mean coupling_inter_contrast1_display_mean coupling_inter_contrast2_display_mean
            clear coupling_intra_contrast1_display_prctl coupling_intra_contrast2_display_prctl coupling_inter_contrast1_display_prctl coupling_inter_contrast2_display_prctl
        end  
    end
    clear bootSample bootSample1 bootSample2 boot_ix coupling_inter_contrast1 coupling_inter_contrast1_mean coupling_inter_contrast2 coupling_inter_contrast2_mean  
    clear coupling_intra_contrast1 coupling_intra_contrast1_mean coupling_intra_contrast2 coupling_intra_contrast2_mean  coupling_mean_intraNoInter_contrast1 
    clear coupling_mean_intraNoInter_contrast2 effect effect_ix grandMean_intra max_inter max_inter_contrast1 max_inter_contrast1_ix max_inter_contrast2 max_inter_contrast2_ix
    clear max_inter_ix max_intra max_intra_contrast1 max_intra_contrast1_ix max_intra_contrast2 max_intra_contrast2_ix max_intra_ix physio_chan result_str meanFct
    
   
    % ** INTRA/INTER average intersubject coupling: display bar plot with effect size and deviation
    result_str = [physio_sig '_' cfgANA.coupling_measure '_average_IntraInter'];
    coupling_display_mean   = [zMean(coupling_intra_allCond_mean(:)), zMean(coupling_inter_allCond_mean(:))];
    coupling_display_std    = [zStd(coupling_intra_allCond_mean(:)), zStd(coupling_inter_allCond_mean(:))];
    figure('Name',['INTRA vs INTER: average intersubject ' coupling_str],'WindowStyle','docked','NumberTitle','off'),
    barPlot(coupling_display_mean,RESULTS.(result_str).p,[.05 .01 .001],1), hold on, 
    errorb(coupling_display_mean,coupling_display_std,'top');
    set(gca,'XTickLabel',{'intra  ', 'inter  '},'fontsize',10), 
    rotateticklabel(gca,50);  set(gca,'OuterPosition',[0 0.1 1 .9],'FontSize',fontSize);
    ylim([0 1.2]), ylabel(['average ' coupling_str],'FontSize',fontSize), colormap summer,
    clear fontSize result_str coupling_display_mean coupling_display_std coupling_intra_allCond_mean coupling_inter_allCond_mean    
  
    
    % *** EDA derived measure (amplitude of ER-SCR): analysis and display
    if (cfgANA.derived_measures == 1) && strcmp(cfgANA.physio_sig, 'EDA')
        disp('Computing randomization test on amplitude of ER-SCR...');  
        for effect_ix = 1:length(cfgANA.EFFECT_ANALYSIS) 
            effect = cfgANA.EFFECT_ANALYSIS{effect_ix}; 
            amp_contrast1 = squeeze(EXT.EDA_allCouples.AMP_cond(effect_ix,1,:));
            amp_contrast2 = squeeze(EXT.EDA_allCouples.AMP_cond(effect_ix,2,:));
            [RESULTS.EDA_AMP_stats.(cfgANA.EFFECT_ANALYSIS{effect_ix}).h, RESULTS.EDA_AMP_stats.(cfgANA.EFFECT_ANALYSIS{effect_ix}).p] = ...
                    pairedPermTest(amp_contrast1,amp_contrast2,cfgANA.tTest_alpha,cfgANA.nbPerms);
            EDA_AMP_pVal_display(effect_ix) = RESULTS.EDA_AMP_stats.(cfgANA.EFFECT_ANALYSIS{effect_ix}).p;
            disp(['Effect of ' effect ' on amplitude of ER-SCR: p=' num2str(EDA_AMP_pVal_display(effect_ix))]);
        end   
        % display mean EDA AMP for 2 contrasts 
        EDA_AMP_mean_display    = mean(EXT.EDA_allCouples.AMP_cond,3);
        EDA_AMP_std_display     = std(EXT.EDA_allCouples.AMP_cond,[],3);
        figure('Name','EDA, amplitude of ER-SCR','WindowStyle','docked','NumberTitle','off'),
        barPlot(EDA_AMP_mean_display,EDA_AMP_pVal_display,[.05 .01 .001],.8), hold on, 
        errorb(EDA_AMP_mean_display,EDA_AMP_std_display,'top');
        xtick = sort([(1:length(cfgANA.EFFECT_ANALYSIS))-.2,(1:length(cfgANA.EFFECT_ANALYSIS))+.2]);
        set(gca,'Xtick',xtick,'XTickLabel',cfgANA.all_effect_str,'fontsize',8), 
        rotateticklabel(gca,50);  set(gca,'OuterPosition',[0 0.3 1 .7]);
        ylim([0 .9]), ylabel('ER-SCR normalized amplitude'), colormap summer,
        clear EDA_AMP_mean_display EDA_AMP_std_display EDA_AMP_pVal_display effect_ix amp_contrast1 amp_contrast2
    end
    
    % *** PRV derived measure (MEAN, SD, LF, HF): analysis and display
    if (cfgANA.derived_measures == 1) && strcmp(cfgANA.physio_sig, 'PRV')
        disp('Computing randomization test on PRV derived measures...');  
        PRV_measure_str = {'MEAN','SD','LF','HF'};
        for measure_ix = 1:length(PRV_measure_str)
            curr_meas = PRV_measure_str{measure_ix};
            for effect_ix = 1:length(cfgANA.EFFECT_ANALYSIS) 
                effect = cfgANA.EFFECT_ANALYSIS{effect_ix}; 
                measure_contrast1 = squeeze(EXT.PRV_allCouples.PRV_measures_effect.(curr_meas)(effect_ix,1,:));
                measure_contrast2 = squeeze(EXT.PRV_allCouples.PRV_measures_effect.(curr_meas)(effect_ix,2,:));
                [RESULTS.(['PRV_stats_' curr_meas]).(cfgANA.EFFECT_ANALYSIS{effect_ix}).h, RESULTS.(['PRV_stats_' curr_meas]).(cfgANA.EFFECT_ANALYSIS{effect_ix}).p] = ...
                        pairedPermTest(measure_contrast1,measure_contrast2,cfgANA.tTest_alpha,cfgANA.nbPerms);
                pVal_display(effect_ix) = RESULTS.(['PRV_stats_' curr_meas]).(cfgANA.EFFECT_ANALYSIS{effect_ix}).p;
                disp(['Effect of ' effect ' on PRV ' curr_meas ': p=' num2str(pVal_display(effect_ix))]);
            end   
            % display average derived PRV measure for both contrasts 
            measure_mean_display    = mean(EXT.PRV_allCouples.PRV_measures_effect.(curr_meas),3);
            measure_std_display     = std(EXT.PRV_allCouples.PRV_measures_effect.(curr_meas),[],3);
            figure('Name',['PRV, ' curr_meas],'WindowStyle','docked','NumberTitle','off'),
            barPlot(measure_mean_display,pVal_display,[.05 .01 .001]), hold on, 
            errorb(measure_mean_display,measure_std_display,'top');
            xtick = sort([(1:length(cfgANA.EFFECT_ANALYSIS))-.2,(1:length(cfgANA.EFFECT_ANALYSIS))+.2]);
            set(gca,'Xtick',xtick,'XTickLabel',cfgANA.all_effect_str,'fontsize',8), 
            rotateticklabel(gca,50);  set(gca,'OuterPosition',[0 0.3 1 .7]);
            ylabel(['PRV ' curr_meas]), colormap summer,
        end
        clear measure_str measure_ix curr_meas effect_ix effect measure_contrast1 measure_contrast2 pVal_display measure_mean_display measure_std_display xtick
    end
    
    
    % *** INTRA COUPLE: VARIOUS DISPLAY ***
    for couple_ix = 1:cfgEEG.nbCouple
        % ** display physio measure of both subjects for whole experiment
        if (cfgDISP.physio_allExpe == 1)
            figure('name', ['Couple n°' int2str(couple_ix) ' ' physio_sig ' on all experiment'],'WindowStyle','docked','NumberTitle','off'),
            sigToPlot = [EXT.([physio_sig '_allCouples']).all(2*couple_ix-1,:); EXT.([physio_sig '_allCouples']).all(2*couple_ix,:)];
            if (strcmp(cfgANA.coupling_measure,'PLV')==1) % display PLV
                sigToPlot = [sigToPlot ; (EXT.([physio_sig '_allCouples']).PLV_all_intra(couple_ix,:) -.5)]; % HERE we remove .5 to PLV for display ease
            end
            color_rest = [100 100 100]./255; %[120 120 120]./255;
            color_unused = [200 200 200]./255; % [200 200 200]./255; [51 204 0]
            color_run   = [0 51 0]./255;
            yPos = -1; h = 2.2;
            % draw rectangles for rest periods
            for rest_ix = 1:size(TRIG.trigSub(1,1).restVal,1)
                xPos = TRIG.trigSub(1,1).restVal(rest_ix,1) ./ cfgEEG.Fs;   
                w = (TRIG.trigSub(1,1).restVal(rest_ix,2) - TRIG.trigSub(1,1).restVal(rest_ix,1)) ./ cfgEEG.Fs;
                rectangle('Position',[xPos,yPos,w,h],'FaceColor',color_rest,'EdgeColor',color_rest); hold on, 
            end
            % draw rectangles for unused and used periods within runs 
            if (cfgANA.trialLim(1)~=0)
                for run_ix = 1:cfgEEG.nbCondTot
                    xPos = (TRIG.trigSub(1,1).restVal(run_ix,2)/cfgEEG.Fs) - cfgANA.trialLim(1);   
                    w = cfgANA.trialLim(1);
                    rectangle('Position',[xPos,yPos,w,h],'FaceColor',color_unused,'EdgeColor',color_unused); hold on, 
                    xPos = TRIG.trigSub(1,1).restVal(run_ix,2)/cfgEEG.Fs;   
                    w = cfgANA.trialLim(2)-cfgANA.trialLim(1);
                    rectangle('Position',[xPos,yPos,w,h],'FaceColor','w','EdgeColor','w'); hold on, 
                end
            end
            xlim([0 max(cfgEEG.t_axis_allExpe)]), 
            ylim(1*[min(min(sigToPlot)), max(max(sigToPlot))]),
            % display run labels
            for run_ix = 1:cfgEEG.nbCondTot
                xPos = TRIG.trigSub(1,1).restVal(run_ix,2)/cfgEEG.Fs;
                yPos = .98*max(get(gca,'YLim'));
                text(xPos,yPos,TRIG.trigSub(2*couple_ix-1).dispTrig2{run_ix},'FontSize',6);%,'FontWeight','bold');
            end
            set(gca,'layer','top'),
            plot(cfgEEG.t_axis_allExpe,sigToPlot(1,:),'linewidth',1.5,'color','k'), % [1 .4 0]
            plot(cfgEEG.t_axis_allExpe,sigToPlot(2,:),'linewidth',1.5,'color','b'), % [0 .5 .4]
            legend('subject #1','subject #2'),
            if (strcmp(cfgANA.coupling_measure,'PLV')==1) % display PLV
                plot(cfgEEG.t_axis_allExpe,sigToPlot(3,:),'linewidth',1.5,'color','r'), % [0 .5 .4]
                legend('subject #1','subject #2','PLV (-.5)'),
            end
            xlabel('time (s)'),ylabel(cfgANA.physio_str2),
            title( ['Couple n°' int2str(couple_ix) ]), 
        end
        clear h rest_ix run_ix sigToPlot w xPos yPos color_rest color_run color_unused couple_ix

        % ** display physiological variability trend for all experiment
        if (cfgDISP.physioTrend_allExpe ==1)
            figure('name', ['Couple n°' int2str(couple_ix) ' ' physio_sig ' TREND on all experiment'],'WindowStyle','docked','NumberTitle','off'),
            sigToPlot = [EXT.([physio_sig '_allCouples']).all_trend(2*couple_ix-1,:); EXT.([physio_sig '_allCouples']).all_trend(2*couple_ix,:)];
            plot(cfgEEG.t_axis_allExpe, sigToPlot(1,:),'color','b','linewidth',2),
            hold on,
            plot(cfgEEG.t_axis_allExpe, sigToPlot(2,:),'color','r','linewidth',2),
            xlim([0 1500]), grid on, ylabel(cfgANA.physio_str2), %xlabel('time (s)'), % legend('subject 1','subject 2'),
            title( ['Couple n°' int2str(couple_ix) ' ' physio_sig ' trend']),
        end

        % ** display physio auto- and cross-correlation on whole experiment
        if (cfgDISP.xcorr_sub_allExpe == 1) && (strcmp(cfgANA.coupling_measure,'xcorr'))
            figure('Name',['Couple n°' int2str(couple_ix) ' ' physio_sig ' total xcorr'],'WindowStyle','docked','NumberTitle','off'),
            sigToPlot = [EXT.([physio_sig '_allCouples']).xcorr_all_intra(couple_ix,:); EXT.([physio_sig '_allCouples']).autocorr_all(2*couple_ix-1,:);  EXT.([physio_sig '_allCouples']).autocorr_all(2*couple_ix,:)];
            plot(EXT.coupling_timeAxis_all, sigToPlot(1,:), 'linewidth',2),
            hold on,
            plot(EXT.coupling_timeAxis_all, sigToPlot(2,:),'g'),
            plot(EXT.coupling_timeAxis_all, sigToPlot(3,:),'c'),
            legend('xcorr 1-2','autocorr 1','autocorr 2'),
            line([0 0],[-1 1],'color','r','linewidth',1),
            xlim([-1500 1500]), grid on,  ylim([-1 1]),  % ylim([-.5 .5]),
            xlabel('lag (s)'), ylabel('correlation coefficient'),
            title(['Couple n°' int2str(couple_ix) ', ' physio_sig ' cross-correlation on whole experiment']), 
        end
        
    end

  	% *** INTER COUPLE: VARIOUS DISPLAY ***
    
   
end   % end procedure B/
    
    
   

%% *** A/ Group BSS and analysis of component spectra
if cfgANA.procedureA_groupBSS == 1

if strcmp(cfgANA.dataType,'EXT_only')
    error('...to perform group BSS on EEG, dataType must comprise EEG!...')
end
    
% *** TEMP DEBUG
DEBUG = 0;
if ~DEBUG
% *** TEMP DEBUG
    
% --- A1) Estimate average ICA components for all subjects -> B_global
% --- Process average cospectra across subjects
% this should be done in a new loop, since we don't know each subject's total EEG power in advance!
nbConcat_contrast1      = 0; % this index increments for each new cospectrum added to cosp_global
nbConcat_contrast2      = 0; % this index increments for each new cospectrum added to cosp_global2
nbLoops                 = cfgEEG.nbSub*cfgEEG.nbCouple*cfgEEG.nbCondTot; % for UI purpose only
currUI_ix               = 0;  % index for UI display only
cfgEEG.cospWinSize      = 2^nextpow2(cfgEEG.Fs/cfgEEG.FreqRes);
cfgEEG.Faxis            = [cfgEEG.FreqRangeCosp(1): cfgEEG.FreqRes : cfgEEG.FreqRangeCosp(2)];
cfgEEG.nbFreq           = length(cfgEEG.Faxis);
EEG.cosp_global = zeros(cfgEEG.nbChanSub,cfgEEG.nbChanSub, cfgEEG.nbFreq);  % memory allocation
if cfgBSS.jDiagContrast == 1
    temp_cosp_global2   = zeros(cfgEEG.nbChanSub,cfgEEG.nbChanSub, cfgEEG.nbFreq); % create another cospectra set for concatenation with contrast conditions
    nbConcat2           = cfgEEG.nbSub*cfgEEG.nbCouple*length(cfgANA.effect_contrast2); % 2 x (nb couples) x (nb conditions) 
end
disp('Cospectra estimation...');

for couple_ix = cfgANA.coupleID  
    for sub_ix = 1:cfgEEG.nbSub  
        sub_ix_total = cfgEEG.nbSub*(find(ismember(cfgANA.coupleID,couple_ix))-1)+sub_ix; % this index increments for each new subject
        
        % compute cospectra for all conditions
        for cond_ix = 1:cfgEEG.nbCondTot
            % UI display
            currUI_ix = currUI_ix + 1;
            prog = fix(100*currUI_ix/nbLoops);
            if mod(prog,10) == 0
                disp(['progression: ' int2str(prog) '%...']); 
            end
            
            % select data on which to process cospectra
            curr_cond       = cfgEEG.CondTrigVal(cond_ix);
            curr_cond_str   = ['cond_' int2str(curr_cond)];
            dataCurrentTrial = EEG.data_sub(sub_ix_total).(curr_cond_str); 
   
            % data standardization ACROSS COUPLES / IN EACH CONDITION
            if strcmp(cfgANA.standardizeEEG,'subject')
                trialMeanPow = mean(std(double(dataCurrentTrial),0,2));
                dataCurrentTrial = dataCurrentTrial ./ ((trialMeanPow/mean(EEG.meanPow))*ones(cfgEEG.S(sub_ix),size(dataCurrentTrial,2))); 
            end

            % Compute real part of cross-spectra only (i.e. cospectra),
            EEG.cosp_sub(sub_ix_total).(curr_cond_str) = real(procCosp( ...
                dataCurrentTrial', cfgEEG.WinType, cfgEEG.cospWinSize, ...
                cfgEEG.WinOverlap, [cfgEEG.Fs cfgEEG.FreqRangeCosp])); 

            % cospectra standardization for each subject / condition / FREQUENCY
            if strcmp(cfgANA.standardizeEEG,'frequency')           
                % normalization to have AVERAGED cosp over frequency with unitary diagonal
                inFreqs = find((cfgEEG.Faxis>=cfgBSS.FreqRangeBSS(1))&(cfgEEG.Faxis<=cfgBSS.FreqRangeBSS(2)));
                meanCosp = mean(EEG.cosp_sub(sub_ix_total).(curr_cond_str)(:,:,inFreqs),3);  
                for f_ix = 1:cfgEEG.nbFreq
%                  	tempCosp = EEG.cosp_sub(sub_ix_total).(curr_cond_str)(:,:,f_ix);
                    EEG.cosp_sub(sub_ix_total).(curr_cond_str)(:,:,f_ix) = ...
                        diag(diag(meanCosp.^(-1/2)))*EEG.cosp_sub(sub_ix_total).(curr_cond_str)(:,:,f_ix)*diag(diag(meanCosp.^(-1/2))); 
                end
            end

            % AVERAGING COSPECTRA over subjects and/or conditions of interest, 
            curr_sub_ix = 2*couple_ix-2+sub_ix;   % subject index in the experiment (e.g. couple 9, sub 1 -> curr_sub_ix = 17)
            % average only if right condition, subject role (empathizer, or transmitter, or both), and result to subjective reports  
            if ismember(curr_cond,cfgANA.cond_contrast1) & ismember(curr_cond,cfgEEG.SR_select_sub{sub_ix_total}) & ismember(curr_sub_ix,cfgANA.sub_contrast1);
                EEG.cosp_global = EEG.cosp_global + EEG.cosp_sub(sub_ix_total).(curr_cond_str);          
                nbConcat_contrast1 = nbConcat_contrast1 + 1;
            elseif (cfgBSS.jDiagContrast==1) & ismember(curr_cond,cfgANA.cond_contrast2) & ismember(curr_cond,cfgEEG.SR_select_sub{sub_ix_total}) & ismember(curr_sub_ix,cfgANA.sub_contrast2)
                temp_cosp_global2 = temp_cosp_global2 + EEG.cosp_sub(sub_ix_total).(curr_cond_str);
                nbConcat_contrast2 = nbConcat_contrast2 + 1;
            end
        end
    end
end

% --- Normalize averaged global cospectrum with correct number added conditions
EEG.cosp_global = (1/nbConcat_contrast1) * EEG.cosp_global;
if (cfgBSS.jDiagContrast==1)
    temp_cosp_global2 = (1/nbConcat_contrast2) * temp_cosp_global2;
end

% --- Spectrum calculation as diag part of cospectra
EEG.spectrum_global = zeros(cfgEEG.nbChanSub, cfgEEG.nbFreq); % here only for conditions from first contrast 
for f_ix = 1 : cfgEEG.nbFreq
   EEG.spectrum_global(:,f_ix) = diag(EEG.cosp_global(:,:,f_ix));
end
clear dataCurrentTrial curr_cond curr_cond_str currUI_ix meanCosp nbLoops nbConcat_contrast1 nbConcat_contrast2
clear couple_ix_total couple_ix effect_ix f_ix sub_ix sub_ix_total Fmin_ix Fmax_ix curr_sub_ix trialMeanPow   

% --- Concatenation of cospectra when simultaneously diagonalize cospectra from both contrast conditions  
if cfgBSS.jDiagContrast == 1
    EEG.cosp_global = cat(3,EEG.cosp_global,temp_cosp_global2);
    clear temp_cosp_global2
end

% --- Weight cospectra with non-diagonality function
if (cfgBSS.ApplyNonDiagFct == 1)
    EEG.diagFct = NonDiagonality(EEG.cosp_global);
    for freq_ix = 1 : size(EEG.cosp_global,3)
        EEG.cosp_global(:,:,freq_ix) = EEG.cosp_global(:,:,freq_ix).* EEG.diagFct(freq_ix);
    end
    disp('Cospectra weighted with non-diagonality function.');
end

% --- Withening
if (cfgBSS.ApplyWhitening == 1)
    % Compute a whitening matrix H on sum of cospectra in given frequency range
    inFreqs = find((cfgEEG.Faxis>=cfgBSS.FreqRangeBSS(1))&(cfgEEG.Faxis<=cfgBSS.FreqRangeBSS(2)));
    if cfgBSS.jDiagContrast == 1    % in case two contrast conditions are concatenated
        inFreqs = [inFreqs (inFreqs+cfgEEG.nbFreq)];
    end
    [EEG.H_global, EEG.Q_global, cfgEEG.nbIC, EEG.eigValsW] =  simpleWhitening(sum(EEG.cosp_global(:,:,inFreqs),3),cfgBSS.PercentVariance);	% sum over freq range 
    cfgEEG.ICset   = [1:cfgEEG.nbIC ; cfgEEG.nbIC+1:2*cfgEEG.nbIC];  % subject1 = 1st row, subject2 = 2nd row
    
    % Whithening cospectra
    for freq_ix = 1 : size(EEG.cosp_global,3)
        EEG.cospW_global(:,:,freq_ix) = EEG.H_global * EEG.cosp_global(:,:,freq_ix) * EEG.H_global'; 
    end
    disp(['Cospectra whitened in range ' int2str(cfgBSS.FreqRangeBSS(1)) ' to ' int2str(cfgBSS.FreqRangeBSS(2)) 'Hz.']);
    clear Fmin_ix Fmax_ix freq_ix
else
    disp('No data pre-whitening.');
    EEG.cospW_global    = EEG.cosp_global;
    cfgEEG.nbIC         = cfgEEG.nbChanSub;
end


% --- Process B_global using U-WEDGE 
inFreqs = (cfgEEG.Faxis>=cfgBSS.FreqRangeBSS(1))&(cfgEEG.Faxis<=cfgBSS.FreqRangeBSS(2));
cospW_global_BSS = EEG.cospW_global(:,:,inFreqs);
EEG.B_global = uwedge(cospW_global_BSS(:,:));  % BSS using UWEDGE algorithm
disp('B_global estimated using U-WEDGE algorithm.');                    
% Calculation of corresponding mixing matrix A as pseudo-inverse of B
EEG.A_global = pinv(EEG.B_global); 
% When cospectra were pre-whithened, undo whithening transform
if (cfgBSS.ApplyWhitening == 1)
    EEG.A_global = EEG.Q_global * EEG.A_global;
    EEG.B_global = EEG.B_global * EEG.H_global; 
end
% Normalization of mixing and demixing matrix 
EEG.A_global = EEG.A_global * diag(sum( EEG.A_global.^2,1).^(-1/2));     % on the columns
EEG.B_global = diag(sum( EEG.B_global.^2,2).^-(1/2)) *  EEG.B_global;    % on the rows
clear cospW_global_BSS inFreqs


% ---  Compute average ICA components for all subjects
% Processing source cospectra and spectra
for f_ix = 1 : cfgEEG.nbFreq
    EEG.cospIC_global(:,:,f_ix) = EEG.B_global * EEG.cosp_global(:,:,f_ix) * EEG.B_global';
end
inFreqs = (cfgEEG.Faxis>=cfgEEG.FreqRangeNorm(1))&(cfgEEG.Faxis<=cfgEEG.FreqRangeNorm(2));
% Normalization at each frequency on the rows and on the columns so as to have 'cospIC_mean' with unitary main diagonal
% meanCosp = mean(EEG.cospW_global(:,:,inFreqs),3); % normalization within frequency range used for JBSS
% for f_ix = 1 : cfgEEG.nbFreq
%     EEG.cospW_global(:,:,f_ix) = diag(diag(meanCosp.^(-1/2))) * EEG.cospW_global(:,:,f_ix) * diag(diag(meanCosp.^(-1/2)));
% end
for IC_ix = 1 : cfgEEG.nbIC
    EEG.spectrumIC_global(IC_ix,:) = EEG.cospIC_global(IC_ix,IC_ix,:);
    EEG.spectrumIC_global(IC_ix,:) = EEG.spectrumIC_global(IC_ix,:) ./ sum(EEG.spectrumIC_global(IC_ix,inFreqs));    
end

% Source separation on EEG and cospectra for each subject / effect
for couple_ix = cfgANA.coupleID  
    for sub_ix = 1:cfgEEG.nbSub  
        sub_ix_total = cfgEEG.nbSub*(find(ismember(cfgANA.coupleID,couple_ix))-1)+sub_ix; % this index increments for each new subject
        for cond_ix = 1:cfgEEG.nbCondTot
            curr_cond       = cfgEEG.CondTrigVal(cond_ix);
            curr_cond_str   = ['cond_' int2str(curr_cond)];
            % source separation on EEG
            EEG.data_sub_IC(sub_ix_total).(curr_cond_str) = EEG.B_global * EEG.data_sub(sub_ix_total).(curr_cond_str);
            % source separation on cospectra for each subject / effect
            for f_ix = 1 : cfgEEG.nbFreq
               EEG.cospIC_sub(sub_ix_total).(curr_cond_str)(:,:,f_ix) = EEG.B_global * EEG.cosp_sub(sub_ix_total).(curr_cond_str)(:,:,f_ix) * EEG.B_global';
               EEG.spectrumIC_sub(sub_ix_total).(curr_cond_str)(:,f_ix) = diag(EEG.cospIC_sub(sub_ix_total).(curr_cond_str)(:,:,f_ix));
            end
        end
    end
end
disp('BSS processed on EEG and (co)spectra.');     
clear inFreqs f_ix IC_ix effect_ix curr_cond_str curr_cond sub_ix sub_ix_total couple_ix prog

% TEMP FOR DEBUG (if debug set to 0, save it, else: just load it)
disp('Saving temporary workspace...');
save TempWorkspace
disp('Temporary workspace saved.');




% --- A2) Visualize components information and put IC of interest in first positions
% --- Display EEG averaged cospectra
figure('Name','EEG averaged cospectra','WindowStyle','docked'), imagescMatrixSelection(EEG.cosp_global);


% --- Display Non-diagonality function 
if cfgBSS.ApplyNonDiagFct == 1
    if cfgBSS.jDiagContrast == 0
        figure('Name','Non-diagonality function (global)','WindowStyle','docked'), 
        plot(cfgEEG.Faxis,EEG.diagFct), grid on, xlabel('Frequency (Hz)'),
    else
        figure('Name','Non-diagonality function (global)','WindowStyle','docked'), 
        plot(cfgEEG.Faxis,EEG.diagFct(1:cfgEEG.nbFreq)), grid on, xlabel('Frequency (Hz)'), hold on,
        plot(cfgEEG.Faxis,EEG.diagFct(cfgEEG.nbFreq+1:end),'r'), legend(cfgANA.effect_contrast1_str,cfgANA.effect_contrast2_str),  
    end
end

% --- Display EEG averaged cospectra AFTER whitening and eigenvalues
if (cfgBSS.ApplyWhitening == 1) && (cfgBSS.ApplyNonDiagFct == 1) 
    figure('Name','EEG averaged cospectra AFTER whitening','WindowStyle','docked'), imagescMatrixSelection(EEG.cospW_global);
    figure('Name','Components eigenvalues','WindowStyle','docked'), 
    plot(EEG.eigValsW), grid on, xlabel('IC#'),ylabel('eigen values'), xlim([1 cfgEEG.nbChanSub]),
    line(get(gca,'Xlim'),[EEG.eigValsW(cfgEEG.nbIC) EEG.eigValsW(cfgEEG.nbIC)],'Color','r','linewidth',1,'linestyle','-');
    line([cfgEEG.nbIC cfgEEG.nbIC],get(gca,'Ylim'),'Color','r','linewidth',1,'linestyle','-');
end

% --- Display global components cospectra
figure('Name','IC cospectra','WindowStyle','docked'), imagescMatrixSelection(EEG.cospIC_global);
inFreqs = find((cfgEEG.Faxis>=cfgBSS.FreqRangeBSS(1))&(cfgEEG.Faxis<=cfgBSS.FreqRangeBSS(2)));
figure('Name',['IC mean cospectrum [' int2str(cfgBSS.FreqRangeBSS(1)) '-' int2str(cfgBSS.FreqRangeBSS(2)) 'Hz]'],...
        'WindowStyle','docked'), imagescMatrixSelection(mean(EEG.cospIC_global(:,:,inFreqs),3));
clear inFreqs
% --- Display global components spectra and scalp maps
% nbLoc       = cfgEEG.nbIC;     	% number of scalp maps to plot for each subject
% Fdisp       = [2 30];           % frequency range for spectrum display
% CompScalpMapSpectra(cfgEEG.Faxis,nbLoc,Fdisp,cfgEEG.eloc,EEG.A_global,EEG.spectrumIC_global);    


% --- Display contribution of global components at specified frequencies
%(TODO)


% --- Display sLORETA localization from global mixing matrix
str.transfMatFile = '\JOJOS MATLAB TOOLBOX\Scripts\EEG_20_sLoreta_alpha0.tm';
str.powFile = '\JOJOS MATLAB TOOLBOX\Scripts\EEG_20_subGlobal.lor';
% select mixing matrices and remove mean from each column 
A_loreta = EEG.A_global - repmat(mean(EEG.A_global),cfgEEG.nbChanSub,1);
% sLORETA source localization
sLoreta_loc2(A_loreta, str.transfMatFile, str.powFile); 
clear str A_loreta

% *** TEMP DEBUG
else
    cfgANA_temp = cfgANA;
    disp('Loading temporary workspace...');
    load TempWorkspace.mat
    disp('Temporary workspace loaded.');
    cfgANA = cfgANA_temp; 
    clear cfgANA_temp 
end
% *** TEMP DEBUG




% --- A3) Compute t-test on components spectra estimated in different conditions
% ---  Compute average IC spectra for each effect in contrasts
% 'RESULTS.tTest_freqs_index' contains indexes of frequencies to use for t-test
%       format: (freqToAverage,tTest_nbFreqs)
RESULTS.tTest_freqs_index = find((cfgEEG.Faxis>=cfgANA.tTest_freqRange(1))&(cfgEEG.Faxis<=cfgANA.tTest_freqRange(2)));
RESULTS.tTest_nbFreqs = length(RESULTS.tTest_freqs_index);
if(cfgANA.tTest_freqRes ~= cfgEEG.FreqRes)
    freqAverageFactor       = fix(cfgANA.tTest_freqRes/cfgEEG.FreqRes);
    RESULTS.tTest_nbFreqs   = fix(RESULTS.tTest_nbFreqs/freqAverageFactor);
    RESULTS.tTest_freqs_index = reshape(RESULTS.tTest_freqs_index(1:RESULTS.tTest_nbFreqs*freqAverageFactor),freqAverageFactor,RESULTS.tTest_nbFreqs); 
end
RESULTS.tTest_freqs_axis = mean(cfgEEG.Faxis(RESULTS.tTest_freqs_index),1); % (1,tTest_nbFreqs) vector of frequencies in Hz

spectrumIC_sub_contrast1 = zeros(cfgEEG.nbIC,cfgEEG.nbFreq,cfgEEG.nbSubTot);
spectrumIC_sub_contrast2 = zeros(cfgEEG.nbIC,cfgEEG.nbFreq,cfgEEG.nbSubTot);
for couple_ix = cfgANA.coupleID  
    for sub_ix = 1:cfgEEG.nbSub  
        sub_ix_total = cfgEEG.nbSub*(find(ismember(cfgANA.coupleID,couple_ix))-1)+sub_ix; % this index increments for each new subject
        nbConcat_contrast1 = 0;
        nbConcat_contrast2 = 0;        
        for cond_ix = cfgANA.cond_contrast1
            curr_cond_str = ['cond_' int2str(cond_ix)];
            if ismember(cond_ix,cfgEEG.SR_select_sub{sub_ix_total}) % consider self report on emotion intensity when adding spectra from current condition
                spectrumIC_sub_contrast1(:,:,sub_ix_total) = spectrumIC_sub_contrast1(:,:,sub_ix_total) + EEG.spectrumIC_sub(sub_ix_total).(curr_cond_str);
                nbConcat_contrast1 = nbConcat_contrast1 + 1;
            end
        end
        spectrumIC_sub_contrast1(:,:,sub_ix_total) = (1/nbConcat_contrast1)*spectrumIC_sub_contrast1(:,:,sub_ix_total);
        for cond_ix = cfgANA.cond_contrast2
            curr_cond_str = ['cond_' int2str(cond_ix)];
            if ismember(cond_ix,cfgEEG.SR_select_sub{sub_ix_total}) % consider self report on emotion intensity when adding spectra from current condition
                spectrumIC_sub_contrast2(:,:,sub_ix_total) = spectrumIC_sub_contrast2(:,:,sub_ix_total) + EEG.spectrumIC_sub(sub_ix_total).(curr_cond_str);
                nbConcat_contrast2 = nbConcat_contrast2 + 1;
            end
        end
        spectrumIC_sub_contrast2(:,:,sub_ix_total) = (1/nbConcat_contrast2)*spectrumIC_sub_contrast2(:,:,sub_ix_total);
    end
end

% -- Power spectrum frequency average and standardization
EEG.spectrumIC_sub_contrast1 = zeros(cfgEEG.nbIC,RESULTS.tTest_nbFreqs,cfgEEG.nbSubTot);
EEG.spectrumIC_sub_contrast2 = zeros(cfgEEG.nbIC,RESULTS.tTest_nbFreqs,cfgEEG.nbSubTot);
for couple_ix = cfgANA.coupleID  
    for sub_ix = 1:cfgEEG.nbSub  
        sub_ix_total = cfgEEG.nbSub*(find(ismember(cfgANA.coupleID,couple_ix))-1)+sub_ix; % this index increments for each new subject
        
        % -- Frequency average (when less freqs tested than formerly computed)
        for f_ix = 1:RESULTS.tTest_nbFreqs   
            EEG.spectrumIC_sub_contrast1(:,f_ix,sub_ix_total) = mean(squeeze(spectrumIC_sub_contrast1(:,RESULTS.tTest_freqs_index(:,f_ix),sub_ix_total)),2);
            EEG.spectrumIC_sub_contrast2(:,f_ix,sub_ix_total) = mean(squeeze(spectrumIC_sub_contrast2(:,RESULTS.tTest_freqs_index(:,f_ix),sub_ix_total)),2);
        end
        
        % -- Standardization: "normal" power spectrum
        % --> i.e. each IC is normalized by its total power in considered frequency band
        if strcmp(cfgANA.tTest_powType,'normal') 
            for f_ix = 1:RESULTS.tTest_nbFreqs   
                EEG.spectrumIC_sub_contrast1(:,f_ix,sub_ix_total) = EEG.spectrumIC_sub_contrast1(:,f_ix,sub_ix_total) ./ sum(squeeze(EEG.spectrumIC_sub_contrast1(:,:,sub_ix_total)),2);
                EEG.spectrumIC_sub_contrast2(:,f_ix,sub_ix_total) = EEG.spectrumIC_sub_contrast2(:,f_ix,sub_ix_total) ./ sum(squeeze(EEG.spectrumIC_sub_contrast2(:,:,sub_ix_total)),2);
            end
        % -- Standardization: "relative" power spectrum
        % --> i.e. each IC is normalized by the total power among all ICs in that discrete frequency
        elseif strcmp(cfgANA.tTest_powType,'relative')     
            for f_ix = 1:RESULTS.tTest_nbFreqs   
                EEG.spectrumIC_sub_contrast1(:,f_ix,sub_ix_total) = EEG.spectrumIC_sub_contrast1(:,f_ix,sub_ix_total) ./ sum(squeeze(EEG.spectrumIC_sub_contrast1(:,f_ix,sub_ix_total)),1);
                EEG.spectrumIC_sub_contrast2(:,f_ix,sub_ix_total) = EEG.spectrumIC_sub_contrast2(:,f_ix,sub_ix_total) ./ sum(squeeze(EEG.spectrumIC_sub_contrast2(:,f_ix,sub_ix_total)),1);
            end
        end
    end
end

% -- Power spectrum average over subjects
EEG.spectrumIC_global_constrast1 = mean(EEG.spectrumIC_sub_contrast1,3);
EEG.spectrumIC_global_constrast2 = mean(EEG.spectrumIC_sub_contrast2,3);
clear nbEffect_global couple_ix sub_ix sub_ix_total effect_ix curr_effect spectrumIC_sub_contrast1 spectrumIC_sub_contrast2


% --- Compute t-test on components spectra at each frequency
RESULTS.tTest_sig    	= zeros(cfgEEG.nbIC,RESULTS.tTest_nbFreqs);        % for each IC x Freq, permutation test significance 
RESULTS.tTest_p_OBS    	= zeros(cfgEEG.nbIC,RESULTS.tTest_nbFreqs);        % for each IC x Freq, observed p values    
RESULTS.tTest_pmin_PERM = zeros(cfgEEG.nbIC,cfgANA.nbPerms);  	% for each IC x Perm, p-value minimum across frequencies
tTest_p_PERM_temp       = zeros(cfgEEG.nbIC,RESULTS.tTest_nbFreqs,cfgANA.nbPerms); % temporary array
RESULTS.tTest_prctile_contrast1	= zeros(cfgEEG.nbIC,RESULTS.tTest_nbFreqs,2);        % for each IC x Freq, 5% and 95% percentile across subjects
RESULTS.tTest_prctile_contrast2	= zeros(cfgEEG.nbIC,RESULTS.tTest_nbFreqs,2);            % for each IC x Freq, 5% and 95% percentile across subjects

% -- create shuffling matrix
cfgEEG.nbSubToPermute   = min([length(cfgANA.sub_contrast1) length(cfgANA.sub_contrast2)]);
RESULTS.permMat         = uniqueShuffle(cfgANA.nbPerms,cfgEEG.nbSubToPermute);  % permutation matrix used to shuffle subjects
cfgANA.nbPerms          = size(RESULTS.permMat,1);              

% -- perform correction of t-test significance level according to number of frequencies and IC tested
if cfgANA.tTest_doAlphaCorr == 1
    cfgANA.tTest_alpha = cfgANA.tTest_alpha / cfgEEG.nbIC;  % Bonferroni correction
end
if cfgANA.tTest_alpha < (1/cfgANA.nbPerms)
   disp(['Not enough permutations for this significance level (' num2str(cfgANA.tTest_alpha) ')']);
   cfgANA.tTest_alpha = 1/cfgANA.nbPerms;
   disp(['--> Significance level set to minimum possible (' num2str(cfgANA.tTest_alpha) ')']);
end


% ********* TEMP FOR DEBUG (this allow testing permutation procedure)
% mu1 = 10;
% mu2 = 20;
% sd  = 1:cfgEEG.nbIC;  % sd is increasing for each IC
% 
% for couple_ix = cfgANA.coupleID  
%     for sub_ix = 1:cfgEEG.nbSub  
%         sub_ix_total = cfgEEG.nbSub*(find(ismember(cfgANA.coupleID,couple_ix))-1)+sub_ix; % this index increments for each new subject
%         EEG.spectrumIC_sub_contrast1(:,:,sub_ix_total) = mu1 + (sd'*ones(1,cfgEEG.nbFreq)).*randn(cfgEEG.nbIC,cfgEEG.nbFreq);
%         EEG.spectrumIC_sub_contrast2(:,:,sub_ix_total) = mu2 + (sd'*ones(1,cfgEEG.nbFreq)).*randn(cfgEEG.nbIC,cfgEEG.nbFreq);
%     end
% end
% EEG.spectrumIC_global_constrast1 = mean(EEG.spectrumIC_sub_contrast1,3);
% EEG.spectrumIC_global_constrast2 = mean(EEG.spectrumIC_sub_contrast2,3);
% *************************


for f_ix = 1:RESULTS.tTest_nbFreqs   
    
    % -- UI display
    prog = fix(100*f_ix/RESULTS.tTest_nbFreqs);
    if mod(prog,10) == 0
        disp(['Processing t-test with permutations (' int2str(prog) '%)...']); 
    end
    
    % -- Keep only subjects with specified role for each contrast
    IC_FreqVal_contrast1    = squeeze(EEG.spectrumIC_sub_contrast1(:,f_ix,ismember(cfgEEG.subRoles(1,:),cfgANA.sub_contrast1)))'; % format: nbSub x nbIC
    IC_FreqVal_contrast2    = squeeze(EEG.spectrumIC_sub_contrast2(:,f_ix,ismember(cfgEEG.subRoles(1,:),cfgANA.sub_contrast2)))';
	
    % -- Perform t-test for OBSERVATION SET
    [h,p] = ttest(log10(IC_FreqVal_contrast1),log10(IC_FreqVal_contrast2));
    RESULTS.tTest_p_OBS(:,f_ix)	= p;
    
    % -- Perform t-test for PERMUTATION SET
    IC_FreqVal_perm1        = zeros(cfgEEG.nbSubTot,cfgEEG.nbIC);
    IC_FreqVal_perm2        = zeros(cfgEEG.nbSubTot,cfgEEG.nbIC);
    for perm_ix = 1:cfgANA.nbPerms
        % draw specific permutation using 'permMat' logical indexes
        IC_FreqVal_perm1(RESULTS.permMat(perm_ix,:),:)  = IC_FreqVal_contrast1(RESULTS.permMat(perm_ix,:),:);
        IC_FreqVal_perm1(~RESULTS.permMat(perm_ix,:),:) = IC_FreqVal_contrast2(~RESULTS.permMat(perm_ix,:),:);
        IC_FreqVal_perm2(~RESULTS.permMat(perm_ix,:),:) = IC_FreqVal_contrast1(~RESULTS.permMat(perm_ix,:),:);
        IC_FreqVal_perm2(RESULTS.permMat(perm_ix,:),:)  = IC_FreqVal_contrast2(RESULTS.permMat(perm_ix,:),:);
        % compute t-test on this permutation
        [h,p] = ttest(log10(IC_FreqVal_perm1),log10(IC_FreqVal_perm2));
        tTest_p_PERM_temp(:,f_ix,perm_ix) = p;
    end
    
    % -- Compute spectrum 95% percentile across subjects (for display only)
    RESULTS.tTest_prctile_contrast1(:,f_ix,1) = prctile(IC_FreqVal_contrast1, 5);
    RESULTS.tTest_prctile_contrast1(:,f_ix,2) = prctile(IC_FreqVal_contrast1, 95);
    RESULTS.tTest_prctile_contrast2(:,f_ix,1) = prctile(IC_FreqVal_contrast2, 5);
    RESULTS.tTest_prctile_contrast2(:,f_ix,2) = prctile(IC_FreqVal_contrast2, 95);
end

% -- keep only p-values that are minimum across all frequencies
RESULTS.tTest_pmin_PERM = squeeze(min(tTest_p_PERM_temp,[],2));
RESULTS.tTest_pmin_PERM = sort(RESULTS.tTest_pmin_PERM,2,'ascend');
clear tTest_p_PERM_temp

% -- T-max (p-min) --> process significance level from permutation values
for f_ix = 1:RESULTS.tTest_nbFreqs 
    for IC_ix = 1:cfgEEG.nbIC
        RESULTS.tTest_sig(IC_ix,f_ix) = (1/cfgANA.nbPerms)*length(find(RESULTS.tTest_p_OBS(IC_ix,f_ix) > RESULTS.tTest_pmin_PERM(IC_ix,:)));
    end    
end

% cfgANA.tTest_alpha = 0.1;  % tempo: just to see the change
RESULTS.tTest_nullRej   = RESULTS.tTest_sig <= cfgANA.tTest_alpha; 
RESULTS.tTest_nb_nullRej= length(find(RESULTS.tTest_nullRej)); % number of significant frequencies
clear f_ix couple_ix sub_ix sub_ix_total IC_FreqVal_contrast1 IC_FreqVal_contrast2 contrast_sub_ix
clear h p ci RESULTS.tTest_nbFreqs perm_ix prog IC_FreqVal_perm1 IC_FreqVal_perm2 IC_ix

% -- display results
disp('t-test computed on components spectra.');
if cfgANA.tTest_doAlphaCorr == 1
    disp(['significance level after correction: alpha=' num2str(cfgANA.tTest_alpha) ]);  
else
    disp(['significance level WITHOUT correction: alpha=' num2str(cfgANA.tTest_alpha) ]);  
end
disp(['found: ' int2str(RESULTS.tTest_nb_nullRej) ' significant freqs']);  
disp(['nb IC: ' int2str(cfgEEG.nbIC) ', nb freqs: ' int2str(length(RESULTS.tTest_freqs_index)) ', nb subjects: ' int2str(cfgEEG.nbSubTot)]); 

% This function displays IC scalp maps and spectra computed in two different conditions. 
% It shows significantly different frequencies, as well as confidence intervals.
dispStatSpectraTopoIC(  EEG.spectrumIC_global_constrast1,EEG.spectrumIC_global_constrast2,...
                        EEG.A_global,cfgEEG.eloc,RESULTS.tTest_freqs_axis,cfgANA.tTest_freqRange,...
                        RESULTS.tTest_nullRej,RESULTS.tTest_prctile_contrast1,RESULTS.tTest_prctile_contrast2);


% --- Display spectra averaged in each contrasted effect
% CompScalpMapSpectra(cfgEEG.Faxis,cfgEEG.nbIC,[2 40],cfgEEG.eloc,EEG.A_global,[EEG.spectrumIC_global_constrast1;EEG.spectrumIC_global_constrast2]);

end












