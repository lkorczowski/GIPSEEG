% Use this script to concatenate data from several sessions to get average
% P300 over 

addpath(genpath('D:\Stage\Matlab\eeglab13_4_4b' ))  
addpath(genpath('D:\Stage\Matlab\Pilot_multi'),genpath('D:\Stage\Matlab\eeglab13_4_4b'),genpath('D:\Stage\Matlab\biosig\doc'),genpath('D:\Stage\Matlab\biosig\t200_FileAccess') )
%%
directory=('D:\Stage\Matlab\DATA\GDF\');
slong=[];
Flashlong=[];
Ylong= [];
grouplist= {'G010732_s'; 'G021314_s';  'G031136_s' ;'G042248_s'; 'G051718_s'};
groupi=5;
pp=1;
for sesi=1:4
    file=[directory grouplist{groupi} num2str(sesi) '.gdf'] ;
    [s h] = sload(file); 
    


%
% ACSTP script starts here

% Save (cell array with the following items) :
% (1) AEA
% (2) AEA after weights, latency corrections and CSTP
% (3) CSTP{iter} Bs Bt As At : the temporal and spatials filters for each
% Pz
% (4) Pzind : optimal indice of iter such as Pz is the optimal subspace
%       with respect of the mask
%
% Thus Save{2}{Save{4}) gives P300 with optimal estimation.

% prepare data and load parameters
% 
 %indice of the subject in 'data' (if several)
%%prepare data


EEG = struct;
EEG.Fs = h.SampleRate;
EEG.s = s;
EEG.h = h;

        
% % Use 65th channel to obtain flash-code (instead of header)
% %     StimCodePos = h.EVENT.POS(h.EVENT.TYP>=33024 & h.EVENT.TYP<=33045);



if pp ==1
     EEG.s=EEG.s(:,1:32);
    elseif pp==2
    EEG.s=EEG.s(:,33:64);
end

% Check length of recording
lengthrec= size(s,1)/512/60; %  in min

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
    [EEG.Y, Flash]= trigextr(s(:,65),pp); % trigextr3(s(:,65),1,1);  %(3) Sweeps' indices   % don't change for TA2
     %(4) Class of each sweep % Change for TA2
   
     EEG.Flash=Flash;

%%%%%%%%%%%%%%%%%% PREPROCESSING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decimationfactor=4; %(5) put 1 to do nothing
    NOTCH=1; %(6) put 1 to remove 50Hz
    BANDPASS=1; %(7) put 1 to BANDPASS (see the filter below for the cutoff freq)
    f1=1; %(8) low cutoff freq  (bandpass)
    f2=20; %(9) high cutoff freq  (bandpass)
    N=4; %(10) filter order (bandpass)

%%%%%%%%%%%%%%%%%% CSTP PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Window=1*Fs; %(11) the sweeps window will be 1s
    Delays=1; %(12) +/- nb shifted samples allowed for the jitter correction

%user parameters for the ACSTP to improve converge :
    winElec=[14,21,22,23,24,25,27,28,29]; %(13)* electrodes used for latency calculation and Pz selection
    % exemple: Cz(7),P7(9),P3(10),Pz(11),P4(12),P8(13),O1(14),Oz(15),O2(16)
    winTime=[floor((0.050)*128):ceil((0.550)*128)]; %(14)* time window (in sample) used for latency calculation and Pz selection
    % example: 50ms to 550ms

%   *optional


%%%%%%%%%%%%%%%%%%%%%% PREPROCESSING (no user input needed)%%%%%%%%%%%%%%%%
[Edec Fsd Flashd]=preprocessingEEG(E',Fs,[f1 f2 N decimationfactor],Flash);
EEG.Fs=Fsd;
Window=Window/decimationfactor;
clear E;
EEG.s=Edec;%EEG.s=[EEG.s; Edec]
EEG.Flash=Flashd;%EEG.Flash=[EEG.Flash; Flashd]
%EEG.Y=[EEG.Y; Y]



% Concatenate data

slong=[slong; Edec];
Flashlong=[Flashlong; Flashd];
Ylong= [Ylong ;EEG.Y];

end

%%

EEG.s=slong;
EEG.Flash=Flashlong;
EEG.Y=Ylong;

savename =['G0' num2str(groupi) '_allses_sub' num2str(pp) '.mat']
 save(savename, 'EEG') 
% do again for group2 sub2, check huuuuuuge artefacts!
ACSTPoptions.Epoch_size=Window;
    ACSTPoptions.LatencyCorr_max=Delays;
    ACSTPoptions.Mask_Electrodes=winElec;
    ACSTPoptions.Mask_Time=winTime;
    ACSTPoptions.computeClassLat=1;
    
    ACSTPoptions.SubspaceDim=32:-1:28 %% Added to make it faster
    
    
% ACSTP loop Algorithm

tic 
[Epochs FILTER]=ACSTP(EEG,ACSTPoptions);
toc

