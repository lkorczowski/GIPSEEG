function EEGFT=EEG_LK2FT(EEG, prestim,poststim,Overlap)

%   EEG is a structure with
%                Fs: scalar (sample rate in Hz)
%           Trigger: [nb samples x1 ] Trigger channel of '0' with '1' at the start
%                     of each sweep. There are [nb epochs] '1'.
%        EpochClass: [nb epochs x1] class of the sweeps (0 for Non-TARGET, 1
%                     for TARGET).
%          Channels: [nb samples x nb channels] preprocessed EEG recordings
%     NoiseTrigger*: the equivalent of Trigger but for the sweep of the noise.
%                     By default, it takes the same.
%   ElectrodesName*: {1 x nb channels} the names of the electrodes (usefull
%                    only in case of plot, i.e. ACSTPoptions.DISPLAY=true)
if nargin<4
    Overlap=0;
end

EEGFT.label=EEG.ElectrodesName;
EEGFT.time=repmat({-prestim:1/EEG.Fs:(poststim-1/EEG.Fs)},1,length(EEG.EpochClass));
EEGFT.fsample=EEG.Fs;
EEGFT.sampleinfo=EEG.EpochClass;
StimPos=find(EEG.Trigger);
EEGFT.trialinfo=[StimPos-prestim*EEG.Fs StimPos+poststim*EEG.Fs-1];

if Overlap
    
    
    for indSection=1:size(EEGFT.trialinfo,1)/6 %for every section (row / column)
        Trigger{indSection}=zeros(size(EEG.Trigger));
        Trigger{indSection}(EEGFT.trialinfo((indSection-1)*6+1:(indSection)*6),1)=1;%find the position of the section's trials
        EEGFT.trial{indSection}=meanOverlap(EEG.Channels',Trigger{indSection},length(EEGFT.time{1}));%average
    end
    
else
    for indT=1:size(EEGFT.trialinfo,1)
        EEGFT.trial{indT}=EEG.Channels(EEGFT.trialinfo(indT,1):EEGFT.trialinfo(indT,2),:)';
    end
end