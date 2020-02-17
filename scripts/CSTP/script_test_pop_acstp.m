%% ERROR 1 : incompatibility of format for the event type (string versus interger)
clear all

load('EEG_ex')

eegdata=EEG.Channels';
fs=EEG.Fs;
Trigger=EEG.Trigger;
Classes=EEG.EpochClass
Elecname=EEG.ElectrodesName



EEGsave=pop_importdata('data', eegdata ,'setname','eegtest','dataformat','matlab',...
    'chanlocs','16caps.locs','nbchan',16,'srate',128);


%events=[Classes find(Trigger)/fs]
events=zeros(size(Trigger));
events(find(Trigger))=Classes+1;
EEGsave.data(17,:)=events;

EEG=pop_chanevent(EEGsave,17,'edge','leading');
% pour info, je ne sais pas pourquoi pop_chanevent me supprime mes chanlocs
% les tags sont 
EEG.chanlocs=EEGsave.chanlocs
pop_acstp(EEG)

% ENCOUNTERED ERRORS 

% WHEN ######## quand l'utilisateur a des epochs labels autres que
% 'string' (dans mon cas c'est 1 (NT) ou 2 (TA) en utilisant pop_chanevent
% ERROR ######## 
% Error using cell/strmatch (line 21)
% Requires character array or cell array of strings as inputs.
% 
% Error in pop_acstp (line 123)
%         indEvent = strmatch(res.event1{iEvent}, eventType, 'exact');

%% ERROR 2 : .epoch is needed in the event structure
clear all

load('EEG_ex')

eegdata=EEG.Channels';
fs=EEG.Fs;
Trigger=EEG.Trigger;
Classes=EEG.EpochClass
Elecname=EEG.ElectrodesName



EEGsave=pop_importdata('data', eegdata ,'setname','eegtest','dataformat','matlab',...
    'chanlocs','16caps.locs','nbchan',16,'srate',128);


%events=[Classes find(Trigger)/fs]
events=zeros(size(Trigger));
events(find(Trigger))=Classes+1;
EEGsave.data(17,:)=events;

EEG=pop_chanevent(EEGsave,17,'edge','leading');
% pour info, je ne sais pas pourquoi pop_chanevent me supprime mes chanlocs
% les tags sont 
EEG.chanlocs=EEGsave.chanlocs

events=EEG.event;
for indE=1:length(events)
    if events(indE).type==2,events(indE).type='TA';else,events(indE).type='NT';end
end
EEG.event=events;


pop_acstp(EEG)

% ENCOUNTERED ERRORS 

% Error in pop_acstp (line 126)
%         indEpochs = [ indEpochs EEG.event(indEvent).epoch ];

%% working exemple : 
% clear all

load('EEG_ex')

eegdata=double(EEG.Channels');
fs=EEG.Fs;
Trigger=EEG.Trigger;
Classes=EEG.EpochClass
Elecname=EEG.ElectrodesName



EEGsave=pop_importdata('data', eegdata ,'setname','eegtest','dataformat','matlab',...
    'chanlocs','16caps.locs','nbchan',16,'srate',128);
EEGsave.data=double(EEGsave.data);

%events=[Classes find(Trigger)/fs]
events=zeros(size(Trigger));
events(find(Trigger))=Classes+1;
EEGsave.data(17,:)=events;

EEG=pop_chanevent(EEGsave,17,'edge','leading');
% pour info, je ne sais pas pourquoi pop_chanevent me supprime mes chanlocs
% les tags sont 
EEG.chanlocs=EEGsave.chanlocs

events=EEG.event;
for indE=1:length(events)
    if events(indE).type==2,events(indE).type='TA';else,events(indE).type='NT';end
    events(indE).epoch=events(indE).latency;
end
EEG.event=events;


pop_acstp(EEG)