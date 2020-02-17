function data = load_EEG_data(Directory,user,ElectrodesName)
% INPUTS: BASEDIR and user 
% such as the data you are loading are in
% BASEDIR\user\SessionX\online*.gdf

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
if exist('sload')~=2
    error('Please install sload.m from ''biosig\t200_FileAccess'' Toolbox')    
end
if nargin<3
    FILT=0;
end
USERDIR = fullfile(Directory,user);

listdir = dir(USERDIR);listdir=listdir([listdir.isdir]);
Nsession = length(listdir)-2;
%l = 0;
for k=1:Nsession
    
       
        file = fullfile(USERDIR,listdir(k+2).name,'*.gdf');
        try
            [s h] = sload(file);
        catch e
            continue;
        end
        try
        session = struct;
        session.Fs = h.SampleRate;
        % Flash
        StimCodePos = h.EVENT.POS(h.EVENT.TYP>=33024 & h.EVENT.TYP<=33045);
        session.StimPos=StimCodePos;
        session.StimCode = h.EVENT.TYP(h.EVENT.TYP>=33024 & h.EVENT.TYP<=33045);
        session.Trigger = zeros(size(s,1),1);
        session.Trigger(StimCodePos)=1;

        %Target
        TargetPos =  h.EVENT.POS(h.EVENT.TYP==33285);
        session.Target = zeros(size(s,1),1);
        session.Target(TargetPos)=1;

        
        % Tag for the classes
        session.EpochClass = session.Target(logical(session.Trigger));
    
        session.Channels = s;
        session.h = h;
        session.name=listdir(k+2).name;
        if nargin>2
        session.ElectrodesName=ElectrodesName;
        end
        data.session(k) = session;
        
        catch
            data.session(k).error=['user' num2str(user) 'session' num2str(k) ' has known an unexpected error during flash convertion. Maybe the recording was aborded.'];
            disp(data.session(k).error)
           continue; 
        end
        
    
end
data.user = user;
data.Nsession = length(data.session);