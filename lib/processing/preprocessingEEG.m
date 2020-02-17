function [Edec Fsdec Flashdec]=preprocessingEEG(s,Fs,FILT,Flash)
%
% OPTION 1: separated inputs
% [Edec Fs Flash]=preprocessingEEG(s,Fs,FILT,Flash)
%               s: [nb samples x nb channels] raw EEG recordings
%              Fs: scalar (sample rate in Hz)
%            FILT: Parameters of the filter
%                  FILT(1) lowcut freq
%                  FILT(2) hight cut freq
%                  FILT(3) Order bandpass filter
%                  FILT(4) Decimation factor (put 1 or less for nothing)
%           Flash: [nb samples x1 ] Trigger channel of '0' with '1' at the start
%                   of each sweep. There are [nb epochs] '1'.
%
% OPTION 2: EEG structure
% [EEG]=preprocessingEEG(EEG,FILT)
% EEG is a structure which will be filtered according to FILT with
%              Fs: scalar (sample rate in Hz)
%         Trigger: [nb samples x1 ] Trigger channel of '0' with '1' at the start
%                   of each sweep. There are [nb epochs] '1'.
%      EpochClass: [nb epochs x1] class of the sweeps (i.g. 0 for Non-TARGET, 1
%                   for TARGET).
%          Channels: [nb samples x nb channels] raw EEG recordings
%                  (working also with old EEG structures when called EEG.Channels)
%   NoiseTrigger*: the equivalent of Trigger but for the sweep of the noise.
%                   By default, it takes the same.
% ElectrodesName*: {1 x nb channels} the names of the electrodes (usefull
%                  only in case of plot, i.e. ACSTPoptions.DISPLAY=true)
%              Fs: scalar (sample rate in Hz)
% FILT: Parameters of the filter
%                  FILT(1) lowcut freq
%                  FILT(2) hight cut freq
%                  FILT(3) Order bandpass filter
%                  FILT(4) Decimation factor (put 1 or less for nothing)
%
%         WARNING: BE CAREFULL AT THE POSSIBLE FILTER DISTORSION, please
%                  consider the good pratices for EEG filtering.
%                  "Digital filter design for electrophysiological data
%                  – a practical approach" (Widmann et al. 2014)
%                  Exemple to limit distorsion for ERP: 0.75-40Hz band pass
%                  filter. use FILTFILT matlab function for zero-phase
%                  distorsion.
%
% *** History: 2015-03-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
% see also: script_test_filter_distorsion
Edec=[];

    
if isempty(s)
    error('Empty Input')
end
if isstruct(s)
    Fsdec=[];
    Flashdec=[];
    EEG=s;
    FILT=Fs;
    if length(FILT)~=4
    error('preprocessingEEG: FILT input size incorrect, it should be 4')
end
    if ~isfield(EEG,'Channels') & isfield(EEG,'Signal')
        EEG.Channels=EEG.Signal; %convert for all version of the EEG struct
    end
    [EEG.Channels EEG.Fs EEG.Trigger]=preprocessingEEG(EEG.Channels,EEG.Fs,FILT,EEG.Trigger);
    
    if isfield(EEG,'NoiseTrigger')
        if ~isempty(EEG.NoiseTrigger)
             EEG.NoiseTrigger=DecimateTrigger(EEG.NoiseTrigger,FILT(4),size(EEG.Channels,1));
        end
    end
    
    if isfield(EEG,'Stimulations')
        if ~isempty(EEG.Stimulations)
            EEG.Stimulations=DecimateTrigger(EEG.Stimulations,FILT(4),size(EEG.Channels,1));
        end
    end
    Edec=EEG;
    
else
    
    %% DECIMATE ALL THE TRIGGERS CHANNELS
    if nargin<4
        Flash=zeros(1,size(s,1));
    end
    if length(Flash)~=size(s,1)
        warning('WARNING ERROR IN FLASH DECIMATION, ')
        if length(Flash)==length(find(Flash))
            %we assume that Trigger gives the position of the EpochClass instead of
            %being a trigger channel at the sample rate Fs
            tmp=zeros(size(s,2),1);
            tmp(Flash)=1;
            Flash=tmp;
            clear tmp
            warning('Flash has been assumed to be the locations instead of a trigger channel at the sample rate Fs')
        else
            warning('ACSTP ERROR3: EEG INPUT ERROR.')
            
        end
    end
    
    if length(FILT)~=4
        error('preprocessingEEG: FILT input size incorrect, it should be 4')
    end
    if length(FILT)>2 %TO REMOVE
        
        if FILT(3)>0 %if the filter order is mimum 1 %TO REMOVE
            %disp('Bandpassing...')
            f1=FILT(1);
            f2=FILT(2);
            N=FILT(3);
            w1=f1/(Fs/2); % LOW cutoff frequency
            w2=f2/(Fs/2); % HIGHT cutoff frequency
            %N=4; % FILTER ORDER
            %[b,a]=butter(N,[w1 w2]);
            [b,a]=butter(N,[w1 w2]);
            s = filtfilt(b,a,s);
            %{
                    %if needed to check the frequency response of the
                    digital filter
                    figure
                    [h,w] = freqz(b,a);
                    plot(w/pi*Fs/2,20*log10(abs(h)))
                    ax = gca;
                    ax.YLim = [-100 20];
                    ax.XTick = 0:.5:2;
                    xlabel('Normalized Frequency (\times\pi rad/sample)')
                    ylabel('Magnitude (dB)')
            %}
            %[b,a]=butter(N,[w1],'high');
            %noise = filtfilt(b,a,s);
            %s=s-noise;
            
            %disp('...done')
        end
    else
        Edec=s; %return s without any filtering
    end
    
    
    %notch filter
    if 1 % NOTCH should be avoided for ERP because of distorsions (Widmann and al. 2014)
        %disp('50Hz notch...')
        wo=50/(Fs/2);bw=wo/35;
        [b,a]=iirnotch(wo,bw);
        s = filtfilt(b,a,s);
        %disp('...done')
    end
    %
    if length(FILT)>3 %TO REMOVE
        %disp('Decimation...')
        decimationfactor=FILT(4);
        if decimationfactor>0 & decimationfactor~=1
            for j=1:size(s,2); %for every channel
                Edec(:,j)=decimate(s(:,j),decimationfactor);
            end
        else
            Edec=s;
        end
        Fsdec=Fs/decimationfactor;
        
        %decimation of the Trigger channel
        Flashdec=DecimateTrigger(Flash,decimationfactor,size(Edec,1));
        %disp('...done')
    else
        %do not decimate
        Edec=s;
        Fsdec=Fs;
        Flashdec=Flash;
    end
end


end
