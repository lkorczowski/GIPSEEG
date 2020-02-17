function h=plotEEG2(EEG,varargin)
%h=plotEEG2(EEG)
% 
% Quick visualization of EEG normalized (unit variance)
%%%% INPUTS :
% ------
% EEG is a structure with
%              Fs: scalar (sample rate in Hz)
%         Trigger: [nb samples x1 ] Trigger channel of '0' with '1' at the start
%                   of each sweep. There are [nb epochs] '1'.
%      EpochClass: [nb epochs x1] class of the sweeps (0 for Non-TARGET, 1
%                   for TARGET).
%        Channels: [nb samples x nb channels] preprocessed continuous EEG recordings
%   NoiseTrigger*: the equivalent of Trigger but for the sweep of the noise.
%                   By default, it takes the same.
% ElectrodesName*: {1 x nb channels} the names of the electrodes (useful
%                  for figures)
%         offset*: integer, number of sample
%          Epoch*: [nb channels x nb samples x nb epochs] epoched signal
%
% *** History: 15-Fev-2018
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2018
%
% ------------- Check for options -------------
% initialization
options = struct('type','averageTA');
   LineWidth=1 ;
   LineStyle='-' ;
   Marker='none' ;
    Color=[.1 .1 .1];
    FontSize=12;
optionNames = fieldnames(options);

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    error('plotEEG2 needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
    inpName = lower(pair{1}); % make case insensitive
    if any(strcmp(inpName,optionNames))
        options.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name',inpName)
    end
end


LABELS=flipud(EEG.ElectrodesName);
switch options.type
    case 'averageTA'
      sig=mean(EEG.Epoch(:,:,EEG.EpochClass==1),3);
       case 'averageNT'
      sig=mean(EEG.Epoch(:,:,EEG.EpochClass==0),3);
    case 'average'
      sig=mean(EEG.Epoch,3);
          case 'all'
             sig=EEG.Epoch(:,:);
                 case 'allTA'
                     tmp=EEG.Epoch(:,:,EEG.EpochClass==1);
             sig=tmp(:,:);clear tmp;
                 case 'allNT'
                               tmp=EEG.Epoch(:,:,EEG.EpochClass==0);
             sig=tmp(:,:);clear tmp;

end
    T=size(sig,2);

VAR=var(sig-repmat(mean(sig,2),[1,T]),[],2);sig=sig./repmat(sqrt(VAR),[1,T]);%unit var
    Scale = quantile(abs(max(sig,[],1)),0.95);
mi=repmat(-Scale,[size(sig,1) 1]);
    ma=repmat(Scale,[size(sig,1) 1]);
    T=size(sig,2);

if ~isfield(EEG,'offset')
EEG.offset=0;
end
    timeoffset=(EEG.offset)/EEG.Fs;
    t=timeoffset:1/EEG.Fs:(timeoffset+T/EEG.Fs)-1/EEG.Fs;


% calculate shift
shift = cumsum([0; abs(ma(1:end-1))+abs(mi(2:end))]);
shift = repmat(shift,1,T);

%plot 'eeg' data
h=plot(t,sig-shift,'LineWidth',1.2,'Color',Color,'Marker',Marker,'LineStyle',LineStyle,'linewidth',LineWidth);


% edit axes
set(gca,'ytick',sort(mean(sig-shift,2)),'yticklabel',LABELS,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic');
%xlhand = get(gca,'ylabel')
%set(xlhand,'fontsize',20)

grid on
ylim([ min(min(-shift-Scale*2)) max(max(-shift+Scale*2))])
xlim([min(t),max(t)]);
startepoch=timeoffset:size(EEG.Epoch,2)/EEG.Fs: max(t);
vline(startepoch,'k');
vline(startepoch-timeoffset,'g:');
if length(startepoch)==length(EEG.EpochClass)
vline(startepoch(EEG.EpochClass==1)-timeoffset+0.3,'r:')
end
%drawline