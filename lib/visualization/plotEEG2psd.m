function h=plotEEG(EEG,varargin)
%h=plotEEG2(EEG)
% 
% Quick visualization of EEG
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
    FontSize=24;
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

T=size(EEG.Epoch,2); %epoch size (temporary)
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
             case 'bartlett'
                winEpoch= hamming(T);winEpoch=repmat(winEpoch',[size(EEG.Epoch,1),1,size(EEG.Epoch,3)]);
                sig=winEpoch.*EEG.Epoch;
end

T=size(sig,2);
sig=abs(fft(sig,2^nextpow2(T*2),2));
sig=mean(sig,3);
sig=log10(sig(:,1:size(sig,2)/4+1));
for indP=1:size(sig,1)
sig(indP,:)=smooth(sig(indP,:),max(round(size(sig,2)/EEG.Fs),3));
sig(indP,:)=sig(indP,:)-min(sig(indP,:));
end
    Scale = max(max(sig,[],2));
    mi=repmat(min(min(sig)),[size(sig,1) 1]);
    ma=repmat(Scale,[size(sig,1) 1]);
    
    T=size(sig,2);
    
t = (0:(T-1))./T*EEG.Fs/4;%frequency range

% calculate shift
shift = cumsum([0; abs(ma(1:end-1))-mi(2:end)]);
shift = repmat(shift,1,T);


%plot 'eeg' data

for indP=1:size(sig,1)

x=t;y=sig(indP,:)-shift(indP,:);
plot(x,y,'linewidth',1.2);
hold on;

bv = min(y);
%bv=-shift(indP,1);
f = y>=bv;
ybv = repmat(bv,1,sum(f));
fill([flipud(x(f)');x(f)'],[ybv,y(f)],'r',...
    'FaceAlpha',0.5,...
    'LineStyle','none');

% h = area([x(f);x(f)]',[ybv;y(f)]',... 
%     'ShowBaseLine','off',...
%     'FaceColor','red',...
%     'FaceAlpha',0.5,...
%     'LineStyle','none');
% delete(h(1));
end

%h=plot(t,sig-shift,'LineWidth',1.2,'Color',Color,'Marker',Marker,'LineStyle',LineStyle,'linewidth',LineWidth);


% edit axes
set(gca,'ytick',sort(mean(sig-shift,2)),'yticklabel',LABELS,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic');
%xlhand = get(gca,'ylabel')
%set(xlhand,'fontsize',20)

grid on
ylim([ min(min(-shift)) max(max(-shift+max(Scale)))])
xlim([min(t),max(t)]);
%vline(0,'g:');
vline(5:5:30);
%drawline