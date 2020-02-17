function Pxx=cospectra_welch(data,pnts,nbFreqs,nover)
% INPUTS
%  data : [electrodes x samples]
%  pnts : number of points per epochs
%  NFFT : number of points for the fft (output length will be NFFT/2+1)
% nover : number of overlapping samples (default:0, Barlett's estimator).
%         If negative, there will be gap between windows (non-consecutive samples).
% OUTPUTS
%   Pxx : cospectra matrix
nbchan = size(data,1);
T=size(data,2);
if nargin<4,nover=0;end
if nargin<3,nbFreqs=2^8;end
dt=pnts-nover;
if dt<1;error('number of overlapping samples is bigger than the window length');end
% Barlett's method
Window=repmat(hamming(pnts),1,nbchan);
nbSections=floor(T/(pnts-nover));
% for each epoch compute the fft

% estimate the fft for each epoch
niquistF=nbFreqs/2+1;
allfft=zeros(niquistF,nbchan,nbSections);
indS=0;
while ((indS)*dt+pnts)<=T
    indS=indS+1;
    indEpoch=(indS-1)*dt+1:(indS-1)*dt+pnts;
%     max(max(data(:,indEpoch));
    tmp=fft(Window.*data(:,indEpoch)',nbFreqs);
    allfft(:,:,indS)=tmp(1:niquistF,:); %stack the respective fft for each epoch
end
if indS<1;error('window''length too big or data in the wrong dimension');end

% compute the cospectra the fft (bartlett) for
Pxx=zeros(nbchan,nbchan,niquistF);
for indF=1:niquistF
    % the barlett method average the amplitude spectrum over al the
    % trials
        Pxx(:,:,indF)=abs(squeeze(allfft(indF,:,:))*squeeze(allfft(indF,:,:))')/(indS*(indS+1))*2;
end
