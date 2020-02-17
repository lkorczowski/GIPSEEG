function Pxx=cospectra_bartlett(data,pnts,nbFreqs )
% INPUTS
%  data : [electrodes x samples]
%  pnts : number of points per epochs
%  NFFT : number of points for the fft (output length will be NFFT/2+1)
% OUTPUTS
%   Pxx : cospectra matrix
nbchan = size(data,1);
T=size(data,2);

if nargin<3,nbFreqs=2^8;end

% Barlett's method
Window=repmat(hamming(pnts),1,nbchan);
nbSections=floor(T/pnts);
% for each epoch compute the fft

% estimate the fft for each epoch
niquistF=nbFreqs/2+1;
allfft=zeros(niquistF,nbchan,nbSections);
for indS=1:nbSections;
    tmp=zeros(nbFreqs,nbchan);
    indEpoch=(indS-1)*pnts+1:(indS)*pnts;
    tmp=fft(Window.*data(:,indEpoch)',nbFreqs);
    allfft(:,:,indS)=tmp(1:niquistF,:); %stack the respective fft for each epoch
end


% compute the cospectra the fft (bartlett) for
Pxx=zeros(nbchan,nbchan,niquistF);
for indF=1:niquistF
    % the barlett method average the amplitude spectrum over al the
    % trials
        Pxx(:,:,indF)=abs(squeeze(allfft(indF,:,:))*squeeze(allfft(indF,:,:))')/(nbSections*(nbSections+1))*2;
end
