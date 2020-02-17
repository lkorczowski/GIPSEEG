function [S,T,w,b] = train_bci(EEG,Fs,mrk,wnd,f,nof,n)
% [S,T,w,b] = train_bci(Raw-Signal, Sample-Rate, Markers, 
%     Epoch-Wnd, Spectral-Flt, Flt-Number, Flt-Length)
% In:  
%       Raw-Signal   :   raw EEG/MEG calibration session data [Samples x Channels]
%       Sample-Rate  :   sampling frequency of the signal in Hz, e.g. 100
%       Markers      :   marker channel: entries in {1,2} indicate different mental conditions,
%                        others are ignored, e.g. sparse(1,[1132 6462 8234 2342],[1 2 2 1])
%       Epoch-Wnd    :   window, in seconds, relative to each condition marker to 
%                        indicate the duration of the conditions, e.g. [0.5 3.5]
%       Spectral-Flt :   spectral filter to be used; function of frequency in Hz, e.g. @(f) f>7 & f<30
%       Flt-Number   :   number of spatial filters to compute per condition via CSP, e.g. 3
%       Flt-Length   :   length, in samples, of the spectral filter, e.g. 200
% Out:
%       S,T,w,b     :   parameters for test_bci
% Example:
%       a = <load via EEGLAB, with markers for condition A as '1' and markers for condition B as '2'>
%       [S,T,w,b] = train_bci(a.data, a.srate, sparse(1,[a.event.latency], strcmp({a.event.type},'1') + 2*strcmp({a.event.type},'2'), ...
%                             [0.5 3.5], @(f) f>7 & f<30, 3, 200);
 
% frequency filtering and temporal filter estimation
[t,c] = size(EEG); idx = reshape(1:t*c-mod(t*c,n),n,[]);
FLT = real(ifft(fft(EEG).*repmat(f(Fs*(0:t-1)/t)',1,c)));
T = FLT(idx)/EEG(idx);
 
% data epoching, class-grouping and CSP
wnd = round(Fs*wnd(1)):round(Fs*wnd(2));
for k = 1:2
    EPO{k} = FLT(repmat(find(mrk==k),length(wnd),1) + repmat(wnd',1,nnz(mrk==k)),:);
end
[V,D] = eig(cov(EPO{2}),cov(EPO{1})+cov(EPO{2}));
S = V(:,[1:nof end-nof+1:end]);
 
% log-variance feature extraction and LDA
for k = 1:2
    X{k} = squeeze(log(var(reshape(EPO{k}*S, length(wnd),[],2*nof))));
end
w = ((mean(X{2})-mean(X{1}))/(cov(X{1})+cov(X{2})))';
b = (mean(X{1})+mean(X{2}))*w/2;
