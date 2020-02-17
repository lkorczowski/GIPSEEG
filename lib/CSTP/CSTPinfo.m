function CSTPinfo
% please use script_ACSTP for full implementation
%
% [Bs Bt As At Class eigV]=ACSTP(X,Y,P,W,winElec,winTime,Pzf,Xnoise,Ynoise)
% Xw=applyCSTP(X,Bs,Bt,As,At,Y)
% [P Class]=EnsembleAverage(E,Y,Flash,Window,W)
% [Emean Class]=meanOverlap(E,Flash,Window,TAG,Weights)
% [X isBad] = epoch_p300(Signal,Target,window,offset,Artefacts)
% [W Xhatk]=WeightsEstimation(X,Y,Bs,Bt,As,At)
% [Latency Crit_TraceN]=LatencyEstimation(Xall,X,Y,W,Bs,Bt,As,At,winElec,winTime)
% FlashCorrected=CorrectLatency(Flash,Latencies,RND)
% Conv=ConvergenceLatencies(Lat1,Lat2)
% [Pzind SNR]=best_Pz(Zbarz,Electrodes,Window)
%
% see also : script_ACSTP, ACSTP, applyCSTP, EnsembleAverage, meanOverlap, epoch_p300, WeightsEstimation,
% LatencyEstimation, CorrectLatency, ConvergenceLatencies, best_Pz
%
% *** History: 2015-03-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)

help CSTPinfo