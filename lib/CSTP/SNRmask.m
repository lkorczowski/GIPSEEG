function SNR=SNRmask(Zbarz,winElec,winTime)
%   
%  *** History: 2015-03-19
%  *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
%  *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)

SNR=norm(Zbarz(winElec,winTime),'fro').^2/norm(Zbarz,'fro').^2;