function Xw=applyWeights(X,W)
% Xw=applyWeights(X,W)
% simple loop Xw(:,:,k)=X(:,:,k)*W(k) when X are epoch-EEG and W the
% respective weights
%
% *** History: 2015-03-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
%  see also: WeightsEstimation

for k=1:length(W)
   Xw(:,:,k)=X(:,:,k)*W(k); 
end