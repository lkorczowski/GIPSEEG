function Xw=applyCSTP(X,Bs,Bt,As,At,Y)
%  Xw=applyCSTP(X,Bs,Bt,As,At,Y)
%   Applies to X an epoch-EEG the CSTP filter such as:
%   Xw(:,:,k)=As*Bs'X(:,:,k)*Bt*At'
%   Y are the class of each epoch (optional).
%
% *** History: 2015-03-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
%   see also: CSTP, ACSTP, WeightsEstimation
if nargin<6
    Y=ones(1,size(X,3));
end
Class=unique(Y);
for k=1:size(X,3)
    c=find(Y(k)==Class);
    Xw(:,:,k)=As{c}*Bs{c}'*X(:,:,k)*Bt{c}*At{c}'; 
end