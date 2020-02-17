function [P Class]=EnsembleAverage(E,Y,Flash,Window,W)
% [P Class]=EnsembleAverage(E,Y,Flash,Window,W)
%
% Compute the average sweep for each class
% E EEG signals
% Y class for each sweep
% Flash indice of each sweep
% Window size (in sample) of each sweep
%
% P=EnsembleAverage(E) with E a 3D EEG signal with each sweep on the
% third dimension. P is simply the average of each sweep : mean(E,3)
%
% [P Class]=EnsembleAverage(E,Y) with E a 3D EEG signal with each sweep on the
% third dimension with respect of his class in Y. P will be a 3D matrix
% with the averaged sweep of each class in the third dimension and Class,
% the vector of class tagged in the order of P.
%
% [P Class]=EnsembleAverage(E,Y,Flash,Window) with E a full EEG signal. Y is the class each sweep,
% Flash is the indice of each sweep and Window is the length (in sample) of each
% sweep.
% P will be a 3D matrix with the averaged sweep of each class in the third dimension and Class,
% the vector of class tagged in the order of P.
%
% *** History: 2015-03-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
% see also : meanOverlap

if nargin>2
    if isempty(Y)
        Y=ones(length(find(Flash)),1);
    end
    if isempty(Flash) || isempty(Window) || length(find(Flash))~=length(Y) || nargin<4
        error('Invalid Inputs in EnsembleAverage please check Flash or Window')
    end
    X=epoch_p300(E,Flash,Window);
    
end

if nargin<3
        Window=size(E,2);
       X=E; 
       if nargin<2
          Y=ones(1,size(X,3)); 
       end
end

if nargin>4
    X=applyWeights(X,W);
end




Class=unique(Y);
P=zeros(size(E,1),Window,length(Class));
for z=1:length(Class)
   c=find(Y==Class(z));
   P(:,:,z)=mean(X(:,:,c),3);
end

end