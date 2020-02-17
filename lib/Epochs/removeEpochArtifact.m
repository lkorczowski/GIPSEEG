function [Flashc Yc]=removeEpochArtifact(Xall,Flash,Y,MaxVoltage)
% *** History: 19-Mar-2015
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
if nargin<4
MaxVoltage=250; %in µV
end
        isBad=zeros(size(Y));
        Flashc=zeros(size(Flash));
       for k=1:size(Xall,3)
           Maxk=max(max(Xall(:,:,k)));
           isBad(k)=Maxk>MaxVoltage;%remove all signals that exceed 500µV
           
       end
       posF=find(Flash);
       Flashc(posF(~isBad))=1;
       Yc=Y(~isBad);
       if count(Yc)~=count(Y)
          disp(['WARNING, ' num2str(count(Y)-count(Yc)) 'epochs as been remove because Voltage exceed ' num2str(MaxVoltage) 'µV'])
       end
       %count(Flashc);
