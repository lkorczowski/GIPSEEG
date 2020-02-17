function [GradFO Xcoor R C Pval]=gradient_frontooccipital(X,LABELS)

load('BImapping32_electrodes_rows')
Rows=[LABELSrow{:,2}];
NbR=max(Rows);
Xcoor=[];
GradFO=[];
ind=1;
for i=1:NbR
    SelectedElectrodes=LABELSrow(Rows==i,1);
    [IndicesElectrodes,LocsElectrodes]=findElectrodes(LABELS, SelectedElectrodes);
    if ~isempty(IndicesElectrodes)
        Xcoor(ind)=median([LocsElectrodes.X]);
        GradFO(ind)=norm(X(IndicesElectrodes,:),'fro');
        ind=ind+1;
    end
end
%GradFO=squeeze(GradFO);

[R,P]  =corrcoef(log(GradFO),Xcoor);
R=R(2);
Pval=P(2);

%C=(Xcoor-mean(Xcoor))'\(log(GradFO)-mean(log(GradFO)))';
 C=(Xcoor(Xcoor>=0)-mean(Xcoor(Xcoor>=0)))'\((GradFO(Xcoor>=0))-mean((GradFO(Xcoor>=0))))';

%figure;plot((Xcoor-mean(Xcoor)),log(GradFO)-mean(log(GradFO)))