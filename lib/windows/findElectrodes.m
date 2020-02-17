function [IndicesElectrodes,LocsElectrodes]=findElectrodes(ElectrodesName, SelectedElectrodes)
% [IndicesElectrodes LocsElectrodes]=findElectrodes(ElectrodesName, SelectedElectrodes)
% From a file of electrodes in a randomized order ElectrodesName, find the
% indices of the wanted electrodes SelectedElectrodes
%
% exemple
% ElectrodesName={'Fp1','Fp2','Cz','O1','O2'}
% SelectedElectrodes= {'Cz','O1'}
if nargout<2
    ComputeLocs=0;
else
    ComputeLocs=1;
end

if exist('readlocs') && ComputeLocs
alllocs=readlocs('Standard-10-20-Cap81.locs');
elseif ~ComputeLocs
else
error('WARNING FINDELECTRODES, function READLOCS not found, please check your eeglab install')
end
% LocsElectrodes=struct;
for ind=1:length(SelectedElectrodes)
    tmp(ind,:)=strcmp(SelectedElectrodes{ind},ElectrodesName);
    if exist('alllocs')
        tmp2=alllocs(strcmp({alllocs.labels},SelectedElectrodes{ind}));
        if ~isempty(tmp2)
        LocsElectrodes(ind)=alllocs(strcmp({alllocs.labels},SelectedElectrodes{ind})); %need eeglab installed

        end
    else
        error('findElectrodes : EEGLAB required, please check if "alllocs" function in Matlab path')
    end
end
IndicesElectrodes=find(sum(tmp,1));

