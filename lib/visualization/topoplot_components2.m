function topoplot_components2(TemporalComponents, SpatialComponents,Class,Fs,nb_rows,nb_columns,localization,savePath,Artefacts)
%topoplot_components2(TemporalComponents, SpatialComponents,Class,Fs,nb_rows,nb_columns,localization,savePath,Artefacts)
%TemporalComponents, TemporalComponents, cells [nbClasses x 1] are the
% ACSTP filter (Bt and Bs respectively). The components are in column.
% localization is either the localization structure from eeglab (see
% readlocs.m) or the electrodes names (cells).
%
% *** History: 2015-04-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)

if nargin<8
    savePath=[]; %do not save the figure
end

if nargin<7
    error('Please consider to give the localization of electrodes or at least the electrodes names')
elseif iscell(localization) %localization is a cell of electrodes names i.g. {'Fz','Cz','Pz'} 
   [~, tmpo]=findElectrodes(localization, localization);
   clear localization
   localization=tmpo;
end

%% TOPOPLOT FOR THE COMPONENTS
CompLABELS=Generate_Components_Numbers(1:size(SpatialComponents{1},2))';
for i=1:length(CompLABELS);
    if(floor(Artefacts(i))==1);CompLABELS{i}=[CompLABELS{i} '+'];end;
    if(mod(Artefacts(i),1)==0.1);CompLABELS{i}=[CompLABELS{i} '~'];end;
end
for classIND=1:length(Class)
    FontSize=24;
    %[Bs]=FILTER(indSubject).fullBs{FILTER(2).BestPz(classIND)}{classIND};
    [Bs]=(SpatialComponents{classIND});
    Bs=normEEG(Bs','space');
    %[Bt]=FILTER(indSubject).fullBt{FILTER(2).BestPz(classIND)}{classIND};
    [Bt]=(TemporalComponents{classIND});
    Bt=(normEEG(Bt(:,1:size(Bs,1))'));
    
    TemporalCells=(1:nb_columns:nb_columns*nb_rows);
    TemporalCells=[TemporalCells TemporalCells+1];
    subplot(nb_rows,nb_columns,TemporalCells);plotEEG(Bt,4,Fs,CompLABELS);xlabel('Time (s)')
    vline([0.250 1.5 2.75],'g'); %0ms
    vline([0.550 1.80 3.05],'r'); %300ms
    vline([1.25 2.5],'k'); % -250ms
    title('NT, TA, pinv(D)')
    %text(0.35,1.05,['ss ' num2str(indSubject)], 'FontSize',FontSize*1.5,'fontname','times new roman','FontAngle','italic','units','normalized')
    for i=1:size(Bs,1)
        SpatialCells=(1:nb_columns*nb_rows);
        SpatialCells=SpatialCells(~ismember(SpatialCells,TemporalCells));
        %Count=[2 4 6 8 11 13 15 17 20 22 24 26 29 31 33 35];
        subplot(nb_rows,nb_columns,[SpatialCells(i)]);
        topoplot((Bs(i,:)),localization,'plotrad',0.55,'colormap',colormap(jet(20)));title(CompLABELS{i},'FontSize',FontSize);
        set(gcf, 'PaperPosition', [0 0 30 20]);
        %caxis([-1 1])
        %caxis([quantile(Bs(:),0.1) quantile(Bs(:),0.9)]);
    end
    hb=colorbar;
    set(hb, 'position', [0.92 0.3 0.02 0.6],'FontSize',FontSize);
    %colormap(flipud(colormap));
    
    

    
    %spaceplots
    if nargin>6
        if ~isempty(savePath)
            if ~exist(fileparts(savePath))
                mkdir(fileparts(savePath));
            end
            print(gcf,savePath,'-dpng','-r450')
        end
    end
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NESTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function AUsers=Generate_Components_Numbers(Numbers)
        for nb=Numbers
            if nb<10
                AUsers{nb}=['c0' num2str(nb)];
            else
                AUsers{nb}=['c' num2str(nb)];
            end
        end
        AUsers=AUsers(find(~cellfun(@isempty,AUsers)));
    end

    function Xnorm=normPOW(X)
        Xnorm=abs(X);
        Xnorm=Xnorm-repmat((min(Xnorm,[],2)), [1,size(Xnorm,2)]);
        Xnorm=Xnorm./repmat((max(Xnorm,[],2)), [1,size(Xnorm,2)]);
    end

    function [IndicesElectrodes LocsElectrodes]=findElectrodes(ElectrodesName, SelectedElectrodes)
        % From a file of electrodes in a randomized order ElectrodesName, find the
        % indices of the wanted electrodes SelectedElectrodes
        %
        % exemple
        % ElectrodesName={'Fp1','Fp2','Cz','O1','O2'}
        % SelectedElectrodes= {'Cz','O1'}
        if exist('readlocs')
            alllocs=readlocs('Standard-10-20-Cap81.locs');
        else
            error('ERROR FINDELECTRODES, READLOCS.m not found, please check your eeglab install')
        end
        %LocsElectrodes=struct;
        for ind=1:length(SelectedElectrodes)
            tmp(ind,:)=strcmpi(SelectedElectrodes{ind},ElectrodesName);
            if exist('alllocs')
                tmp2=alllocs(strcmpi({alllocs.labels},SelectedElectrodes{ind}));
                if ~isempty(tmp2)
                    LocsElectrodes(ind)=alllocs(strcmpi({alllocs.labels},SelectedElectrodes{ind})); %need eeglab installed
                    
                end
            end
        end
        IndicesElectrodes=find(sum(tmp,1));
    end
end