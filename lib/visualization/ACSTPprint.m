function ACSTPprint(ACSTPstruct,PrintOptions)
% outputfig=['D:\Mes Documents GIPSA\MATLAB\figures\CSPT\ERWAN\EAvsCSTP\Subject_' num2str(Subjects) '.tiff'];
% PrintOptions is a structure with
%     'outputfolder' : output directory and file name without extension (default none)
%     'scale' : scaling for the figures (in µV for EEG) (default 4)
%     'fontsize' : size of the font in the output figures (default 20)
%     'fs' : sample rate (Hz) (default 128 Hz)
%     'offset' : offset in samples (ex: -128 for 1s prestim at 128Hz),
%     (default 0)
%     'label' : the channels names (default none)
%     'paperfactor' : a factor to increase or decrease the printing size
%     (default 1)
%     'indclass' : the indice of the class to print (default : all)

% ------------- Check for options -------------
% initialization (you can copy-past this line to setup your own parameters)
options = struct('outputfolder','','scale',4,'fontsize',20,'fs',128,'offset',0,'label',[],'paperfactor',1,'indclass',size(ACSTPstruct.EA,3),'title','');

% read the acceptable names
optionNames = fieldnames(options);

if nargin<2
    PrintOption=[];
end
% count arguments
if ~isempty(PrintOptions)
    nameArgs=fieldnames(PrintOptions);
    nArgs = length(nameArgs);
    
    for indO = 1: nArgs % pair is {propName;propValue}
        inpName = lower(nameArgs{indO}); % make case insensitive
        if any(strcmp(inpName,optionNames))
            options.(inpName) = PrintOptions.(inpName);
        else
            error('%s is not a recognized parameter name',inpName)
        end
    end
end

t = (options.offset:(size(ACSTPstruct.EA,2)-1+options.offset))./options.fs;

%% AEA versus ACSTP
indClass=options.indclass;
for IND=1:2
    nbLINE=9;
    nbCOLUMN=5;
    
    
    switch IND
        case 1
            %PLOTX=mean(X(:,:,Y==1),3);
            PLOTX=ACSTPstruct.EA(:,:,indClass);
            info=[options.title ' AEA'];
            % options.scale=FindScaling(PLOTX,ACSTPoptions.Mask_Electrodes,ACSTPoptions.Mask_Time);
        case 2
            %PLOTX=Xcstp(:,:,Y==1);
            PLOTX=ACSTPstruct.EAcstp(:,:,indClass);
            info=[options.title ' ACSTP(' num2str(ACSTPstruct.BestPz(indClass)) ')'];
            
    end
    xticklabels = 0:.5:1;
    
    set(gcf, 'PaperPosition', [0 0 12.5 30])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT EEG %%%%%%%%%%%%%%%
    subplot(nbLINE,nbCOLUMN,[1 2 6 7 11 12 16 17 21 22 26 27 31 32 36 37]+2*(IND-1));plotEEG(PLOTX,options.scale,options.fs,options.label); % Xbar TA
    if IND==1
        set(get(gca,'YLabel'),'Rotation',0)
    else
        set(gca,'YtickLabel',[])
    end
    
    if IND>1
        set(gca, 'color', [0.95 .95 .95])
    end
    set(gcf, 'color', [1 1 1])
    
    title(info,'FontSize',options.fontsize)
    %xlabel('Time (s)')
    xticks = linspace(0, 1.01, numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',options.fontsize,'fontname','times new roman','FontAngle','italic')
    set(gcf, 'InvertHardCopy', 'off');
    % global field power on average
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT GFP mean %%%%%%%%%%%%%%%
    subplot(nbLINE,nbCOLUMN,[21 22]+20+2*(IND-1))
    area(t,global_field_power(PLOTX'),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
    ylim([0 6]);grid on;
    if IND==1
        ylabel('GFP','FontSize',options.fontsize,'fontname','times new roman','FontAngle','italic');
        %             set(get(gca,'YLabel'),'Rotation',0)
    else
        set(gca,'YtickLabel',[])
        
    end
    xticks = linspace(0, 1, numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',options.fontsize,'fontname','times new roman','FontAngle','italic')
    set(gca, 'color', [0.95 .95 .95])
    
    set(gcf, 'PaperPosition', [0 0 20 30]*options.paperfactor,'units','normalized','outerposition',[0.4 0.1 0.5 .9])
    

    
end
    if ~isempty(options.outputfolder)
        
        print(gcf, [options.outputfolder '_class' num2str(options.indclass)],'-dtiff','-r450')
    end
%% included function
%%%%%%%%%%%%%%%%%%%%%%%%%%% NESTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function AUsers=Generate_Components_Numbers(Numbers)
        for i=Numbers
            if i<10
                AUsers{i}=['c0' num2str(i)];
            else
                AUsers{i}=['c' num2str(i)];
            end
        end
        AUsers=AUsers(find(~cellfun(@isempty,AUsers)));
    end

    function Xnorm=normPOW(X)
        Xnorm=abs(X);
        Xnorm=Xnorm-repmat((min(Xnorm,[],2)), [1,size(Xnorm,2)]);
        Xnorm=Xnorm./repmat((max(Xnorm,[],2)), [1,size(Xnorm,2)]);
    end


    function Xnorm=normEEG(X,Method,Param)
        if nargin<2
            Method='';
        end
        if strcmp(Method,'fro')
            
            for k=1:size(X,3)
                coef(k)=norm(X(:,:,k),'fro')/length(X(:,:,k));
                Xnorm(:,:,k)=X(:,:,k)/coef(k);
            end
            
        elseif strcmp(Method,'baseline')
            if nargin<2 || isempty(Param)
                winBaseline=96:128; % HARDCODED NOT GOOD !!!!!!!!!!!!! ! ! !
            else
                winBaseline=Param;
            end
            for k=1:size(X,3)
                coef(:,:,k)=repmat(mean(X(:,winBaseline,k),2),1,size(X,2)); %compute the baseline
            end
            Xnorm=X-coef;
            
        else
            Xnorm=X./repmat(sqrt(var(X,[],2)), [1,size(X,2)]);
            
            for k=1:size(Xnorm,3)
                Xnorm(:,:,k)=Xnorm(:,:,k)-repmat(mean(Xnorm(:,:,k),2),[1 size(Xnorm,2)]);
            end
            
        end
    end
    function [gfp]=global_field_power(P,param)
        %%
        % *** History: 19-Mar-2015
        % *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
        % *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
        %
        %close all
        for k=1:size(P,3);
            u=P(:,:,k)';
            gfp(:,k)=sqrt(mean(u.^2,1));
        end
    end
end
