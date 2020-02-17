function ACSTPcspctp(Epochs,ACSTPstruct,ACSTPoptions,PrintOptions)
% ACSTPcspctp(Epochs,ACSTPstruct,ACSTPoptions,PrintOptions)
% compute the simple CSP, CTP and CSTP and print it (no weight, no
% latencies)
%
% Epochs is a structure with
%     X : [nb chan x nb samples x nb trials] unfiltered trials
%     EpochClass : [nb chan x nb samples x nb trials] trials class
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

N=size(ACSTPstruct.EA,1);
T=size(ACSTPstruct.EA,2);

t = (options.offset:(T-1+options.offset))./options.fs;

%% compare CSTP CSP CTP

for IND=1:4
    nbLINE=12
    nbCOLUMN=5
%     Mapping={'Fp1';%1
%         'Fp2';%2
%         'F5';%3
%         'AFz';%4
%         'F6';%5
%         'T7';%6
%         'Cz';%7
%         'T8';%8
%         'P7';%9
%         'P3';%10
%         'Pz';%11
%         'P4';%12
%         'P8';%13
%         'O1';%14
%         'Oz';%15
%         'O2'}%16
    
    % find best Pz for each filter
    clear Xcsp Xcstp Xctp Zcstp Zcsp Zctp Bs Bt As At
    
    %such best Pz for every filter
    Pzs=[15:-1:12];
    for indPz=1:length(Pzs)
        [Bs{indPz}]=cellfun(@(x) x(:,1:Pzs(indPz)),ACSTPstruct.fullBs,'UniformOutput',0);
        [Bt{indPz}]=cellfun(@(x) x(:,1:Pzs(indPz)),ACSTPstruct.fullBt,'UniformOutput',0);
        [As{indPz}]=cellfun(@(x) x(:,1:Pzs(indPz)),ACSTPstruct.fullAs,'UniformOutput',0);
        [At{indPz}]=cellfun(@(x) x(:,1:Pzs(indPz)),ACSTPstruct.fullAt,'UniformOutput',0);
        Zcstp{indPz}=applyCSTP(ACSTPstruct.EA,Bs{indPz},Bt{indPz},As{indPz},At{indPz},ACSTPstruct.Class);
        Zcsp{indPz}=applyCSTP(ACSTPstruct.EA,Bs{indPz},{eye(T) eye(T)},As{indPz},{eye(T) eye(T)},ACSTPstruct.Class);
        Zctp{indPz}=applyCSTP(ACSTPstruct.EA,{eye(N) eye(N)},Bt{indPz},{eye(N) eye(N)},At{indPz},ACSTPstruct.Class);

    end
    
    Pz1=best_Pz(Zcstp,ACSTPoptions.Mask_Electrodes,ACSTPoptions.Mask_Time); %CSTP
    Pz2=best_Pz(Zcsp,ACSTPoptions.Mask_Electrodes,ACSTPoptions.Mask_Time); %CSP
    Pz3=best_Pz(Zctp,ACSTPoptions.Mask_Electrodes,ACSTPoptions.Mask_Time); %CTP
    
%         Xcstp=Epochs.Xhat;%applyCSTP(Epochs.X,Bs{Pz1},Bt{Pz1},As{Pz1},At{Pz1},Epochs.EpochClass);
        Xcsp=applyCSTP(Epochs.X,Bs{Pz2},{eye(T) eye(T)},As{Pz2},{eye(T) eye(T)},Epochs.EpochClass);
        Xctp=applyCSTP(Epochs.X,{eye(N) eye(N)},Bt{Pz3},{eye(N) eye(N)},At{Pz3},Epochs.EpochClass);

    switch IND
        case 1
            PLOTX=Epochs.X(:,:,Epochs.EpochClass==1);
            PLOTXm=ACSTPstruct.EA(:,:,2);
            info='AEA';
            
        case 2
            PLOTX=Epochs.Xhat(:,:,Epochs.EpochClass==1);
            PLOTXm=ACSTPstruct.EAcstp(:,:,2);
            info=['ACSTP(' num2str(ACSTPstruct.BestPz(2)) ')'];
            
        case 3
            PLOTX=Xcsp(:,:,Epochs.EpochClass==1);
            PLOTXm=mean(PLOTX,3);
            info=['ACSP(' num2str(Pzs(Pz2)) ')'];
            
        case 4
            PLOTX=Xctp(:,:,Epochs.EpochClass==1);
            PLOTXm=mean(PLOTX,3);
            info=['ACTP(' num2str(Pzs(Pz3)) ')'];
            
    end
    xticklabels = 0:.5:1;
    
    clims=[0 options.scale*2];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT EEG %%%%%%%%%%%%%%%
    subplot(nbLINE,nbCOLUMN,[1 2 3 6 7 8 11 12 13 16 17 18 21 22 23 26 27 28 31 32 33 36 37 38]);plotEEG(PLOTXm,options.scale,128,options.label,options.fontsize); % Xbar TA
    title('(A)')
    %xlabel('Time (s)')
    xticks = linspace(0, 1, numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', [])
    
    text(-0.00,6*options.scale,[options.title ' ' info], 'FontSize',options.fontsize,'fontname','times new roman','FontAngle','italic')
%     text(-0.05,4.5*options.scale,[info], 'FontSize',options.fontsize,'fontname','times new roman','FontAngle','italic')
    
    % global field power on average
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT GFP mean %%%%%%%%%%%%%%%
    subplot(nbLINE,nbCOLUMN,[21 22 23]+20)
    area(t,global_field_power(PLOTXm'),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
    ylim([0 options.scale*3]);grid on;
    ylabel('(B)','FontSize',options.fontsize,'fontname','times new roman','FontAngle','italic');set(gca,'YtickLabel',[])
    xticks = linspace(0, 1.01, numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', [],'FontSize',options.fontsize,'fontname','times new roman','FontAngle','italic')
    
    %single sweep global field power
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT GFP single trial %%%%%%%%%%%%%%%
    subplot(nbLINE,nbCOLUMN,[26 27 28 31 32 33]+20)
    [GFP1]=global_field_power(permute(PLOTX,[2 1 3]));
    imagesc(GFP1',clims); %mean GFP Xbar TA
    colormap(gray)
    colormap(flipud(colormap))

    %set(gca,'XtickLabel',[]);
    %axis([1 1 1 1]);
    set(gca,'YtickLabel',[])
    ylabel('(C)','FontSize',options.fontsize,'fontname','times new roman','FontAngle','italic');
    xticks = linspace(1, 129, numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', [],'FontSize',options.fontsize,'fontname','times new roman','FontAngle','italic')
    
    % average of single sweep global field power
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT mean GFP %%%%%%%%%%%%%%%
    subplot(nbLINE,nbCOLUMN,[36 37 38]+20)
    area(t,mean(global_field_power(permute(PLOTX,[2 1 3])),2),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
    ylim([0 options.scale.^2]);grid on;
    ylabel('(D)','FontSize',options.fontsize,'fontname','times new roman','FontAngle','italic','fontname','times new roman','FontAngle','italic');set(gca,'YtickLabel',[],'FontSize',options.fontsize)
    xlabel('Time (s)','FontAngle','italic','fontname','times new roman');
    xticks = linspace(0, 1, numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'FontSize',options.fontsize,'fontname','times new roman','FontAngle','italic')
        set(gcf, 'PaperPosition', [0 0 20 30]*options.paperfactor,'units','normalized','outerposition',[0.4 0.1 0.5 .9])
%         spaceplots
    if ~isempty(options.outputfolder)
        print(gcf, [options.outputfolder 'class' num2str(options.indclass) '_' info],'-dtiff','-r450')
    end
end