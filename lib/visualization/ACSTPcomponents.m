function ACSTPcomponents(ACSTPstruct,PrintOptions)
% ACSTPcomponents(ACSTPstruct,PrintOptions)
% print the ACSTP components
%
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

%% Plot components
 CompLABELS=Generate_Components_Numbers(1:N)';

[As]=ACSTPstruct.fullAs{2};
As=normPOW(As');
[At]=ACSTPstruct.fullAt{2};
At=normEEG(At');


nbColumns=5;%always 5
nbRows=N/(nbColumns-1);
tempocoord=[1:nbColumns:N];
subplot(nbRows,nbColumns,tempocoord);plotEEG((At(1:N,:)),options.scale,options.fs,CompLABELS);xlabel('Time (s)');
for i=1:size(As,1)
    spatiacoord=1:nbColumns*nbRows;
    spatiacoord=spatiacoord(~ismember(spatiacoord,tempocoord));
    %Count=[2 4 6 8 11 13 15 17 20 22 24 26 29 31 33 35];
    subplot(nbRows,nbColumns,[spatiacoord(i)])
topoplot(abs(As(i,:)),'Brain_Invaders_Erwan_16.locs','plotrad',0.55);title(CompLABELS{i},'FontSize',options.fontsize)
caxis([0 1])
end
hb=colorbar
set(hb,'Units','normalized', 'position', [0.92 0.3 0.02 0.6],'FontSize',options.fontsize);
colormap(gray)
colormap(flipud(colormap))
text(-5.8,5.5,options.title, 'FontSize',options.fontsize*1.5,'fontname','times new roman','FontAngle','italic')

        set(gcf, 'PaperPosition', [0 0 30 25]*options.paperfactor,'units','normalized','outerposition',[0.4 0.1 0.5 .9])


  if ~isempty(options.outputfolder)
        print(gcf, [options.outputfolder],'-dtiff','-r450')
  end
  
