function [h1 h2 h3]=plotEEGvariability(Zbarz,varargin)


% this function is made for variability analysis.
%  INPUT :
% ------
%           Zbarz: is 3D matrice [NxTxk] where N is the number of channels, T is the
%                  number of samples and k is the number of trials
%
%  OPTIONS :
% ------
% the options are called by their name such as fun(..., 'optionName', value)
%              fs : scalar (sample rate in Hz)
%         channel : [scalar/vector] the choosen channel among N to plot. If
%                   several, it will average over the channels.
%           color : [1 x3 double] color in RGB. Default: [0.5 0.5 0.5] (grey)
%           alpha : [scalar] ratio for the solid area.
%                   Default: 0.05 (only 5% of the extremum will be presented
%                            as butterfly plot, the others as solid area).
%           type  :   1: plot area and butterfly (default, slow)
%                     0: plot area only (faster)
%                    -1: plot butterfly only (slow)
%
% *** History: 16-Apr-2015 (updated 04-Jun-2018)
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO et al. "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis"
%
% see also : gp_median_quantile, plotEEG

% ------------- Check for options -------------
% initialization
options = struct('color',[0.5 0.5 0.5],'alpha',0.05,'channel',1:size(Zbarz,1),'fs',1,'offset',0,'type',1);

% read the acceptable names
optionNames = fieldnames(options);

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    error('plotEEGvariability needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
    inpName = lower(pair{1}); % make case insensitive
    if any(strcmp(inpName,optionNames))
        options.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name',inpName)
    end
end

Zbarz=mean(Zbarz(options.channel,:,:),1);% average over the channels
% ------------- plot -------------
hold on
T=size(Zbarz,2);
time = (options.offset:(options.offset+T-1))./options.fs;

meanZbarz=mean(Zbarz,3);

%plot butterfly (can be slow if many observations)
h1=[];
if abs(options.type)==1
    for i=1:size(Zbarz,3)
        h1=plot(time,Zbarz(:,:,i),'color',options.color);
    end
end

h2=[];h3=[];
% plot area
if options.type>=0
X=[time,fliplr(time)];
Zsorted=squeeze(sort(Zbarz,3));
nb=size(Zsorted,2);
Y=[Zsorted(:,ceil(nb*(options.alpha/2)))',fliplr(Zsorted(:,floor(nb*(1-options.alpha/2)))')];
h3=fill(X,Y,options.color);
end

%plot mean
h2=plot(time,meanZbarz,'Color', options.color*0.5,'LineWidth',3);


hold off