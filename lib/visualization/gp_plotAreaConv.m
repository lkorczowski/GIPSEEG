function [h1,h2,h3]=gp_plotAreaConv(Tab,varargin)
% [h1,h2]=gp_plotAreaConv(Tab,varargin)
% convergence plot with variability (area) and outliers
%  INPUT :
% ------
%           Tab: is matrix [KxT] K is the number of repetition, T the
%           number of conditions (or iteration).
%
%  OPTIONS :
% ------
% the options are called by their name such as fun(..., 'optionName', value)
%              fs: scalar (sample rate in Hz)
%         channel: [scalar/vector] the choosen channel among N to plot. If
%                   several, it will average over the channels.
%           color: [1 x3 double] color in RGB. Default: [0.5 0.5 0.5] (grey)
%           alpha: [scalar] ratio for the solid area.
%                  Default: 0.05 (only 5% of the extremum will be presented
%                  as butterfly plot, the others as solid area)
%  OUTPUT :
% ------
%           h1: median line, h2: variability area
%
%  EXAMPLE :
% ------
% x=[1e-3,5e-3,1e-2,5e-2,1e-1];Tab=randn(1000,5).*.5.*repmat(x,[1000 ,1])+repmat(1-x,[1000 ,1]);
% gp_plotAreaConv(10*log(abs(Tab)),'X',x);xlabel('signal2noise ratio');ylabel('convergence(dB)');set(gca,'XScale','log')
%
% *** History: 09-Apr-2018
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2018
%
% see also : plotEEGvariability, gp_median_quantile, boxplot

% ------------- Check for options -------------
% initialization
T=size(Tab,2);
if T<2; warning('gp_plotAreaConv: not enough conditions, consider using boxplot instead');end
options = struct('color',[0.5 0.5 0.5],'alpha',0.05,'x',[]);

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

% ------------- plot -------------
hold on
if isempty(options.x)
    X=1:T;
else
    X=options.x;
end

% butterfly plot
for i=1:size(Tab,1)
    plot(X,Tab(i,:),'.','color',options.color); 
end

% preparing area
Xfill=[X,fliplr(X)]; % X coordinate for area
Zsorted=squeeze(sort(Tab,1)); % sorting data
nb=size(Zsorted,1);

% area removing p-quantile & (1-p)-quantile
Y=[Zsorted(ceil(nb*(options.alpha/2)),:),fliplr(Zsorted(floor(nb*(1-options.alpha/2)),:))];
h2=fill(Xfill,Y,options.color,'LineStyle','none');
set(h2,'facealpha',.25) 

% plot median
h1=plot(X,mean(Tab),'Color', options.color*0.5,'LineWidth',3);
hold off
