function plotSpectraMosaic(data, xaxis, xlab, ylab, pattern, tlab)
% This function enables to watch simultaneously a number of spectra (the
% more, the tinier...)

% inputs: 
% - data:   set of channels (channels*samples)
% - xaxis:  (optional) specifies axis for plotting data
% - xlab:   (optional) specified xlabel for plotting data
% - ylab:   (optional) specified ylabel for plotting data
% - pattern: (optional) specified nb of subplot lines and columns
%                       should be in the format: [nbLines nbCol]
% - titles: (optional) specified titles for each spectrum
%                       should be a structure with same num of titles as data
%                      
% history:
% Last version: 2011
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com



% assertions on inputs
if nargin < 2  xaxis = [1: size(data,2)];  end
if nargin < 3  xlab  = [];  end
if nargin < 4  ylab  = [];  end
if nargin < 6  tlab  = 0;  end

% calculate number of subplot lines and columns
NbSpectra = size(data,1);
if (nargin < 5) || (prod(pattern) < NbSpectra)
    NbLines = round(sqrt(NbSpectra));
    NbCols = NbLines;
    while((NbLines*NbCols) < NbSpectra) NbLines=NbLines+1; end
else
    NbLines = pattern(1);
    NbCols  = pattern(2);
end


% plot spectra
for i=1:NbSpectra
    subplot(NbLines,NbCols,i);
    plot(xaxis,data(i,:),'blue');
    xlabel(xlab)
    ylabel(ylab)
    if iscellstr(tlab)
        title(tlab{i}); 
    end
    hold off
    axis tight;
end
end