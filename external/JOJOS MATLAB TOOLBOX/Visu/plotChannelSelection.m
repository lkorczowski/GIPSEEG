function plotChannelSelection(Data, xaxis, xlab, ylab)
% This script enables to browse a set of channels, using a simple slider.
% *** Inputs ***
% - Data:   set of channels (channels*samples)
% - xaxis:  (optional) specified axis for plotting data
% - xlab:   (optional) specified xlabel for plotting data
% - ylab:   (optional) specified ylabel for plotting data
%
%  *** History *** 
% Last version: 2011
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

NbChan      = size(Data,1);
NbSample    = size(Data,2);
ChanSel     = NbChan;
yLims       = [ min(min(Data)) max(max(Data))];
h           = gcf;
set(h,'Toolbar','figure')


if nargin < 4  ylab  = [];  end
if nargin < 3  xlab  = [];  end
if nargin < 2  xaxis = [1: NbSample];  end
if NbChan == 1
    warning('function PlotChannelSelection: data to browse has just one row!');
    plot(xaxis, Data(ChanSel,:))
    xlabel(xlab), ylabel(ylab), ylim(yLims), grid on
    return 
end


figLoc = get(h, 'Position');
sliderSel = uicontrol(  'style','slider','String','IC#','Min',1,'Max',NbChan, ...
                        'SliderStep',[1/NbChan 5/NbChan],'Value',ChanSel, ...
                        'units','normalized','Position',[.94 .25 .05 .5], ...
                        'BackgroundColor','w','Callback',@redraw_cb);
redraw_cb(h, ChanSel);


function redraw_cb(h, ChanSel)
    ChanSel = NbChan-floor(get(sliderSel, 'value'))+1 ;
    plot(xaxis, Data(ChanSel,:))
    xlabel(xlab)
    ylabel(ylab)
    ylim(yLims);
    title(['Channel # ' int2str(ChanSel)]) 
    grid on
end


end