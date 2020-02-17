function plot2ChannelSelection(Data1, Data2, xaxis, xlab, ylab1, ylab2, smoothEnabled)
% This script enables to browse TWO sets of channels, using a simple slider.
% inputs: 
% - Data:   set of channels (channels*samples)
% - xaxis:  (optional) specified axis for plotting data
% - xlab:   (optional) specified xlabel for plotting data
% - ylab1/2:   (optional) specified ylabel for plotting data
% - smoothEnabled: (optional) specifies is smoothing is enabled

%  *** History *** 
% Last version: 2011
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

NbChan1         = size(Data1,1);
NbChan2         = size(Data2,1);
NbSample        = size(Data1,2);
SmoothFactor    = 0;            % default smooth factor
ChanSel1        = NbChan1;
ChanSel2        = NbChan2;
FlagBrowse1     = 1;
FlagBrowse2     = 1;
yLims          = [ min(min(Data1)) max(max(Data1))  min(min(Data2)) max(max(Data2))];
yLims = [min(yLims) max(yLims)];
h               = gcf;
set(h,'Toolbar','figure');

% *** Assertions on inputs
if nargin < 7  smoothEnabled  = 1;  end
if nargin < 6  ylab2  = '';  end
if nargin < 5  ylab1  = '';  end
if nargin < 4  xlab  = [ ];  end
if nargin < 3  xaxis = [1: NbSample];  end

if NbChan1 == 1
    plot(xaxis, Data1(ChanSel1,:))
    xlabel(xlab), ylabel(ylab1), ylim(yLims), grid on
    FlagBrowse1 = 0;
end

if NbChan2 == 1
    plot(xaxis, Data2(ChanSel2,:))
    xlabel(xlab), ylabel(ylab2), ylim(yLims), grid on
    FlagBrowse2 = 0;
end

if (NbChan1 == 1) && (NbChan2 == 1)
    warning('function PlotChannelSelection: data to browse has just one row!');
    return 
end


% *** GUI parametrization
figLoc = get(h, 'Position');

if FlagBrowse1
sliderSel1 = uicontrol(  'style','slider','String','IC#','Min',1,'Max',NbChan1, ...
                        'SliderStep',[1/NbChan1 5/NbChan1],'Value',ChanSel1, ...
                        'units','normalized','Position',[.95 .25 .04 .5], ...
                        'BackgroundColor','w','Callback',@redraw_cb);
end

if FlagBrowse2
sliderSel2 = uicontrol(  'style','slider','String','IC#','Min',1,'Max',NbChan2, ...
                        'SliderStep',[1/NbChan2 5/NbChan2],'Value',ChanSel2, ...
                        'units','normalized','Position',[.90 .25 .04 .5], ...
                        'BackgroundColor','w','Callback',@redraw_cb);
end                    

if smoothEnabled    
    PushPlus = uicontrol( 'style','pushbutton', 'String', '+', ...
                            'units','normalized','Position',[.96 .15 .03 .03], ...
                            'BackgroundColor','w','Callback',@incSmooth_cb);

    PushMinus = uicontrol( 'style','pushbutton', 'String', '-', ...
                            'units','normalized','Position',[.92 .15 .03 .03], ...
                        'BackgroundColor','w','Callback',@decSmooth_cb);
end

redraw_cb(h, ChanSel1, ChanSel2);


function redraw_cb(h, ChanSel1, ChanSel2)
    if(FlagBrowse1) ChanSel1 = NbChan1-floor(get(sliderSel1, 'value'))+1 ;
    else            ChanSel1=1;     end
    if(FlagBrowse2) ChanSel2 = NbChan2-floor(get(sliderSel2, 'value'))+1 ;
    else            ChanSel2=1;     end
    
    D1 = Data1(ChanSel1,:);
    D2 = Data2(ChanSel2,:);
    
    % Smooth using moving average
    if smoothEnabled
        for i = 1: SmoothFactor
            D1 = filtfilt([1/4 1/2 1/4],1,D1); 
            D2 = filtfilt([1/4 1/2 1/4],1,D2); 
        end
    % text('Units','normalized','Position',[1.05 .02],'string',[ 'Smooth: ' int2str(SmoothFactor)]);
    end
    
%     plotyy(xaxis, Data1(ChanSel1,:),xaxis, Data2(ChanSel2,:))
    plot(xaxis, D1,'b','LineWidth',2)
    hold on
    plot(xaxis, D2,'r','LineWidth',2)
    xlabel(xlab)
    legend(ylab1,ylab2)
    ylim(yLims);
    title(['Channel # ' int2str(ChanSel1) ' ,Channel # ' int2str(ChanSel2)]) 
    grid on
    hold off
end


function incSmooth_cb(h, PushPlus)
    if SmoothFactor < 20
        SmoothFactor = SmoothFactor+1;
        redraw_cb(h, ChanSel1, ChanSel2)
    end 
end

function decSmooth_cb(h, PushMinus)
    if SmoothFactor > 0
        SmoothFactor = SmoothFactor-1;
        redraw_cb(h, ChanSel1, ChanSel2)
    end 
end
    

end