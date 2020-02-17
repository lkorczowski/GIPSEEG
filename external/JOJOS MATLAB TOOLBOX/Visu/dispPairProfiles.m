function dispPairProfiles(DATA, xaxis, xlab, ylab, smoothEnabled)
% This function enables to browse in a 3D matrix (profiles are in 3rd dim)
%
% INPUTS:
% --------
% - DATA:   set of channels (channels*channels*samples)
% - xaxis:  (optional) specified axis for plotting data
% - xlab:   (optional) specified xlabel for plotting data
% - ylab1:   (optional) specified ylabel for plotting data
% - smoothEnabled: (optional) specifies is smoothing is enabled

% HISTORY:
% --------
% First version: 2013-06-20
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

NbChan1         = size(DATA,1);
NbChan2         = size(DATA,2);
NbSample        = size(DATA,3);
SmoothFactor    = 0;            % default smooth factor
ChanSel1        = NbChan1;
ChanSel2        = NbChan2;
FlagBrowse1     = 1;
FlagBrowse2     = 1;
yLims           = [min(min(min(DATA))) max(max(max(DATA))) ];
h               = gcf;
set(h,'Toolbar','figure');

% *** Assertions on inputs
if nargin < 5  smoothEnabled  = 1;  end
if nargin < 4  ylab     = [];   end
if nargin < 3  xlab     = [];   end
if nargin < 2  xaxis    = [1: NbSample];  end

if NbChan1 == 1
    FlagBrowse1 = 0;
end

if NbChan2 == 1
    FlagBrowse2 = 0;
end

if (NbChan1 == 1) && (NbChan2 == 1)
    warning('[DispPairProfiles] data to browse has just one dimension !');
    return 
end


% *** GUI parametrization
figLoc = get(h, 'Position');

if FlagBrowse1
sliderSel1 = uicontrol(  'style','slider','String','IC#','Min',1,'Max',NbChan1, ...
                        'SliderStep',[1/NbChan1 5/NbChan1],'Value',ChanSel1, ...
                        'units','normalized','Position',[.90 .25 .04 .5], ...
                        'BackgroundColor','w','Callback',@redraw_cb);
end

if FlagBrowse2
sliderSel2 = uicontrol(  'style','slider','String','IC#','Min',1,'Max',NbChan2, ...
                        'SliderStep',[1/NbChan2 5/NbChan2],'Value',ChanSel2, ...
                        'units','normalized','Position',[.95 .25 .04 .5], ...
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
    
    D = squeeze(DATA(ChanSel1,ChanSel2,:));
    
    % Smooth using moving average
    if smoothEnabled
        for i = 1: SmoothFactor
            D = filtfilt([1/4 1/2 1/4],1,D); 
        end
    % text('Units','normalized','Position',[1.05 .02],'string',[ 'Smooth: ' int2str(SmoothFactor)]);
    end
    
    plot(xaxis, D,'b','LineWidth',2)
    xlabel(xlab)
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