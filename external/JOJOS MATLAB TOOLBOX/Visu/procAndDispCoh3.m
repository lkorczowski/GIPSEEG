function procAndDispCoh3(Data1, Data2, Faxis, Linearity, WindowSize, Overlap, FreqRange, WindowType)
% function ProcAndDispCoh3(Data1, Data2, Faxis, WindowType, WindowSize, Overlap, FreqRange)

% This function processes, displays dependence measures between two
% datasets. 
% Specific channels are chosen by means of a slider.
% Specific coherence measure displayed (Ordinary, Instantaneous and Lagged) 
% are chosen  by means of check boxes.
%
% Inputs:
% - Data1/2         --> signal matrix (N channels by M time samples)
% - Faxis           --> x axis for plotting coherence
% - Linearity       --> (optionnal) specifies if dependence measure is
%                       linear (1) or non-linear (0)
% - FreqRange       --> (optionnal) range of frequencies of interest, 
%                       should be [sample_rate min_frec max_frec]
% - WindowType      --> (optionnal) string with the type of the window
% - WindowSize      --> (optionnal) size of the window in number of samples 
% - Overlap         --> (optionnal) percentage of the overlapping window (from 0 to 0.99)
%
% History:
% --------
% *** 2011-11-07
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

% default windows type
if(nargin < 8)                          WindowType  = 'hanning';    end
 % default frequency range
if(nargin < 7)                          FreqRange   = [];           end
% default overlapping
if(nargin < 6) || (isempty(Overlap))    Overlap     = .75;          end
% default windows size
if(nargin < 5) || (isempty(WindowSize)) WindowSize  = 128;          end
% default measure is linear
if(nargin < 4) || (isempty(Linearity))  Linearity   = 1;            end


NbChan1     = size(Data1,1);
NbChan2     = size(Data2,1);
ChanSel1    = NbChan1;  
ChanSel2    = NbChan2;  
global SmoothFactor
global rho_ordinary
global rho_instant
global rho_lagged

% default smooth factor
SmoothFactor = 1;

% GUI configuration
h = gcf;
set(h,'Toolbar','figure');
figLoc = get(h, 'Position');

checkBox1 = uicontrol(  'style','checkbox','TooltipString','Ordinary', ...
                        'String','O', 'Value', true,...
                        'units','normalized','Position',[.92 .90 .05 .03], ...
                        'BackgroundColor','w','Callback',@reproc_cb);
                    
checkBox2 = uicontrol(  'style','checkbox','TooltipString','Instantaneous', ...
                        'String',' I', 'Value', false, ...
                        'units','normalized','Position',[.92 .85 .05 .03], ...
                        'BackgroundColor','w','Callback',@reproc_cb);    
                    
checkBox3 = uicontrol(  'style','checkbox','TooltipString','Lagged', ...
                        'String',' L', 'Value', true, ...
                        'units','normalized','Position',[.92 .80 .05 .03], ...
                        'BackgroundColor','w','Callback',@reproc_cb);
                    
checkBox4 = uicontrol(  'style','checkbox','TooltipString','Linearity', ...
                        'String',' Lin', 'Value', true, ...
                        'units','normalized','Position',[.92 .70 .10 .03], ...
                        'BackgroundColor','w','Callback',@reproc_cb);                
                    
                    
sliderSel1 = uicontrol( 'style','slider','Min',1,'Max',NbChan1, ...
                        'SliderStep',[1/NbChan1 5/NbChan1],'Value',ChanSel1, ...
                        'units','normalized','Position',[.93 .17 .02 .5], ...
                        'BackgroundColor','w','Callback',@reproc_cb);
                    
sliderSel2 = uicontrol( 'style','slider','Min',1,'Max',NbChan2, ...
                        'SliderStep',[1/NbChan2 5/NbChan2],'Value',ChanSel2, ...
                        'units','normalized','Position',[.96 .17 .02 .5], ...
                        'BackgroundColor','w','Callback',@reproc_cb);
                    
                    
PushPlus = uicontrol( 'style','pushbutton', 'String', '+', ...
                        'units','normalized','Position',[.96 .13 .03 .03], ...
                        'BackgroundColor','w','Callback',@incSmooth_cb);
                    
PushMinus = uicontrol( 'style','pushbutton', 'String', '-', ...
                        'units','normalized','Position',[.92 .13 .03 .03], ...
                        'BackgroundColor','w','Callback',@decSmooth_cb);

reproc_cb(h, ChanSel1, ChanSel2);


% This callback is called every time a slider is changed
function reproc_cb(h, ChanSel1, ChanSel2)
    ChanSel1 = NbChan1-floor(get(sliderSel1, 'value'))+1 ;
    ChanSel2 = NbChan2-floor(get(sliderSel2, 'value'))+1 ;
    Linearity = get(checkBox4, 'value');
    
    % Dependence measure calculation
    [rho_ordinary, rho_instant, rho_lagged] = procDependence( ...
        [Data1(ChanSel1,:);Data2(ChanSel2,:)], [], Linearity, ...
        WindowSize, Overlap, FreqRange, WindowType);
     
    % Smoothing coherence using moving average
    for i = 1: SmoothFactor
        rho_ordinary = filtfilt([1/4 1/2 1/4],1,rho_ordinary); 
        rho_instant = filtfilt([1/4 1/2 1/4],1,rho_instant); 
        rho_lagged = filtfilt([1/4 1/2 1/4],1,rho_lagged);
    end

    % Display coherence
    disp_cb(h, checkBox1, checkBox2, checkBox3);
    
    % This title is specific to current analysis (ssVEP)
    if(Linearity==1) LinStr='Coherence'; else LinStr='Phase synchrony'; end
    title(  [ LinStr ', Subject1 IC # ' int2str(ChanSel1) ', Subject2 IC # ' int2str(ChanSel2)]);
    text('Units','normalized','Position',[1.05 .01],'string',[ 'Smooth: ' int2str(SmoothFactor)]);
    xlabel('Frequency (Hz)');
end

    
% Coherence display
function disp_cb(h, checkBox1, checkBox2, checkBox3)
    clear gcf
    
    if(get(checkBox1, 'value'))
        plot(Faxis, rho_ordinary, 'color', 'b','LineWidth',3)
        hold on
        [legend_h,object_h,plot_h,text_strings] = legend();
        legend(cat(2,text_strings,{'Ordinary'}))

    end
    if(get(checkBox2, 'value'))
        plot(Faxis, rho_instant, 'color', 'r','LineWidth',3)
        hold on
        [legend_h,object_h,plot_h,text_strings] = legend();
        legend(cat(2,text_strings,{'Instantaneous'}))
    end
    if(get(checkBox3, 'value'))
        plot(Faxis, rho_lagged, 'color', 'g','LineWidth',3)
        [legend_h,object_h,plot_h,text_strings] = legend();
        legend(cat(2,text_strings,{'Lagged'}))
    end
    grid on
    hold off
end

function changeLin_cb(h, checkBox4)
    disp('coucou')
    get(checkBox4, 'value')
%     disp(int2str(get(checkBox4, 'value')))
    if Lin ~= get(checkBox4, 'value')
        Lin = get(checkBox4, 'value');
        reproc_cb(h, ChanSel1, ChanSel2);
    end 
end

function incSmooth_cb(h, PushPlus)
    if SmoothFactor < 30
        SmoothFactor = SmoothFactor+1;
        reproc_cb(h, ChanSel1, ChanSel2);
    end 
end

function decSmooth_cb(h, PushMinus)
    if SmoothFactor > 0
        SmoothFactor = SmoothFactor-1;
        reproc_cb(h, ChanSel1, ChanSel2);
    end 
end
end