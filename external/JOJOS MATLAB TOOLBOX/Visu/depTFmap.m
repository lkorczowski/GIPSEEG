function depTFmap(DATA, cfg_DEP, cfg_DISP)        
%   function depTFmap(DATA,FS,LABELS,TRIGGER,WIN_LENGTH,COLOR,LINE_WIDTH)    
% 
% This function displays T/F map of (non)-linear uni/multivariate
% dependence measure applied to input data. This measure is processed
% and can be tuned online.
%
% *** Inputs *** 
% - DATA            --> 2D-data [channels x time samples]
% - cfg_DEP         --> data  and dependence configuration structure with:
%       - FS:           (optional)  sampling frequency in Hz
%       - S_DATA:       (optional)  size of each dataset (S must be a column vector)
%                                   if not set each row is considered as a single dataset
%       - LINEARITY :   (optional)  specifies if DFT is normalized (non-linear 
%                                   case, value 0) or not (linear case, value 1)
%       - AVERAGING :   (optional)  vector with 2 positive integers [av_t av_f],
%                                   specifies number of time windows 'av_t' and/or frequency bins 'av_f' 
%                                   on which cospectra are averaged for a better estimation, i.e.
%                                   output frequency resolution is: (av_f * sample_rate) / WINDOW_SIZE
%                                   output time resolution is:      WINDOW_SIZE / (av_t * sample_rate)
%       - DO_REG :      (optional)  when set to 1, regularize cross-spectra (default is 0)
%       - WINDOW_SIZE:  (optional)  size of the window in number of samples
%       - ZERO_PAD:     (optional)  zero-padding factor (adds ZERO_PAD*WINDOW_SIZE zeros to fft)
%       - OVERLAP:      (optional)  percentage of the overlapping window (from 0 to 0.99)
%                                   all samples are considered uniformely only if OVERLAP has value 1 - 2^(-n)
%       - FREQ_RANGE:   (optional)  average output over specified frequency range (only one freq then)
%                                   format must be [sample_rate min_freq_in_Hz max_freq_in_Hz]
%       - WINDOW_TYPE:  (optional)  string with the type of the window
%
% - cfg_DISP         --> display configuration structure with:
%       - LABELS:       (optional)  name of channels in DATA (string array)
%       - TRIGGER:      (optional)  trigger, can be a vector with time samples: 
%                                       [t1,t2,...,tn]
%                                   or a 2D-cell array including trigger name:
%                                       TRIGGER{1:n,1}=(t1,...,tn) 
%                                       TRIGGER{1:n,2}=(name1,...name_n)
%       - WIN_LENGTH:   (optional)  length of window to be displayed (in seconds)
%       - COLOR:        (optional)  '0' for black and white, '1' for color,
%                                   or directly a matrix of [channels , 3] RGB values
%                                   or directly a cell structure of {channels , 3 } RGB values
%       - LINE_WIDTH:   (optional)  vector with line widths to display (channels elements) 
%       - FREQ_RANGE:   (optional)  display output only in range [min_freq_in_Hz max_freq_in_Hz]
%       - TIME_AXIS:    (optional)  time axis (vector with as many time samples as DATA) 
%
%
%  *** History *** 
% First version: 25/03/2013
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


% assess display configuration 
if (nargin < 3)  cfg_DISP = []; end
if (~isfield(cfg_DISP, 'LABELS'))
    cfg_DISP.LABELS = [repmat('data',size(DATA,1),1) int2str((1:size(DATA,1))')];
elseif length(cfg_DISP.LABELS) ~= size(DATA,1)
    error('Not the same number of labels and data channels');
end
if (~isfield(cfg_DISP, 'TRIGGER'))
	display_trigger  = 0;  
else
    display_trigger = 1; 
    if iscell(cfg_DISP.TRIGGER)   
       TrigVal = cell2mat(cfg_DISP.TRIGGER(:,1)); 
    else
       TrigVal = cfg_DISP.TRIGGER;
    end
end   
if (~isfield(cfg_DISP, 'COLOR'))         cfg_DISP.COLOR = 1; end
if (~isfield(cfg_DISP, 'WIN_LENGTH'))    cfg_DISP.WIN_LENGTH  = 10;end

% assess data and dependence configuration  
if (nargin < 2)   cfg_DEP = []; end
if (~isfield(cfg_DEP, 'FS'))            cfg_DEP.FS          = 128; end    
if (~isfield(cfg_DEP, 'S_DATA'))        cfg_DEP.S_DATA      = []; end      
if (~isfield(cfg_DEP, 'LINEARITY')) 	cfg_DEP.LINEARITY   = []; end     
if (~isfield(cfg_DEP, 'AVERAGING')) 	cfg_DEP.AVERAGING   = []; end     
if (~isfield(cfg_DEP, 'WINDOW_SIZE'))   cfg_DEP.WINDOW_SIZE = []; end     
if (~isfield(cfg_DEP, 'ZERO_PAD'))      cfg_DEP.ZERO_PAD    = []; end  
if (~isfield(cfg_DEP, 'DO_REG'))        cfg_DEP.DO_REG      = []; end 
if (~isfield(cfg_DEP, 'OVERLAP'))   	cfg_DEP.OVERLAP     = []; end     
if (~isfield(cfg_DEP, 'FREQ_RANGE'))	cfg_DEP.FREQ_RANGE  = []; end     
if (~isfield(cfg_DEP, 'WINDOW_TYPE'))   cfg_DEP.WINDOW_TYPE = []; end     
if (~isfield(cfg_DISP, 'FREQ_RANGE'))   cfg_DISP.FREQ_RANGE  = [0 cfg_DEP.FS/2]; end  


% various init
global CURFIG;  
global DEP;         % defines DEP as a global variable
Ns  = size(DATA,2);
if (isfield(cfg_DISP, 'TIME_AXIS'))    
    t_axis  =  cfg_DISP.TIME_AXIS; 
else
    t_axis      = [0:Ns-1]./cfg_DEP.FS;
end
total_win   = round(2*t_axis(end)/cfg_DISP.WIN_LENGTH);  % total number of time windows
win_sel     = 1;        % select first time window
isGrid      = 0;        % default grid display
CURFIG= figure('Toolbar','figure','Position', get(0,'Screensize')); % Maximize figure.
pos = get(gca,'Position');
t1 = get(gca,'TightInset');
set(gca,'Position',[2.05*t1(1) pos(2) pos(3)+.65*pos(1)-1.5*t1(1) 1.05*pos(4)]);


sliderSel1  = uicontrol(    'style','slider','Min',1,'Max',total_win+1,'SliderStep',[.05/total_win 3/total_win],...
                            'units','normalized','Position',[.35 .02 .33 .03], 'BackgroundColor','w', ...
                            'Value',win_sel, 'Callback',@redraw_cb);
                    
winSizePlus = uicontrol(     'style','pushbutton', 'String', '+', ...
                            'units','normalized','Position',[.96 .15 .03 .03], ...
                            'BackgroundColor','w','Callback',@incWinSize_cb);

winSizeMinus = uicontrol(    'style','pushbutton', 'String', '-', ...
                            'units','normalized','Position',[.92 .15 .03 .03], ...
                            'BackgroundColor','w','Callback',@decWinSize_cb);
                        
buttonTimeMinus = uicontrol(     'style','pushbutton', 'String', '+', ...
                            'units','normalized','Position',[.78 .02 .03 .03], ...
                            'BackgroundColor','w','Callback',@decTimeScale_cb);

buttonTimePlus = uicontrol(    'style','pushbutton', 'String', '-', ...
                            'units','normalized','Position',[.73 .02 .03 .03], ...
                            'BackgroundColor','w','Callback',@incTimeScale_cb);

buttonGridSwitch = uicontrol(    'style','pushbutton', 'String', 'Grid', ...
                            'units','normalized','Position',[.92 .22 .03 .03], ...
                            'BackgroundColor','w','Callback',@switchGrid_cb);
                        
buttonTrigSwitch = uicontrol(    'style','pushbutton', 'String', 'Trigger', ...
                            'units','normalized','Position',[.96 .22 .03 .03], ...
                            'BackgroundColor','w','Callback',@switchTrig_cb);
                        
ZPplus      = uicontrol(    'style','pushbutton', 'String', '+', ...
                            'units','normalized','Position',[.96 .29 .03 .03], ...
                            'BackgroundColor','w','Callback',@incZP_cb);
                        
ZPminus    	= uicontrol(    'style','pushbutton', 'String', '-', ...
                            'units','normalized','Position',[.92 .29 .03 .03], ...
                            'BackgroundColor','w','Callback',@decZP_cb);

buttonRegSwitch = uicontrol(    'style','pushbutton', 'units','normalized',...
                                'Position',[.92 .50 .03 .03], 'BackgroundColor','w','Callback',@switchReg_cb);
if cfg_DEP.DO_REG        set(buttonRegSwitch,'String', 'REG');
else                     set(buttonRegSwitch,'String', 'no REG'); end
    
buttonLinSwitch = uicontrol(    'style','pushbutton', 'String', 'LIN', 'units','normalized',...
                                'Position',[.96 .50 .03 .03], 'BackgroundColor','w','Callback',@switchLin_cb);
if cfg_DEP.LINEARITY    set(buttonLinSwitch,'String', 'LIN');
else                    set(buttonLinSwitch,'String', 'NON-LIN'); end     

buttonUpdateProc = uicontrol(    'style','pushbutton', 'String', 'Update', ...
                            'units','normalized','Position',[.92 .55 .06 .03], ...
                            'BackgroundColor','w','Callback',@procDep_cb);
                        
TimeAveragePlus = uicontrol(     'style','pushbutton', 'String', '+', ...
                            'units','normalized','Position',[.96 .40 .03 .03], ...
                            'BackgroundColor','w','Callback',@incTav_cb);

TimeAverageMinus = uicontrol(    'style','pushbutton', 'String', '-', ...
                            'units','normalized','Position',[.92 .40 .03 .03], ...
                            'BackgroundColor','w','Callback',@decTav_cb);
                                        
WinSizeText = annotation('textbox', 'Position', [.93 .10 .03 .03], 'FitBoxToText', 'on', ...
                        'String',['Win: ' int2str(cfg_DEP.WINDOW_SIZE)], 'EdgeColor','red');  
                    
ZPText = annotation('textbox', 'Position', [.93 .34 .03 .03], 'FitBoxToText', 'on', ...
                        'String',['ZP: ' int2str(cfg_DEP.ZERO_PAD)], 'EdgeColor','red');  
                    
TavText = annotation('textbox', 'Position', [.93 .45 .03 .03], 'FitBoxToText', 'on', ...
                        'String',['t-AV: ' int2str(cfg_DEP.AVERAGING(1))], 'EdgeColor','red');  

procDep_cb();               
                    
% callback function called when changing slider value
function redraw_cb(h,temp)
    set(0, 'currentfigure', CURFIG);
    
    % get new value from slider 
    win_sel = get(sliderSel1, 'value');
    if win_sel > 0
        t2 	= clamp(round(win_sel*cfg_DISP.WIN_LENGTH*cfg_DEP.FS),[1 Ns]);
        t1  = clamp(t2-cfg_DISP.WIN_LENGTH*cfg_DEP.FS,[1 Ns]);
    else
        t1 = 1; 
        t2 = Ns;
    end
    
    NbFreq = size(DEP,1);
    f_axis = [0:NbFreq-1] * .5*cfg_DEP.FS/NbFreq;

    % plot time / frequency map
    cla,            % clear previous plots
    legend off,
    axis manual,    % prevent automatic scaling   
    surf_in = {t_axis(t1:t2)', f_axis', DEP(:,t1:t2)};     % 3D inputs for surf function
%     surf(surf_in{:},'EdgeColor','none');
    imagesc(surf_in{:}),
    set(gca,'YDir','normal'), % reverse y-axis from bottom to top
    colormap(jet),
	colorbar('location','EastOutside'),  
   	xlabel('Time (s)');
    ylim(cfg_DISP.FREQ_RANGE);
%     xlim([t1 t2]/cfg_DEP.FS);
    caxis([0 1]), % dependance is always comprised in range [0 1]
    box on,
    hold off,
    if isGrid == 1
       grid on 
    end
    
    % draw triggers
    if display_trigger
        trig_inWin = find((TrigVal>t1)&(TrigVal<t2)); % is there triggers in selected window?
        for trig_ix = 1:length(trig_inWin) 
            trig_Xval = TrigVal(trig_inWin(trig_ix))  ; % trigger value in time sample
            line([t_axis(trig_Xval) t_axis(trig_Xval)], get(gca,'YLim'),'color','k','LineWidth', 3, 'LineStyle', '--')
            % plot trigger name
            if iscell(cfg_DISP.TRIGGER) 
                text(t_axis(trig_Xval),max(get(gca,'YLim'))+.1,cfg_DISP.TRIGGER(trig_inWin(trig_ix),2),'FontSize',14,'FontWeight','bold');
            end
        end
    end
end


% callback: process mutlivariate dependence measure
function procDep_cb(h, buttonUpdateProc)
    DEP = procDependence3(DATA, cfg_DEP.S_DATA, cfg_DEP.LINEARITY, ...
                    cfg_DEP.AVERAGING, cfg_DEP.WINDOW_SIZE, cfg_DEP.ZERO_PAD, cfg_DEP.DO_REG, ...
                    cfg_DEP.OVERLAP, cfg_DEP.FREQ_RANGE, cfg_DEP.WINDOW_TYPE,0,1);
    if(DEP~=-1)     % this happens when user cancels computation
        redraw_cb() 
    end
end


% callback: increase window size for fft calculation
function incWinSize_cb(h, winSizePlus)
    cfg_DEP.WINDOW_SIZE = clamp(cfg_DEP.WINDOW_SIZE*2,[2 256]);
    set(WinSizeText,'String',['Win: ' int2str(cfg_DEP.WINDOW_SIZE)]),  
    procDep_cb();
end

% callback: decrease window size for fft calculation
function decWinSize_cb(h, winSizeMinus)
    cfg_DEP.WINDOW_SIZE = clamp(cfg_DEP.WINDOW_SIZE/2,[2 256]);
    set(WinSizeText,'String',['Win: ' int2str(cfg_DEP.WINDOW_SIZE)]),
    procDep_cb();
end

% callback: increase zero padding factor for fft calculation
function incZP_cb(h, ZPplus)
    cfg_DEP.ZERO_PAD = clamp(cfg_DEP.ZERO_PAD*2,[1 128]);
    set(ZPText,'String',['ZP: ' int2str(cfg_DEP.ZERO_PAD)]),  
    procDep_cb();
end

% callback: decrease zero padding factor for fft calculation
function decZP_cb(h, ZPminus)
    cfg_DEP.ZERO_PAD = clamp(cfg_DEP.ZERO_PAD/2,[1 128]);
    set(ZPText,'String',['ZP: ' int2str(cfg_DEP.ZERO_PAD)]),
    procDep_cb();
end


% callback: increase zero padding factor for fft calculation
function incTav_cb(h, TimeAveragePlus)
    cfg_DEP.AVERAGING(1) = clamp(cfg_DEP.AVERAGING(1)*2,[1 64]);
    set(TavText,'String',['t-AV: ' int2str(cfg_DEP.AVERAGING(1))]),  
    procDep_cb();
end

% callback: decrease zero padding factor for fft calculation
function decTav_cb(h, TimeAverageMinus)
    cfg_DEP.AVERAGING(1) = clamp(cfg_DEP.AVERAGING(1)/2,[1 64]);
    set(TavText,'String',['t-AV: ' int2str(cfg_DEP.AVERAGING(1))]),
    procDep_cb();
end


% callback: increase time scale
function incTimeScale_cb(h, buttonTimePlus)
    cfg_DISP.WIN_LENGTH = 2^nextpow2(cfg_DISP.WIN_LENGTH * 2);
    if cfg_DISP.WIN_LENGTH > Ns/cfg_DEP.FS
        cfg_DISP.WIN_LENGTH = Ns/cfg_DEP.FS;
        win_sel = 0;
        set(sliderSel1,'Min',0,'Max',1,'value',win_sel,'visible','off');
    else
        total_win   = ceil(t_axis(end)/cfg_DISP.WIN_LENGTH);  % total number of time windows
        win_sel = clamp(get(sliderSel1, 'value'),[1 total_win]);
        set(sliderSel1,'Max',total_win,'SliderStep',[.05/total_win 3/total_win],'value',win_sel,'visible','on');
    end
    redraw_cb()
end

% callback: decrease time scale
function decTimeScale_cb(h, buttonTimeMinus)
    cfg_DISP.WIN_LENGTH = 2^nextpow2(cfg_DISP.WIN_LENGTH / 2);
    if cfg_DISP.WIN_LENGTH < .5
       cfg_DISP.WIN_LENGTH = .5; 
    end
    total_win   = ceil(t_axis(end)/cfg_DISP.WIN_LENGTH);  % total number of time windows
    win_sel = clamp(get(sliderSel1, 'value'),[1 total_win]);
    set(sliderSel1,'Min',1,'Max',total_win,'SliderStep',[.05/total_win 3/total_win],'value',win_sel,'enable','on','visible','on');
    set(buttonTimePlus,'enable','on','visible','on');
    redraw_cb()
end

% callback: switch regularization
function switchReg_cb(h, buttonRegSwitch)
    cfg_DEP.DO_REG = ~cfg_DEP.DO_REG;
    if cfg_DEP.DO_REG
        set(h,'String', 'REG');
    else
        set(h,'String', 'no REG');
    end
    procDep_cb();
end

% callback: switch linearity
function switchLin_cb(h, buttonLinSwitch)
    cfg_DEP.LINEARITY = ~cfg_DEP.LINEARITY;
 	if cfg_DEP.LINEARITY
        set(h,'String', 'LIN');
    else
        set(h,'String', 'NON-LIN');
    end
    procDep_cb();
end

% callback: switch grid display
function switchGrid_cb(h, buttonGridSwitch)
    grid
    isGrid = ~isGrid;
end

function switchTrig_cb(h, buttonTrigSwitch)
    if ~isempty(cfg_DISP.TRIGGER) 
        display_trigger = ~display_trigger;
        redraw_cb()
    end
end

function [y] = clamp(x, range)
    y = x;
    y(x<range(1)) = range(1);
    y(x>range(2)) = range(2);
end

end
