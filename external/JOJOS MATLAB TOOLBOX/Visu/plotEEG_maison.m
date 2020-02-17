function plotEEG_maison(DATA,FS,LABELS,TRIGGER,WIN_LENGTH,COLOR,LINE_WIDTH,YLABEL)    
% function PlotEEG_maison(DATA,MASK,FS,LABELS,WIN_LENGTH)    
% 
% This function plots channels side by side or all in one row ('butterfly' mode).
% It features a slider control on time windows. 
% It is a (very) light version of 'eegplot' function from eeglab toolbox.  
%
% *** Inputs *** 
% - DATA            --> 2D-data [channels x time samples]
% - FS              --> (optional) sampling frequency in Hz
% - LABELS          --> (optional) name of channels in DATA (string array)
% - TRIGGER         --> (optional) trigger, can be a vector with time samples: 
%                                       [t1,t2,...,tn]
%                                  or a 2D-cell array including trigger name:
%                                       TRIGGER{1:n,1}=(t1,...,tn) 
%                                       TRIGGER{1:n,2}=(name1,...name_n)
% - WIN_LENGTH      --> (optional) length of window to be displayed (in seconds)
% - COLOR           --> (optional) '0' for black and white, '1' for color,
%                                   or directly a matrix of [channels , 3] RGB values
%                                   or directly a cell structure of {channels , 3 } RGB values
% - LINE_WIDTH    	--> (optional) vector with line widths to display (channels elements) 
% - YLABEL          --> (optional) string with ylabel
%
%  *** History *** 
% Last version: 16/11/2012
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

% default ylabel
if(nargin < 8)||(isempty(YLABEL))           YLABEL  = [];  end
% default line widths
if(nargin < 7)||(isempty(LINE_WIDTH))      	LINE_WIDTH  = ones(1,size(DATA,1));  end
% default color mode
if(nargin < 6)||(isempty(COLOR))            COLOR  = 1;  end
% default window length
if(nargin < 5)||(isempty(WIN_LENGTH))       WIN_LENGTH  = 10;  end
% default trigger
if(nargin < 4)||(isempty(TRIGGER))          
    display_trigger  = 0;  
else
   display_trigger = 1; 
   if iscell(TRIGGER)   
       TrigVal = cell2mat(TRIGGER(:,1)); 
   else
       TrigVal = TRIGGER;
   end
end
% default labels
if(nargin < 3)||(isempty(LABELS))           
    LABELS = [repmat('data',size(DATA,1),1) int2str((1:size(DATA,1))')];
elseif length(LABELS) ~= size(DATA,1)
    error('Not the same number of labels and data channels');
end
% default frequency sampling
if(nargin < 2)||(isempty(FS))               FS = 128 ;  end


% various init
[Nchan Ns]  = size(DATA);
t_axis      = [0:Ns-1]./FS;
if iscell(COLOR) 
    colors = cell2mat(COLOR);
    if size(COLOR,1) ~= Nchan
        error('Number of elements in color cell must be the same as number of channels');
    end
elseif length(COLOR)~=1
    colors = COLOR;
    if size(COLOR,1) ~= Nchan
        error('Number of rows in color matrix must be the same as number of channels');
    end
elseif COLOR == 0
    colors      = ones(Nchan,1)* [0.4792 0.5625 0.5625];
elseif COLOR == 1
    colors      = hsv(Nchan);
end
min_data    = min(min(DATA));
max_data    = max(max(DATA));
total_win   = round(2*t_axis(end)/WIN_LENGTH);  % total number of time windows
h           = gcf;
set(h,'Toolbar','figure')
pos = get(gca,'Position');
t1 = get(gca,'TightInset');
set(gca,'Position',[1.6*t1(1) pos(2) pos(3)+pos(1)-1.5*t1(1) 1.05*pos(4)]);

win_sel     = 1;        % select first time window
amp_scale   = 1;        % default amplitude scale  
% chan_space  = mean(max(DATA)-min(DATA));        % space between each channel plot
chan_space  = 1/Nchan;  % space between each channel plot
isGrid      = 0;        % default grid display
display_mode= 0;        % '1'->on row per channel , '0' -> all channels in one row

sliderSel1  = uicontrol(    'style','slider','Min',1,'Max',total_win+.01,'SliderStep',[.05/total_win 3/total_win],...
                            'units','normalized','Position',[.35 .02 .33 .03], 'BackgroundColor','w', ...
                            'Value',win_sel, 'Callback',@redraw_cb);
                    
buttonAmpPlus = uicontrol(     'style','pushbutton', 'String', '+', ...
                            'units','normalized','Position',[.96 .15 .03 .03], ...
                            'BackgroundColor','w','Callback',@incAmpScale_cb);

buttonAmpMinus = uicontrol(    'style','pushbutton', 'String', '-', ...
                            'units','normalized','Position',[.92 .15 .03 .03], ...
                            'BackgroundColor','w','Callback',@decAmpScale_cb);
                        
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
                        
buttonAllinOne = uicontrol(    'style','pushbutton', 'String', 'Display Mode', ...
                            'units','normalized','Position',[.92 .29 .07 .03], ...
                            'BackgroundColor','w','Callback',@switchMode_cb);
                                        
                        
% first draw
redraw_cb(h, win_sel);
                    
                    
% callback function called when changing slider value
function redraw_cb(h, win_sel)
    % get new value from slider 
    win_sel = get(sliderSel1, 'value');
    if win_sel > 0
        t2 	= clamp(round(win_sel*WIN_LENGTH*FS),[1 Ns]);
        t1  = clamp(t2-WIN_LENGTH*FS,[1 Ns]);
    else
        t1 = 1; 
        t2 = Ns;
    end
    chan_offsets = chan_space*(Nchan-(1:Nchan)+.5); % offset to display channels on separate rows
    
    cla,            % clear previous plots
    legend off,
    axis manual,    % prevent automatic scaling
    ylim([0 1]),
    hold on,        % save this for all subsequent plots

    % plot signals    
    for chan_ix = 1:Nchan
        % display_mode: '1'-> one row per channel , '0' -> all channels in one row
        if display_mode == 1   
            plot(t_axis(t1:t2),amp_scale*DATA(chan_ix,t1:t2)+chan_offsets(chan_ix),'color',colors(chan_ix,:),'LineWidth',LINE_WIDTH(chan_ix));
            line(get(gca,'XLim'),[chan_offsets(chan_ix) chan_offsets(chan_ix)],'color',[.7 .7 .7],'LineWidth', .5, 'LineStyle', ':')
        else
            plot(t_axis(t1:t2),amp_scale*DATA(chan_ix,t1:t2),'Color',colors(chan_ix,:),'LineWidth',LINE_WIDTH(chan_ix));
            ylim(1.2*[min_data max_data]);
        end
    end

    % draw triggers
    if display_trigger
        trig_inWin = find((TrigVal>t1)&(TrigVal<t2)); % is there triggers in selected window?
        for trig_ix = 1:length(trig_inWin) 
            trig_Xval = TrigVal(trig_inWin(trig_ix))  ; % trigger value in time sample
            line([t_axis(trig_Xval) t_axis(trig_Xval)], get(gca,'YLim'),'color','k','LineWidth', 2, 'LineStyle', '--')
            % plot trigger name
            if iscell(TRIGGER) 
                xPos = t_axis(trig_Xval) - .01*(max(get(gca,'XLim'))-min(get(gca,'XLim')));
                yPos = 1.02*max(get(gca,'YLim'));
                text(xPos,yPos,TRIGGER(trig_inWin(trig_ix),2),'FontSize',8,'FontWeight','bold');
            end
        end
    end
    
    % set axis and labels
    xlabel('Time (s)');
    xlim([t1 t2]/FS);
    if display_mode == 1 
        set(gca,'YTick',chan_offsets(end:-1:1)),
        set(gca,'YTickLabel',LABELS(end:-1:1)),
    else
        set(gca,'YTickMode','auto','YTickLabelMode','auto'),
        scaled_labels = str2num(get(gca,'YTickLabel'))/amp_scale;
        set(gca,'YTickLabel',num2str(scaled_labels))   % correct Ylabel with current amplitude scale
        legend(LABELS,'Location','BestOutside'),
        line(get(gca,'XLim'),[0 0],'color',[.7 .7 .7],'LineWidth', .5,'LineStyle',':'),
        ylabel(YLABEL),
    end    
    box on,
    hold off,
    if isGrid == 1
       grid on 
    end
end


% callback: increase amplitude scale
function incAmpScale_cb(h, buttonAmpPlus)
    amp_scale = amp_scale*2;
    redraw_cb(h, win_sel)
end

% callback: decrease amplitude scale
function decAmpScale_cb(h, buttonAmpMinus)
    amp_scale = amp_scale/2;
    redraw_cb(h, win_sel)
end

% callback: increase time scale
function incTimeScale_cb(h, buttonTimePlus)
    WIN_LENGTH = 2^nextpow2(WIN_LENGTH * 2);
    if WIN_LENGTH > Ns/FS
        WIN_LENGTH = Ns/FS;
        win_sel = 0;
        set(sliderSel1,'Min',0,'Max',1,'value',win_sel,'visible','off');
    else
        total_win   = ceil(t_axis(end)/WIN_LENGTH);  % total number of time windows
        win_sel = clamp(get(sliderSel1, 'value'),[1 total_win]);
        set(sliderSel1,'Max',total_win,'SliderStep',[.05/total_win 3/total_win],'value',win_sel,'visible','on');
    end
    redraw_cb()
end

% callback: decrease time scale
function decTimeScale_cb(h, buttonTimeMinus)
    WIN_LENGTH = 2^nextpow2(WIN_LENGTH / 2);
    if WIN_LENGTH < .5
       WIN_LENGTH = .5; 
    end
    total_win   = ceil(t_axis(end)/WIN_LENGTH);  % total number of time windows
    win_sel = clamp(get(sliderSel1, 'value'),[1 total_win]);
    set(sliderSel1,'Min',1,'Max',total_win,'SliderStep',[.05/total_win 3/total_win],'value',win_sel,'enable','on','visible','on');
    set(buttonTimePlus,'enable','on','visible','on');
    redraw_cb()
end

% callback: switch display mode 
function switchMode_cb(h, buttonPlus)
    display_mode = ~display_mode;
    redraw_cb(h, win_sel)
end

% callback: switch grid display
function switchGrid_cb(h, buttonGridSwitch)
    grid
    isGrid = ~isGrid;
end

function switchTrig_cb(h, buttonTrigSwitch)
    if ~isempty(TRIGGER) 
        display_trigger = ~display_trigger;
        redraw_cb(h, win_sel)
    end
end

function [y] = clamp(x, range)
    y = x;
    y(x<range(1)) = range(1);
    y(x>range(2)) = range(2);
end

end
