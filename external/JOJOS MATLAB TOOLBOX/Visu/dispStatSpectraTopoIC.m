function dispStatSpectraTopoIC(SPECTRA_COND1,SPECTRA_COND2,A,ELEC_LOC,F_AXIS,F_RANGE,TTEST_H,PERCENTILE_COND1,PERCENTILE_COND2,DISPLAY_MODE)                
% function dispStatSpectraTopoIC(SPECTRA_COND1,SPECTRA_COND2,A,ELEC_LOC,F_AXIS,TTEST_H,TTEST_VAR,DISPLAY_MODE)  
%
% This function displays IC scalp maps and spectra computed in two
% different conditions. Using prior t-test results, it shows frequencies 
% that are significantly different, as well as variance of IC spectra.
%
% *** Inputs *** 
% - SPECTRA_COND1	--> IC spectra in condition 1 (nbIC,nbFreq)
% - SPECTRA_COND2  	--> IC spectra in condition 2 (nbIC,nbFreq)
% - A               --> demixing matrix (Nelec,nbIC)
% - ELEC_LOC        --> eeg-lab electrode location 
% - F_AXIS         	--> frequency axis (1,nbFreq)
% - F_RANGE         --> frequency range for spectrum display [Fmin Fmax]
% - TTEST_H         --> (optionnal) H0 t-test rejection for each frequency (nbIC,nbFreq)
% - PERCENTILE_COND1--> (optionnal) inferior and superior percentile of IC spectra in cond 1(nbIC,nbFreq,2)
% - PERCENTILE_COND2--> (optionnal) inferior and superior percentile of IC spectra in cond 1(nbIC,nbFreq,2)
% - DISPLAY_MODE    --> (optionnal) show either all IC at the same time
%                                   ('0' default) or IC one by one ('1')
%
% /!\ if error "Undefined function or method 'topoplot'..."->launch eeglab!
%
%  *** History *** 
% Last version: 24/04/2013
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


% --- assess inputs
if mod(size(SPECTRA_COND1,1),size(A,2))~=0
    error('Not the same number of IC in spectra and demixing matrix!'); 
end
if nargin >= 7
    show_TTEST_H = 1;
else
    show_TTEST_H = 0;
end
if nargin >= 8
    show_SPEC_PER1 = 1;
else
    show_SPEC_PER1 = 0;
end
if nargin >= 9
    show_SPEC_PER2 = 1;
else
    show_SPEC_PER2 = 0;
end
if nargin < 10
    DISPLAY_MODE = 0;
end

% -- various init
FreqsInRange    = find((F_AXIS>=F_RANGE(1))&(F_AXIS<=F_RANGE(2)));
F_AXIS          = F_AXIS(FreqsInRange);
SPECTRA_COND1   = SPECTRA_COND1(:,FreqsInRange);
SPECTRA_COND2   = SPECTRA_COND2(:,FreqsInRange);
nbIC = size(A,2);
ySubPlot = clamp(nbIC, [1 5]);          % max 5 columns
xSubPlot = 2*ceil(nbIC/ySubPlot);       % nb lines (Topo+Spectrum x nbLines)
A = A - repmat(mean(A),size(A,1),1);    % remove mean from each column
ylims1  = 20*log10([min(min(SPECTRA_COND1)) max(max(SPECTRA_COND1))]);
ylims2  = 20*log10([min(min(SPECTRA_COND2)) max(max(SPECTRA_COND2))]);
YLIMS   = [min([ylims1 ylims2]) max([ylims1 ylims2])];
h = figure( 'Name', 'Scalp maps and spectra','Position',get(0,'ScreenSize'), 'NumberTitle','off');
global h_topo
h_topo = [];
% UI objects
sliderSel = uicontrol(  'style','slider','String','IC#','Min',1,'Max',nbIC , ...
                        'SliderStep',[1/nbIC 5/nbIC],'Value',nbIC, ...
                        'units','normalized','Position',[.94 .25 .04 .5], ...
                        'BackgroundColor','w','Callback',@redraw_cb);
buttonSwitchDisplay = uicontrol(    'style','pushbutton', 'String', 'Display Mode', ...
                            'units','normalized','Position',[.93 .15 .06 .05], ...
                            'BackgroundColor','w','Callback',@switchMode_cb);
redraw_cb(h, buttonSwitchDisplay); % first draw

% -- plot source scalp maps and spectra
function redraw_cb(h, buttonSwitchDisplay)
    cla(h_topo),    % allow reseting topoplot
    cla,            % clear previous plots
    legend off,
    axis manual,    % prevent automatic scaling
    hold on,        % save this for all subsequent plots
    
    % show all IC at the same time
    if DISPLAY_MODE == 0
        for IC_ix = 1:nbIC
            posSubPlot = IC_ix + 1*fix((IC_ix-1)./ySubPlot)*ySubPlot;
            % scalp map IC 
            subplot(xSubPlot,ySubPlot,posSubPlot), 
            topoplot(A(:,IC_ix),ELEC_LOC,'electrodes','on');
            % spectrum IC in both conditions
            subplot(xSubPlot,ySubPlot,posSubPlot+ySubPlot), 
            plot(F_AXIS,20*log10(SPECTRA_COND2(IC_ix,:)),'o','MarkerEdgeColor','k','MarkerSize',5,'LineWidth',1),
        	hold on,
            plot(F_AXIS,20*log10(SPECTRA_COND1(IC_ix,:)),'b','linewidth',2),
            title(['IC_{' int2str(IC_ix) '}']), grid on, xlim(F_RANGE), ylim(YLIMS),
            % shows spectrum percentile for first condition
            if show_SPEC_PER1
                plot(F_AXIS,10*log10(squeeze(PERCENTILE_COND1(IC_ix,:,1).^2)),'--b','linewidth',1),
                plot(F_AXIS,10*log10(squeeze(PERCENTILE_COND1(IC_ix,:,2).^2)),'--b','linewidth',1),
            end
            if show_SPEC_PER2
                plot(F_AXIS,10*log10(squeeze(PERCENTILE_COND2(IC_ix,:,1).^2)),':k','linewidth',1),
                plot(F_AXIS,10*log10(squeeze(PERCENTILE_COND2(IC_ix,:,2).^2)),':k','linewidth',1),
            end
            % shows frequencies that are significantly different
            if show_TTEST_H
                Sig_ix = logical(TTEST_H(IC_ix,:));
                plot(F_AXIS(Sig_ix),20*log10(SPECTRA_COND2(IC_ix,Sig_ix)),'o','MarkerEdgeColor','r','MarkerSize',8,'LineWidth',2),
            end
        end
    % show IC one by one
    elseif DISPLAY_MODE == 1
        IC_ix = nbIC-round(get(sliderSel, 'value'))+1;
        % scalp map IC 
        h_topo = subplot(2,1,1);
        topoplot(A(:,IC_ix),ELEC_LOC,'electrodes','on');
        % spectrum IC in both conditions
        subplot(2,1,2);
      	% displays dotted lines for IC spectrum variance at each frequency
%         if show_SPEC_PER1
%             plot(F_AXIS(FreqsInRange),20*log10(SPECTRA_COND1(IC_ix,FreqsInRange)-.5*SPECTRA_VAR(IC_ix,FreqsInRange)),'--k','LineWidth',1),
%             hold on,
%             plot(F_AXIS(FreqsInRange),20*log10(SPECTRA_COND1(IC_ix,FreqsInRange)+.5*SPECTRA_VAR(IC_ix,FreqsInRange)),'--k','LineWidth',1),
%             YLIMS = [	min(min(20*log10(SPECTRA_COND1(:,FreqsInRange))-20*log10(.5*SPECTRA_VAR(:,FreqsInRange)))) , ...
%                      	max(max(20*log10(SPECTRA_COND1(:,FreqsInRange))+20*log10(.5*SPECTRA_VAR(:,FreqsInRange))))];
%         end
        plot(F_AXIS,20*log10(SPECTRA_COND1(IC_ix,:)),'b','linewidth',2),
        hold on ,
        plot(F_AXIS,20*log10(SPECTRA_COND2(IC_ix,:)),'o','MarkerEdgeColor','k','MarkerSize',8,'LineWidth',2),
        title(['IC_{' int2str(IC_ix) '}']), 
%         grid on, 
        xlim(F_RANGE), ylim(YLIMS),
        xlabel('Frequency (Hz)'),
        ylabel('PSD (dB)'),
        % shows spectrum variance in first condition
        if show_SPEC_PER1
            plot(F_AXIS,10*log10(squeeze(PERCENTILE_COND1(IC_ix,:,1).^2)),'--b','linewidth',1),
            plot(F_AXIS,10*log10(squeeze(PERCENTILE_COND1(IC_ix,:,2).^2)),'--b','linewidth',1),
        end
        if show_SPEC_PER2
            plot(F_AXIS,10*log10(squeeze(PERCENTILE_COND2(IC_ix,:,1).^2)),':k','linewidth',1),
            plot(F_AXIS,10*log10(squeeze(PERCENTILE_COND2(IC_ix,:,2).^2)),':k','linewidth',1),
        end
        % shows frequencies that are significantly different
        if show_TTEST_H
            Sig_ix = logical(TTEST_H(IC_ix,:));
            plot(F_AXIS(Sig_ix),20*log10(SPECTRA_COND2(IC_ix,Sig_ix)),'o','MarkerEdgeColor','r','MarkerSize',12,'LineWidth',3),
        end
    end
end

% callback: switch display mode 
function switchMode_cb(h, buttonSwitchDisplay)
    DISPLAY_MODE = ~DISPLAY_MODE;
    if(DISPLAY_MODE == 0)
        set(sliderSel,'Enable','off');
    else
        set(sliderSel,'Enable','on');
        h_topo = [];
    end
    redraw_cb(h, buttonSwitchDisplay)
end


end

    
    