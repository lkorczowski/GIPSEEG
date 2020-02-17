function compScalpMapSpectra(FAXIS,NB_IC,F_RANGE,ELEC_LOC, A_SUB1,SPECTRUM_IC,A_SUB2)    
% function CompScalpMapSpectra(FAXIS,NB_IC,F_RANGE,ELEC_LOC, A_SUB1,SPECTRUM_IC,A_SUB2)   
% 
% This function plots components scalp maps along with their spectra for
% one or two subjects (the more components the tinier!).
% /!\ This function uses EEGLAB's "topoplot" function.
%
% *** Inputs *** 
% - FAXIS           --> frequency axis (1,f)
% - NB_IC           --> number of components to plot for each subject
% - F_RANGE         --> frequency range for spectrum display [Fmin Fmax]
% - ELEC_LOC        --> eeg-lab electrode location 
% - A_SUB1         	--> demixing matrix for subject 1 (Nelec,Ncomp)
% - SPECTRUM_IC     --> matrix of components spectra 
%                       (N*Ncomp,f). Amplitudes must be in linear scale.
% - A_SUB2         	--> (optionnal) demixing matrix for subject 2 (Nelec,Ncomp)
%
% if N>1 and A_SUB2 is not given, then N spectra are displayed for each
% component scalp map (e.g. when estimated IC spectra in different conditions).
%
% /!\ if error "Undefined function or method 'topoplot'..."->launch eeglab!
%
%  *** History *** 
% First version: 24/04/2013
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

if mod(size(SPECTRUM_IC,1),size(A_SUB1,2))~=0
    error('Not the same number of IC in spectrum and demixing matrix!'); 
end

% -- plot source scalp maps and spectra
FreqsInRange = find((FAXIS>=F_RANGE(1))&(FAXIS<=F_RANGE(2)));
N = size(SPECTRUM_IC,1)/size(A_SUB1,2);
ySubPlot = clamp(NB_IC, [1 10]);    % max 10 columns
A_SUB1 = A_SUB1 - repmat(mean(A_SUB1),size(A_SUB1,1),1); % remove mean from each column
YLIMS = 20*log10([min(min(SPECTRUM_IC(:,FreqsInRange))) max(max(SPECTRUM_IC(:,FreqsInRange)))]);
figure( 'Name', 'Scalp maps and spectra','Position',get(0,'ScreenSize'), 'NumberTitle','off');

% only one subject
if nargin <= 6 
    xSubPlot = ceil(NB_IC/ySubPlot)*(1+N);  % (1+N) plots (1 map + N spectra/map)
    for IC_ix = 1:NB_IC
        posSubPlot = IC_ix + N*fix((IC_ix-1)./ySubPlot)*ySubPlot;
        % scalp map IC 
        subplot(xSubPlot,ySubPlot,posSubPlot), 
        topoplot(A_SUB1(:,IC_ix),ELEC_LOC,'electrodes','on');
        % spectrum IC 
        for N_ix = 1:N
            subplot(xSubPlot,ySubPlot,posSubPlot+(N_ix*ySubPlot)), 
            plot(FAXIS(FreqsInRange),20*log10(SPECTRUM_IC((N_ix-1)*NB_IC+IC_ix,FreqsInRange))),
            title(['IC_{' int2str(IC_ix) '}']), grid on, xlim(F_RANGE), ylim(YLIMS),
        end
    end
% two subjects
else
    A_SUB2 = A_SUB2 - repmat(mean(A_SUB2),size(A_SUB1,1),1);
    xSubPlot = 4*ceil(NB_IC/ySubPlot);  % 2 subjects x 2 plots (map + spectrum)
    for IC_ix = 1:NB_IC
        posSubPlot = IC_ix + 3*fix((IC_ix-1)./ySubPlot)*ySubPlot;
        % scalp map IC subject 1
        subplot(xSubPlot,ySubPlot,posSubPlot), 
        topoplot(A_SUB1(:,IC_ix),ELEC_LOC,'electrodes','on');
        % scalp map IC subject 2
        subplot(xSubPlot,ySubPlot,posSubPlot+2*ySubPlot), 
        topoplot(A_SUB2(:,IC_ix),ELEC_LOC,'electrodes','on');
        % spectrum IC subject 1
        subplot(xSubPlot,ySubPlot,posSubPlot+ySubPlot), 
        plot(FAXIS(FreqsInRange),20*log10(SPECTRUM_IC(IC_ix,FreqsInRange))),
        title(['S1 - IC_{' int2str(IC_ix) '}']), grid on, xlim(F_RANGE), ylim(YLIMS),
        % spectrum IC subject 2
        subplot(xSubPlot,ySubPlot,posSubPlot+3*ySubPlot), 
        plot(FAXIS(FreqsInRange),20*log10(SPECTRUM_IC(IC_ix+size(A_SUB1,2),FreqsInRange))),
        title(['S2 - IC_{' int2str(IC_ix) '}']),grid on, xlim(F_RANGE), ylim(YLIMS),
    end
end
    
    