% This script tests the influence of various parameters on Canonical
% Correlation Analysis (CCA) of two datasets band-pass filtered and sliced 
% in a number of time windows. 
% One CCA being computed for each of these time window. 
%
% Tested parameters are:
% - data length (ie. windows size)
% - frequency band width
% - number of variables in each of the two datasets
%
% Last version: 20/07/2012
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

close all
clear all

% Simulation parameters
cfgSIMU.NbChan      = 5;       % number of channels in each dataset
cfgSIMU.WinSize     = 512;      % windows size in number of samples
cfgSIMU.BandWidth   = 4;      % frequency band width in Hz
temp_ix = 1;
disp('Simulation launched.')

% 1) Generation of 2 datasets of independent signals 
cfgSIMU.Ns = 100000;     % total number of samples
cfgSIMU.Fs = 128;       % sampling rate

X1 = randn(cfgSIMU.NbChan, cfgSIMU.Ns);
X2 = randn(cfgSIMU.NbChan, cfgSIMU.Ns);


% 2) Band-pass filtering
H_bp = designBandPassFilter_FIR(32,cfgSIMU.BandWidth); 
X1_filt = filter(H_bp,X1,2); 
X2_filt = filter(H_bp,X2,2); 

channelSel=2;    % plot spectrum of selected channel before/after filtering
figure,
CompareSpectra(X1(channelSel,:),X1_filt(channelSel,:),2048,'blackman',.75,cfgSIMU.Fs); 
clear channelSel H_bp


% 3) Processing CCA on time windows
[cca.A,cca.B,cca.R,cca.stats,cca.t,cca.masked] = proc_CCA_win([X1_filt ; X2_filt]', [cfgSIMU.NbChan cfgSIMU.NbChan], cfgSIMU.WinSize);
Chisq 	= zeros(cfgSIMU.NbChan,length(cca.stats));
pChisq	= zeros(cfgSIMU.NbChan,length(cca.stats));
for win_ix = 1:length(cca.stats)
    if ~isempty(cca.stats{1,win_ix})
        Chisq(:,win_ix)     = cca.stats{1,win_ix}.chisq;
        pChisq(:,win_ix)    = cca.stats{1,win_ix}.pChisq;
    else
        Chisq(:,win_ix) = 0;
        pChisq(:,win_ix) = 0;
    end
end
clear win_ix

% 4) display mean(R) for all time windows along with significance levels
figure,
R_mean          = mean(cca.R(cca.masked,:),1);
pChisq_mean     = mean(pChisq(:,cca.masked),2)';
p_levels        = [0.001 0.01 0.1];
% figure('Name', ['MEAN canonical correlations in couple''s EEG'], 'NumberTitle','off','WindowStyle', 'docked'),
subplot(4,1,temp_ix),
bar(R_mean,pChisq_mean), %barPlot(R_mean,pChisq_mean,p_levels),
hold on,
grid on,
xlabel('CCA components'),
ylabel('canonical correlation')
ylim([0 1.1])
title(['time windows = ' int2str(cfgSIMU.WinSize) ' samples   /   BP = ' int2str(cfgSIMU.BandWidth) 'Hz   /   f_0 = ' int2str(32) 'Hz'], 'FontSize', 10);
disp('Simulation done.')
clear channelSel AX H1 H2


