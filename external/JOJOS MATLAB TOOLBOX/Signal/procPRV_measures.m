function [PRV_measures] = procPRV_measures(PRV, Fs, display_mode)
% function [PRV_measures] = procPRV_measures(PRV, Fs, display_mode)
%
% ************************************************************************
% Process common measures from Pulse Rate Variability (PRV) signal
%
% Input:
% ------
% PRV           --> input timeserie (nbChan x nbSamples) (s)
% Fs            --> (optional) PRV frequency rate (Hz)
% display_mode  --> (optional) if set to '1', display PSD
%
% Output:
% -------
% PRV_measures -> structure with fields:
% MEAN          --> Mean value in considered time interval (nbChan x 1)
% SD            --> Standard deviation of all P-P intervals in ms (nbChan x 1)
% HF            --> Power in high frequency range 0.15–0.4 Hz in ms^2 (nbChan x 1)
% LF            --> Power in low frequency range 0.04-0.15 Hz in ms^2 (nbChan x 1)
% PSD           --> PRV power spectrum density in ms^2/Hz (nbChan x 1)
% PSD_faxis     --> Frequency axis for power spectrum    
%
% History:
% --------
% Last version:  2013-06-10
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com
% based on Task force of the European Society of Cardiology and the North 
% American Society of Pacing and Electrophysiology, 1996.

% *** assess inputs
if (nargin < 3) || isempty(display_mode)
    display_mode = 0;
end
if (nargin < 2) || isempty(Fs)
    Fs = 128;
end
sig_length  = size(PRV,2)/Fs; % signal duration in seconds
nb_chan     = size(PRV,1);
if( sig_length < 60 )
    error('At least 1 minute of PRV is required to process time/freq variables!')
end

% *** process time-domain measure of PRV
LF_min_freq = .04;  % in Hz
HF_min_freq = .15;  % in Hz
HF_max_freq = .4;   % in Hz
PRV_measures.MEAN = mean(PRV,2);
PRV_measures.SD = std(PRV,[],2);

% *** now decimate PRV 
% (this has a direct impact on subsequent frequency definition)
if( sig_length > 150 )
    dec_factor = 16;   
else
    dec_factor = 6;    
end
F_dec   = Fs/(dec_factor);
PRV_dec = resample(PRV',1,dec_factor)'; 

% *** process PSD
win_length  = 1024;
overlap     = .75;
PRV_measures.PSD = zeros(nb_chan, win_length+1);
for chan_ix = 1:nb_chan
    % - parametric spectral estimation: AR
    % [PRV_measures.PSD(chan_ix,:), PRV_measures.PSD_faxis] = pyulear(PRV_dec(chan_ix,:), 10, 2*win_length, F_dec);
    % - non-parametric spectral estimation: Welch periodogram
    [PRV_measures.PSD(chan_ix,:), PRV_measures.PSD_faxis] = pwelch(PRV_dec(chan_ix,:),win_length,overlap*win_length,2*win_length,F_dec);
end

% ***  process frequency-domain measure of PRV
VLF_ix  = find((PRV_measures.PSD_faxis<LF_min_freq));
LF_ix   = find((PRV_measures.PSD_faxis>=LF_min_freq)&(PRV_measures.PSD_faxis<=HF_min_freq));
HF_ix   = find((PRV_measures.PSD_faxis>=HF_min_freq)&(PRV_measures.PSD_faxis<=HF_max_freq));
PRV_measures.VLF = sum(abs(PRV_measures.PSD(:,VLF_ix)),2);
PRV_measures.LF = sum(abs(PRV_measures.PSD(:,LF_ix)),2);
PRV_measures.HF = sum(abs(PRV_measures.PSD(:,HF_ix)),2);


% *** display PSD
if display_mode == 1
    chan_sel = 12;
    figure('WindowStyle','docked'),
    LF_ix = [max(VLF_ix);LF_ix];% only for display purpose
    HF_ix = [max(LF_ix);HF_ix]; % only for display purpose
    area(PRV_measures.PSD_faxis(VLF_ix),PRV_measures.PSD(chan_sel,VLF_ix),'FaceColor',[250 250 250]./255,'EdgeColor','none'), hold on,
    area(PRV_measures.PSD_faxis(LF_ix),PRV_measures.PSD(chan_sel,LF_ix),'FaceColor',[50 50 50]./255,'EdgeColor','none'),
    area(PRV_measures.PSD_faxis(HF_ix),PRV_measures.PSD(chan_sel,HF_ix),'FaceColor',[150 150 150]./255,'EdgeColor','none'),
    legend('VLF','LF','HF'),
    plot(PRV_measures.PSD_faxis,PRV_measures.PSD(chan_sel,:),'color','k','linewidth',2),
    xlabel('Frequency (Hz)'), ylabel('PSD (s²/Hz)'), title('PRV power spectrum'),
    xlim([0 .5]), ylim([0 .05]),
    grid on, set(gca,'box','on'),
end









