function [PRV_norm, PRV_trend] = procPRV(PPG, Fs, PPI_method, PPI_min)
% function  [PRV_norm, PRV_trend] = procPRV(PPG)
%
% ************************************************************************
% Calculate Pulse Rate Variability (PRV) from raw (unfiltered) 
% photoplethysmographic signal (PPG)
%
% Input:
% ------
% PPG           --> input timeseries (nbSig x nbSamples)
% Fs            --> (optional)  frequency rate (Hz)
% PPI_method    --> (optional)  method used to estimate Pulse to Pulse 
%                               Interval (PPI)
%                   'peak' for using beginning of the anacrotic phase
%                   'foot' for using beginning of the catacrotic phase
%                   'derivative1' for using steepest part of the upstroke 
%                   'derivative2' for using max of 2nd PPG derivative 
% PPI_min       --> (optional) min PPI, i.e. max pulsation beat (bpm)
%
% Output:
% -------
% PRV_norm   	--> 0_mean, low-pass filtered PRV (nbSig x nbSamples). unity: s. 
% PRV_trend 	--> PRV trend, obtained with polynomial curve fitting
%
% History:
% --------
% Last version:  2013-05-27
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com
% based on review from (Schäfer and Vagedes, 2012), and work from (Gil et
% al., 2010)


if (nargin < 4) || isempty(PPI_min)
    PPI_min = 120;
end
if (nargin < 3) || isempty(PPI_method)
    PPI_method = 'foot';
end
if (nargin < 2) || isempty(Fs)
    Fs = 128;
end
nbSig = size(PPG,1);
nbSamples = size(PPG,2);
PRV_norm = zeros(nbSig,nbSamples);
if nargout>1
    PRV_trend = zeros(nbSig,nbSamples);
end


% resample PPG at 32Hz
Fs_dec = 32; % frequency rate after decimation
PPG_dec = resample(double(PPG'),1,Fs/Fs_dec)';


% remove mean and high-pass filter PPG 
PPG_hp = PPG_dec - repmat(mean(PPG_dec,2),[1 size(PPG_dec,2)]);
% PPG_hp = filter_HP(PPG_hp,Fs);


% measure PPI
PPI_min = (60*Fs_dec)/PPI_min;  % convert PPI_min from ppm to samples

switch PPI_method 
    case 'peak'
        PPG_detect = PPG_hp;
    case 'foot'
        PPG_detect = -PPG_hp;
    case 'derivative1'
        PPG_detect = diff(PPG_hp,1,2);
    case 'derivative2'
        PPG_detect = diff(PPG_hp,2,2);
    otherwise
        error('[procPRV] Wrong string for input parameter PPI_method')
end        
parfor sig_ix = 1:nbSig
    % find PPI (Pulse-to-Pulse Interval)
    [pks{sig_ix},locs{sig_ix}] = findpeaks(PPG_detect(sig_ix,:),'minpeakdistance',PPI_min);
    pks{sig_ix} = PPG_hp(sig_ix,locs{sig_ix});

    % process PRV from PPI
    interval = diff(locs{sig_ix},1,2) ./ Fs_dec; % unity: seconds. 

    % remove mean and low-pass filter
    interval_norm = interval - repmat(mean(interval,2),[1 size(interval,2)]);
%     interval_norm = interval_norm ./ std(interval_norm);
	interval_norm = filter_LP(interval_norm,Fs*size(interval_norm,2)/nbSamples);
    
    % interpolation to obtain as many time samples as input signal
    axis_interp	= 1:nbSamples;
    axis_rho 	= 1:nbSamples/length(interval):nbSamples;
    PRV_norm(sig_ix,:) = spline(axis_rho,interval_norm,axis_interp); % pchip or spline interpolation
    
    if nargout>1
        [p,S,mu] = polyfit(1:length(interval),interval,12);    % extract trend with polynomial fit of order 10
        interval_trend = polyval(p,1:length(interval),S,mu);
        PRV_trend(sig_ix,:) = spline(axis_rho,interval_trend,axis_interp); % pchip or spline interpolation
        % display tachograms
%         figure, plot(interval,'b'),hold on, plot(interval_trend,'k'),
%         figure, plot(PRV_norm(sig_ix,:),'k'), xlabel('beat #'), hold on, plot(PRV_trend(sig_ix,:),'b'), xlabel('beat #'), ylabel('P-P interval (s)'), 
    end
end
end

% display PPI
% figure,
% sig_plot = 1; 
% t_axis = (1:size(PPG_hp,2)) ./ (Fs);
% plot(t_axis,PPG_hp(sig_plot,:),t_axis(locs{sig_plot}),pks{sig_plot},'rv','MarkerFaceColor','r'),
% xlim([63 65]),


% design high-pass filter
function out = filter_HP(in,fs)
  	Fstop = .05;     	% Stopband Frequency
    Fpass = 2;          % Passband Frequency
    Astop = 80;      	% Stopband Attenuation (dB)
    Apass = 1;      	% Passband Ripple (dB)
    h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, fs);
    Hd = design(h, 'butter', 'MatchExactly', 'stopband');
    out = filtfilt(Hd.sosMatrix,Hd.ScaleValues,in')';
%     fvtool(b,a,'Fs',fs,'Analysis','freq','Legend','on') 
end

% design low-pass filter 
function out = filter_LP(in,fs)
%     Fpass = .1;        % Passband Frequency % test
%     Fstop = .4;        % Stopband Frequency % test
    Fpass = .05;        % Passband Frequency  % values used in 2013 FiBN paper
    Fstop = .2;        % Stopband Frequency   % values used in 2013 FiBN paper
    Astop = 80;          % Stopband Attenuation (dB)
    Apass = 1;           % Passband Ripple (dB)
    h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, fs);
    Hd = design(h, 'butter', 'MatchExactly', 'stopband');
   	out = filtfilt(Hd.sosMatrix,Hd.ScaleValues,in')';
% 	fvtool(b,a,'Fs',fs,'Analysis','freq','Legend','on') 
end
