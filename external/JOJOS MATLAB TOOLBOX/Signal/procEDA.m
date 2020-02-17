function [EDA, EDA_trend] = procEDA(EDA_RAW, Fs, visu_debug)
% function  [EDA, EDA_trend] = procEDA(EDA_RAW, Fs, visu_debug)
%
% ************************************************************************
% Calculate Electrodermal Activity (EDA) from raw (unfiltered) signal 
% obtained using a galvanic sensor.
%
% Mean of EDA time series is removed, then they are normalized for 
% cross-subject comparison.
%
% Input:
% ------
% EDA_RAW          --> input timeseries (nbSig x nbSamples)
% Fs            --> (optional)  frequency rate (Hz)
% visu_debug    --> (optional) various visualizations
%
% Output:
% -------
% EDA           --> 0_mean, unfiltered EDA signal (nbSig x nbSamples)
% EDA_trend     --> EDA trend, obtained with polynomial curve fitting
%
% History:
% --------
% Last version:  2013-06-05
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com
% inspired from methodology used in (Yuan et al., NeuroImage, 2013)

if (nargin < 3) || isempty(visu_debug)
    visu_debug = 0;
end
if (nargin < 2) || isempty(Fs)
    Fs = 128;
end
nbSig       = size(EDA_RAW,1);
nbSamples   = size(EDA_RAW,2);
EDA_trend   = zeros(nbSig,nbSamples);
EDA_RAW     = double(EDA_RAW);

% *** Remove mean and normalize
EDA = EDA_RAW - repmat(mean(EDA_RAW,2),[1 nbSamples]);
% EDA = filter_HP(EDA,Fs);
EDA = EDA ./ repmat(max(abs(EDA),[],2),1,nbSamples); % normalize with higher min or max

% *** Extract trend --> see if we can remove the loop!
if nargout>1
    for sig_ix = 1:nbSig
        [p,S,mu]    = polyfit(1:nbSamples,EDA_RAW(sig_ix,:),12);    % extract trend with polynomial fit of order 12
        EDA_trend(sig_ix,:) = polyval(p,1:nbSamples,S,mu);
    end
    EDA_trend = EDA_trend ./ repmat(max(abs(EDA),[],2),1,nbSamples);   % normalize with higher min or max
end

% % *** Resample and high-pass filter
Fs_dec = 32; % frequency rate after decimation
EDA_dec = resample(double(EDA'),1,Fs/Fs_dec)';
EDA_dec = filter_HP(EDA_dec,Fs_dec);

% interpolation to obtain as many time samples as input signal
axis_interp	= 1:nbSamples;
axis_rho 	= 1:nbSamples/length(EDA_dec):nbSamples;
EDA = pchip(axis_rho,EDA_dec,axis_interp); % pchip or spline interpolation

% *** DEBUG: display EDA, EDA trend and histogram 
if (visu_debug == 1)
 	if nargout>1
        figure('name', 'EDA trend','WindowStyle','docked','NumberTitle','off'),
        t_axis = (1:size(EDA_trend,2)) ./ (Fs);
        plot(t_axis,EDA(sig_ix,:)+mean(EDA_RAW),'color','b','linewidth',2),
%         hold on, plot(t_axis,EDA_trend(sig_ix,:),'color','k','linewidth',2), 
%         legend('EDA','EDA trend'), xlabel('time (s)'),
    end
    figure('name', 'EDA distribution','WindowStyle','docked','NumberTitle','off'),
    hist(EDA(1,:),10000),
end

end


% design high-pass filter
function out = filter_HP(in,fs)
%     Fstop = .0005;        % Stopband Frequency (value used in obsolet analysis)
%     Fpass = .0025;        % Passband Frequency (value used in obsolet analysis)
%     Fstop = .01;        % Stopband Frequency
%     Fpass = .05;        % Passband Frequency
    Fstop = .005;     	% Stopband Frequency  	% values used in 2013 FiBN paper
    Fpass = .02;          % Passband Frequency   	% values used in 2013 FiBN paper
    Astop = 80;      	% Stopband Attenuation (dB)
    Apass = 1;      	% Passband Ripple (dB)
    h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, fs);
    Hd = design(h, 'butter', 'MatchExactly', 'stopband');
    out = filtfilt(Hd.sosMatrix,Hd.ScaleValues,in')';
%     fvtool(Hd,'Fs',fs,'Analysis','freq','Legend','on') 
end



% design low-pass filter
function out = filter_LP(in,fs)
    Fpass = .1;        % Passband Frequency
    Fstop = .2;        % Stopband Frequency
    Astop = 80;       	% Stopband Attenuation (dB)
    Apass = 1;         	% Passband Ripple (dB)
    h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, fs);
    Hd = design(h, 'butter', 'MatchExactly', 'stopband');
    out = filtfilt(Hd.sosMatrix,Hd.ScaleValues,in')';
%     fvtool(Hd,'Fs',fs,'Analysis','freq','Legend','on') 
end



% % home-made function to crop outliers, based on deviation from median
% % outliers are defined and cropped when higher than x deviations.
% % EDA(sig_ix,:) = cropOutliers(EDA(sig_ix,:)); % remove potential outliers due to interpolation
% function cropped_data = cropOutliers(data)
%     med = median(data);
%     max_dev = 3*std(data);
%     max_crop = med + max_dev;
%     min_crop = med - max_dev;
%     cropped_data = data;
%     cropped_data(data>max_crop) = max_crop;
%     cropped_data(data<min_crop) = min_crop;
% end










