function [RVT, RVT_trend] = procRVT(RESP, Fs, BBI_min, visu_debug)
% function  [RVT, RVT_trend] = procRVT(RESP, Fs, BBI_min, visu_debug)
%
% ************************************************************************
% Calculate Respiration Volume per Time (RVT) from raw (unfiltered) signal 
% obtained using a piezo-electric respiration sensor.
%
% This measure is estimated as the difference between the max and min 
% at the peaks of inspiration and expiration. This difference is further 
% divided by the duration of the respiration resulting in the RVT. 
% Finally, the RVT time series is then interpolated to be at the sampling 
% rate of physiological signal
%
% Input:
% ------
% RESP          --> input timeseries (nbSig x nbSamples)
% Fs            --> (optional)  frequency rate (Hz)
% BBI_min       --> (optional) min Breath-to-Breath Interval (in bpm)
% visu_debug    --> (optional) various visualizations
%
% Output:
% -------
% RVT           --> 0_mean, unfiltered RVT signal (nbSig x nbSamples)
% RVT_trend     --> RVT trend, obtained with polynomial curve fitting
%
% History:
% --------
% Last version:  2013-06-05
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com
% inspired from methodology used in (Yuan et al., NeuroImage, 2013)

if (nargin < 4) || isempty(visu_debug)
    visu_debug = 0;
end
if (nargin < 3) || isempty(BBI_range)
    BBI_min = 80;
end
if (nargin < 2) || isempty(Fs)
    Fs = 128;
end
nbSig = size(RESP,1);
nbSamples = size(RESP,2);

% *** resample RESP at 32Hz
Fs_dec = 32; % frequency rate after decimation
RESP_dec = resample(double(RESP'),1,Fs/Fs_dec)';

% *** remove mean 
RESP_hp = RESP_dec - repmat(mean(RESP_dec,2),[1 size(RESP_dec,2)]);

% *** measure minimum and maximum peak positions
BBI_min = (60*Fs_dec)/BBI_min;  % convert BBI_range from ppm to samples       
for sig_ix = 1:nbSig
    % find min/max peaks
    [max_pks{sig_ix},max_locs{sig_ix}] = findpeaks(RESP_hp(sig_ix,:),'minpeakdistance',BBI_min);
    max_pks{sig_ix}     = RESP_hp(sig_ix,max_locs{sig_ix});
  	[min_pks{sig_ix},min_locs{sig_ix}] = findpeaks(-RESP_hp(sig_ix,:),'minpeakdistance',BBI_min);
    min_pks{sig_ix}     = RESP_hp(sig_ix,min_locs{sig_ix});
    
    % treat the particular case when having two successive max/min peaks (remove the 2nd ones)
    all_locs            = sort([min_locs{sig_ix} max_locs{sig_ix}],'ascend');
    loc_members_diff    = diff(ismember(all_locs,min_locs{sig_ix}));
    rep_loc             = (loc_members_diff == 0);
%     max_locs{sig_ix}(ismember(max_locs{sig_ix},all_locs(logical([0 rep_loc]))))/128; % DEBUG: this gives pos of successive maxima
    min_pks{sig_ix}     = min_pks{sig_ix} (~ismember(min_locs{sig_ix},all_locs(logical([0 rep_loc]))));
    max_pks{sig_ix}     = max_pks{sig_ix} (~ismember(max_locs{sig_ix},all_locs(logical([0 rep_loc]))));
    min_locs{sig_ix}    = min_locs{sig_ix}(~ismember(min_locs{sig_ix},all_locs(logical([0 rep_loc]))));
    max_locs{sig_ix}    = max_locs{sig_ix}(~ismember(max_locs{sig_ix},all_locs(logical([0 rep_loc]))));
    
    % crop peak vectors to have same number of maxima and minima
    minimum_length      = min([length(min_pks{sig_ix}),length(max_pks{sig_ix})]);
    min_pks{sig_ix}     = min_pks{sig_ix}(1:minimum_length);
    max_pks{sig_ix}     = max_pks{sig_ix}(1:minimum_length);
    min_locs{sig_ix}    = min_locs{sig_ix}(1:minimum_length);
    max_locs{sig_ix}    = max_locs{sig_ix}(1:minimum_length);
end


% *** process RVT from peak amplitudes and locations
RVT         = zeros(nbSig,nbSamples);
RVT_trend   = zeros(nbSig,nbSamples);

for sig_ix = 1:nbSig
    
    % always perform max-(successive(min)) and not the opposite
	if min_locs{sig_ix}(1) < max_locs{sig_ix}(1)
        min_pks{sig_ix} = min_pks{sig_ix}(2:end); % remove first element
        min_locs{sig_ix}= min_locs{sig_ix}(2:end); % remove first element
        max_pks{sig_ix} = max_pks{sig_ix}(1:end-1); % remove last element
        max_locs{sig_ix}= max_locs{sig_ix}(1:end-1); % remove last element
    end
    
    % PROCESS RVT
    RVT_noInterp = (max_pks{sig_ix}-min_pks{sig_ix}) ./ (-(max_locs{sig_ix}-min_locs{sig_ix}));
        
    % remove mean and outliers
    RVT_filt = RVT_noInterp - repmat(mean(RVT_noInterp,2),[1 size(RVT_noInterp,2)]);
    
    % interpolation to obtain as many time samples as input signal
    axis_interp	= 1:nbSamples;
    axis_rho 	= 1:nbSamples/length(RVT_noInterp):nbSamples;
    RVT(sig_ix,:) = pchip(axis_rho,RVT_filt,axis_interp); % pchip or spline interpolation
    
    % Extract trend
    if nargout>1
        [p,S,mu] = polyfit(1:length(RVT_noInterp),RVT_noInterp,12);    % extract trend with polynomial fit of order 10
        RVT_trend_noInterp  = polyval(p,1:length(RVT_noInterp),S,mu);
        RVT_trend(sig_ix,:) = pchip(axis_rho,RVT_trend_noInterp,axis_interp); % pchip or spline interpolation
    end
end


% *** DEBUG: display peak localization, RVT, RVT trend and histogram 
if (visu_debug == 1)
    figure('name', 'Resp, peaks and RVT','WindowStyle','docked','NumberTitle','off'),
    sig_plot = 1; 
    t_axis = (1:size(RESP_hp,2)) ./ (Fs_dec);
    RVT_dec = resample(double(RVT'),1,Fs/Fs_dec)'; 
    RVT_dec = RVT_dec .* (max(abs(RESP_hp(sig_plot,:)))/max(abs(RVT_dec))); % normalization
    plot(t_axis,RESP_hp(sig_plot,:)),
    hold on, 
    plot(t_axis,RVT_dec(sig_plot,:),'color','k','linewidth',2),
    legend('Resp Signal','RVT'),
    plot(t_axis(max_locs{sig_plot}),max_pks{sig_plot},'rv','MarkerFaceColor','r'),
    plot(t_axis(min_locs{sig_plot}),min_pks{sig_plot},'gv','MarkerFaceColor','g'),
    title('Respiration signal with peaks, and RVT (normalized for better visualization)'),
    xlabel('time (s)'),% xlim([155 180]),
    if nargout>1
        figure('name', 'RVT trend','WindowStyle','docked','NumberTitle','off'),
        t_axis = (1:size(RVT_trend,2)) ./ (Fs);
        plot(t_axis,RVT(sig_ix,:)+mean(RVT_noInterp),'color','b','linewidth',2),
        hold on, plot(t_axis,RVT_trend(sig_ix,:),'color','k','linewidth',2), 
        legend('RVT','RVT trend'), xlabel('time (s)'),
    end
    figure('name', 'RVT distribution','WindowStyle','docked','NumberTitle','off'),
    hist(RVT,10000),
end

end








% % home-made function to crop outliers, based on deviation from median
% % outliers are defined and cropped when higher than x deviations.
% % RVT(sig_ix,:) = cropOutliers(RVT(sig_ix,:)); % remove potential outliers due to interpolation
% function cropped_data = cropOutliers(data)
%     med = median(data);
%     max_dev = 3*std(data);
%     max_crop = med + max_dev;
%     min_crop = med - max_dev;
%     cropped_data = data;
%     cropped_data(data>max_crop) = max_crop;
%     cropped_data(data<min_crop) = min_crop;
% end


% % design low-pass filter
% % RVT_filt = filter_LP(RVT_filt,Fs*size(RVT_filt,2)/nbSamples);
% function out = filter_LP(in,fs)
%     LPfilt_f = [.05, .1]; % stopband and passband frequencies [.02, .05]
%     LPfilt_d = [0.005756, 1e-005];    % stopband attenuation and passband ripple
%     [N, Fo, Ao, W] = firpmord(LPfilt_f/(fs/2), [1 0], LPfilt_d);
%     b  = firpm(N, Fo, Ao, W, {20});
%     a = 1;
%     out = filtfilt(b,a,in')';
% % 	fvtool(b,a,'Fs',fs,'Analysis','freq','Legend','on') 
% end

