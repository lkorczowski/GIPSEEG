function [Data_out] = phaseShiftData(Data_in, Fs, F, Pshift)
% This function shifts a time serie with a given phase shift at a given
% frequency.
%
% INPUTS
% - Data_in 	--> univariate or multivariate time serie (N channels by T time samples)
% - Fs          --> Sampling frequency (Hz)
% - F           --> Frequency at which the phase shift is applied (Hz)
% - Pshift      --> Phase shift (degrees), within range [-360,+360]
%
% OUTPUTS
% - Data_out    --> Phase shifted data (N channels by T time samples)
%                   OUTPUT DATA IS ZERO-PADDED:
%                       - at the end if positive phase shift
%                       - at the begin if negative phase shift
%               
% HISTORY 
% First version:  23/11/2011
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

% default precision on phase shift
eps = 1e-03;

% default UI display
if(nargin < 4) error('[phaseShiftData] wrong arguments');  end

% Phase shift direction
P_dir = sign(Pshift);
Pshift = mod(abs(Pshift),360);
if(Pshift == 0) 
    Data_out = Data_in;
    return
end

% optimal num of samples to shift
d_samples = (Fs*Pshift) / (F*360); 

% if this is not necessary to resample, shift data directly
if(fix(d_samples)==d_samples)
    if(P_dir == 1)
        Data_out = [Data_in(:,d_samples:end) , zeros(size(Data_in,1),d_samples-1)];
    else
        Data_out = [zeros(size(Data_in,1),d_samples) , Data_in(:,1:end-d_samples-1)];
    end
else
    % new sampling rate allowing for a phase shift with precision eps
    d_samples = 0;
    Fs_up = Fs;
    % we don't want to loose information in downsampling data!
    while( Fs_up < Fs  ||  abs((d_samples*F*360/Fs_up)-Pshift) > eps  ||  d_samples == 0)
        d_samples = d_samples+1;
        Fs_up = fix((d_samples*F*360)/Pshift);
    end
    % interpolate data
    Data_temp = resample(Data_in,Fs_up,Fs);
    % shift interpolated data
    if(P_dir == 1)
        Data_temp = [Data_temp(:,1+d_samples:end), zeros(size(Data_temp,1),d_samples-1)];
    else
        Data_temp = [zeros(size(Data_temp,1),d_samples) , Data_temp(:,1:end-d_samples-1)];
    end
    % decimate data to the initial sample rate
    Data_out =  resample(Data_temp,Fs,Fs_up);
end




