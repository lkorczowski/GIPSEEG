function [Hd b a] = designBandPassFilter_FIR(FO,BP)
%DESIGNBANDPASSFILTER_FIR(FO,BP) Returns a discrete-time filter object.
%
% Equiripple Bandpass filter designed using the FIRPM function.

% When BP is not specified, this filter will always be of order 118 with
% band-pass of 4Hz. Otherwise order will increase as BP decreases.

% *** Inputs *** 
% - FO      --> central frequency, in Hz.
% - BP      --> (optional) Band-pass, in Hz. Default value: 4Hz



% default BP
if(nargin < 2)||(isempty(BP))   BP = 4 ;  end


% All frequency values are in Hz.
Fs = 128;  % Sampling Frequency

Fstop1 = FO-1*BP;                 % First Stopband Frequency
Fpass1 = FO-0.5*BP;               % First Passband Frequency
Fpass2 = FO+0.5*BP;               % Second Passband Frequency
Fstop2 = FO+1*BP;                 % Second Stopband Frequency
Dstop1 = 0.001;                 % First Stopband Attenuation
Dpass  = 0.0057563991496;       % Passband Ripple
Dstop2 = 0.001;                 % Second Stopband Attenuation
dens   = 20;                    % Density Factor

if  Fstop1<0     Fstop1 = 0;    end
if  Fstop2>Fs/2  Fstop2 = Fs/2;	end


% Calculate the order from the parameters using FIRPMORD.
if  (Fpass1<=0) || (Fpass2>=Fs/2)    
    Hd = dfilt.scalar;
else    
    [N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
                              0], [Dstop1 Dpass Dstop2]);
    b  = firpm(N, Fo, Ao, W, {dens});
    a = 1;
    Hd = dfilt.df2(b,a);
end

% [EOF]
