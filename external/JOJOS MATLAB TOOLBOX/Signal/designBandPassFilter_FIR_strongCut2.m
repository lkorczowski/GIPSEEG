function [Hd b a] = designBandPassFilter_FIR_strongCut2(Fstop1,Fstop2,Fs)
%DESIGNBANDPASSFILTER_FIR_STRONGCUT Returns a discrete-time filter object.

% Equiripple Bandpass filter designed using the FIRPM function.
% All frequency values are in Hz.

% default Fs
if(nargin < 3)||(isempty(Fs))   Fs = 128 ;  end


BP = Fstop2-Fstop1;
Fpass1 = Fstop1+0.02*BP;             % First Passband Frequency
Fpass2 = Fstop2-0.02*BP;            	% Second Passband Frequency

Dstop1 = 1e-005;            % First Stopband Attenuation
Dpass  = 0.005;             % Passband Ripple
Dstop2 = 1e-005;            % Second Stopband Attenuation
dens   = 20;                % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
                          0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
a = 1;
Hd = dfilt.dffir(b);

% [EOF]
