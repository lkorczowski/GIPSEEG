function [Hd b a] = designLowPassFilter_FIR_strongCut(Fpass,Fs)
%DESIGNLOWPASSFILTER_FIR_STRONGCUT Returns a discrete-time filter object.

%
% M-File generated by MATLAB(R) 7.9 and the Signal Processing Toolbox 6.12.
%
% Generated on: 17-Apr-2013 17:38:56
%

% Equiripple Lowpass filter designed using the FIRPM function.


% All frequency values are in Hz.
% default Fs
if(nargin < 3)||(isempty(Fs))   Fs = 128 ;  end

Fstop = Fpass+1;               % Stopband Frequency
Dpass = 0.057501127785;  % Passband Ripple
Dstop = 1e-005;          % Stopband Attenuation
dens  = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
a = 1;
Hd = dfilt.dffir(b);

% [EOF]