function [Hd b a] = designBandPassFilter_FIR_strongCut(FO,BP)
%DESIGNBANDPASSFILTER_FIR_STRONGCUT Returns a discrete-time filter object.
%
% When BP is not specified, this filter will always be of order 215 with
% band-pass of 4Hz. Otherwise order will increase as BP decreases.

% *** Inputs *** 
% - FO      --> central frequency, in Hz.
% - BP      --> (optional) Band-pass, in Hz. Default value: 4Hz

% Equiripple Bandpass filter designed using the FIRPM function.

% All frequency values are in Hz.

% default BP
if(nargin < 2)||(isempty(BP))   BP = 4 ;  end


Fs = 128;  % Sampling Frequency


% Fstop1 = 3;               % First Stopband Frequency
% Fpass1 = 5;               % First Passband Frequency
% Fpass2 = 23;              % Second Passband Frequency
% Fstop2 = 25;              % Second Stopband Frequency

Fstop1 = FO-1.25*BP;            	% First Stopband Frequency
Fpass1 = Fstop1+0.3*BP;             % First Passband Frequency
Fstop2 = FO+1.25*BP;               	% Second Stopband Frequency
Fpass2 = Fstop2-0.3*BP;            	% Second Passband Frequency

Fstop1 = clamp(Fstop1,[0.5 Fs/2]); 	% First Stopband Frequency
Fstop2 = clamp(Fstop2,[0.5 Fs/2]); 	% Second Stopband Frequency

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
