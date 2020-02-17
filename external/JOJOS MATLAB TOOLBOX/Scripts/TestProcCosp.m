% % TEST procCosp
% This script aims at testing procCosp.m
% Here the user can compute multivariate time series composed of a number 
% of sources (sinus) with different frequencies, phase shift, additive noise.

% Last version: 2013
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

close all
clear all


% Signal parameters
Fe = 128;
SNR = 1;
f1 = 10;     % freq source 1
f2 = 15;    % freq source 2
f3 = 36;    % freq source 3

%% Source 1
t_axis = [0:1/Fe:100];
phi = 0;
s1_0   = sin(f1*(2*pi*t_axis) - phi) ;
if(SNR)   s1_0 = s1_0 + (1/SNR)*randn(size(t_axis));    end

phi = pi/12;
s1_15 =sin(f1*(2*pi*t_axis) - phi)+ (1/SNR)*randn(size(t_axis));
if(SNR)   s1_15 = s1_15 + (1/SNR)*randn(size(t_axis));    end
    
phi = pi/6;
s1_30 =sin(f1*(2*pi*t_axis) - phi)+ (1/SNR)*randn(size(t_axis));
if(SNR)   s1_30 = s1_30 + (1/SNR)*randn(size(t_axis));    end
    
phi = pi/4;
s1_45 =sin(f1*(2*pi*t_axis) - phi)+ (1/SNR)*randn(size(t_axis));
if(SNR)   s1_45 = s1_45 + (1/SNR)*randn(size(t_axis));    end
    
phi = pi/2;
s1_90 =sin(f1*(2*pi*t_axis) - phi)+ (1/SNR)*randn(size(t_axis));
if(SNR)   s1_90 = s1_90 + (1/SNR)*randn(size(t_axis));    end
    
phi = pi/18000; % 0.01 deg
s1_0001 =sin(f1*(2*pi*t_axis) - phi)+ (1/SNR)*randn(size(t_axis));
if(SNR)   s1_0001 = s1_0001 + (1/SNR)*randn(size(t_axis));    end
    
figure('Name', ['Source 1 @' num2str(f1) 'Hz with phase 0, 45 and 90°'],'WindowStyle','docked'); 
plot(t_axis,s1_0,'b'),hold on, 
plot(t_axis,s1_45,'r'),hold on, 
plot(t_axis,s1_90,'g'), hold off,
xlim([0 (2*pi)/10])

%% Source 2
phi = 0;
s2_0   = sin(f2*(2*pi*t_axis) - phi) + (1/SNR)*randn(size(t_axis));
if(SNR)   s2_0 = s2_0 + (1/SNR)*randn(size(t_axis));    end

phi = pi/4;
s2_45 =sin(f2*(2*pi*t_axis) - phi)+ (1/SNR)*randn(size(t_axis));
if(SNR)   s2_45 = s2_45 + (1/SNR)*randn(size(t_axis));    end

phi = pi/2;
s2_90 =sin(f2*(2*pi*t_axis) - phi)+ (1/SNR)*randn(size(t_axis));
if(SNR)   s2_90 = s2_90 + (1/SNR)*randn(size(t_axis));    end
% figure('Name', ['Source 2 @' int2str(f2) 'Hz with phase 0, 45 and 90°'],'WindowStyle','docked'); 
% plot(t_axis,s1_0,'b'),hold on, 
% plot(t_axis,s1_45,'r'),hold on, 
% plot(t_axis,s1_90,'g'), hold off,
% xlim([0 (2*pi)/10])

%% Source 3
phi = 0;
s3_0   = sin(f3*(2*pi*t_axis) - phi) + (1/SNR)*randn(size(t_axis));
if(SNR)   s3_0 = s3_0 + (1/SNR)*randn(size(t_axis));    end

phi = pi/4;
s3_45 =sin(f3*(2*pi*t_axis) - phi)+ (1/SNR)*randn(size(t_axis));
if(SNR)   s3_45 = s3_45 + (1/SNR)*randn(size(t_axis));    end

phi = pi/2;
s3_90 =sin(f3*(2*pi*t_axis) - phi)+ (1/SNR)*randn(size(t_axis));
if(SNR)   s3_90 = s3_90 + (1/SNR)*randn(size(t_axis));    end
% figure('Name', ['Source 3 @' int2str(f3) 'Hz with phase 0, 45 and 90°'],'WindowStyle','docked'); 
% plot(t_axis,s3_0,'b'),hold on, 
% plot(t_axis,s3_45,'r'),hold on, 
% plot(t_axis,s3_90,'g'), hold off,
% xlim([0 (2*pi)/10])

%% Process cospectra on a multiple sets of sources
WindowType  = 'hanning' ;   %'hanning' or 'welch' or ...
WindowSize  = 128; 
Overlaping	= 0.75;
PhaseCorr  	= 0;            % Phase shift correction
FreqRange   = [Fe 5 25];        % eg. [Fe 5 25]
% ("sx_y" --> x: source num, y: phase shift in degree)
Dataset1    = [s1_0 ];
Dataset2    = [s1_45];

% Put source s1 in third dataset only during defined event period
event_t1    = fix(30*Fe); % event begin @30s
event_t2    = fix(70*Fe); % event end @70s
s1_event    = s1_0;
s1_event([1:event_t1,event_t2:end]) = 0;
Dataset3    = [s1_event]; 
MASK        = ones(1,length(t_axis));
MASK(event_t1:event_t2) = 0; % mask on second sinus

% Data to analyse is the concatenation of the three sets
DATA    = [Dataset1;Dataset2;Dataset3];
Sdata   = [size(Dataset1,1);size(Dataset2,1);size(Dataset3,1)];

[S] = procCosp(DATA', WindowType, WindowSize, Overlaping, [], MASK, PhaseCorr);
C = real(S);
Q = imag(S);

%% display results
NbFreq = size(C,3);

% frequency axis calculation
if(isempty(FreqRange))
    f_axis = [0:NbFreq-1] * .5*Fe/NbFreq;
else
    Fmin = FreqRange(2);
    Fmax = FreqRange(3);
    f_axis = [0:NbFreq-1] * (Fmax-Fmin)/(NbFreq-1) + Fmin;
end

% spectrum calculation
spectrum = zeros(3,NbFreq);
for f = 1 : NbFreq
   spectrum(:,f)     = (diag(C(:,:,f)));
end

% normalize cospectra (APPROXIMATE correlation, since sources are not independent)
for f = 1 : NbFreq
    C(:,:,f) = diag(diag(C(:,:,f).^(-1/2))) * C(:,:,f) * diag(diag(C(:,:,f)).^(-1/2));
end

% visualization of spectrum and cospectra
ImagescMatrixSelection(C,[1,1,1]);
figure('NumberTitle','off','WindowStyle', 'docked'), 
PlotChannelSelection(spectrum, f_axis, 'frequency', 'power spectrum')




