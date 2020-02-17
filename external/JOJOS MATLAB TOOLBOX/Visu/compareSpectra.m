function compareSpectra(DATA1,DATA2,WINDOW_SIZE,WINDOW_TYPE,OVERLAP, FS)    
%   function CompareSpectra(data1,data2,window_size,window_type,overlap, Fs)
% 
% This function plot superimposed spectra of two inputs vectors
%
% *** Inputs *** 
% - DATA1           --> input data vector 1
% - DATA2           --> input data vector 2
% ...
%
%
%  *** History *** 
% First version: 12/07/2012
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


% default frequency sampling
if(nargin < 6)||(isempty(FS))                   FS  = true(1,size(DATA,1));  end
% default overlapping
if(nargin < 5)||(isempty(OVERLAP))              OVERLAP = .75;  end
% default windows type
if(nargin < 4)||(isempty(WINDOW_TYPE))          WINDOW_TYPE = 'hanning';  end
% default windows size
if(nargin < 3)||(isempty(WINDOW_SIZE))          WINDOW_SIZE = 1024;  end


switch lower(WINDOW_TYPE)
    case 'hamming'
        win = hamming(WINDOW_SIZE);
    case 'hann'
        win = hann(WINDOW_SIZE);
    case 'hanning'
        win = hanning(WINDOW_SIZE);
    case 'blackman'
        win = blackman(WINDOW_SIZE);
    case 'barthannwin'
        win = barthannwin(WINDOW_SIZE);
    case 'blackmanharris'
        win = blackmanharris(WINDOW_SIZE);
    case 'bohmanwin'
        win = bohmanwin(WINDOW_SIZE);
    case 'chebwin'
        win = chebwin(WINDOW_SIZE);
    case 'gausswin'
        win = gausswin(WINDOW_SIZE);
    case 'kaiser'
        win = kaiser(WINDOW_SIZE);
    case 'nuttallwin'
        win = nuttallwin(WINDOW_SIZE);
    case 'parzenwin'
        win = parzenwin(WINDOW_SIZE);
    case 'tukeywin'
        win = tukeywin(WINDOW_SIZE);
    case 'rectangular'
        win = ones(WINDOW_SIZE,1);
    case 'flattopwin'
        win = flattopwin(WINDOW_SIZE);
    case 'welch'
         %win = welchwin(WINDOW_SIZE, 0);
        Nd2=(WINDOW_SIZE-1)/2;
        for i=1 : WINDOW_SIZE
            win(i)=1-((i-Nd2)/Nd2)^2;
        end
        win = win';
    otherwise
        disp('WARNING - Unknown window type!');
        disp('Switching to default window type : Hamming Window');
        win = hamming(WINDOW_SIZE);
end

[PSD_beforeFilt, f_axis] = pwelch(DATA1,win, OVERLAP*WINDOW_SIZE, 1*WINDOW_SIZE, FS); 
[PSD_afterFilt] = pwelch(DATA2,win, OVERLAP*WINDOW_SIZE, 1*WINDOW_SIZE, FS); 
figure('Name',['Power spectrum before/after filtering'], 'NumberTitle','off','WindowStyle', 'docked'),
plot(f_axis, 10*log10(PSD_beforeFilt), 'b'), hold on,
plot(f_axis, 10*log10(PSD_afterFilt), 'r'),
grid on, xlabel('Frequency (Hz)'), ylabel('Magnitude (dB)'), 
legend('before filt','after filt');