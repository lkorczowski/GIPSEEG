function [X Freq]=Spectra(x,fs)
X=abs(fft(x))/size(x,1);
X=X(1:length(x)/2,:);
if nargin >1
    Freq=0:fs/length(x):fs-1/length(x);

else
    Freq=1:length(x);
end
    Freq=Freq(1:length(x)/2);
