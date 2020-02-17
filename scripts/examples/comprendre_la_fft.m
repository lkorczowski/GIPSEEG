clear all; 
for i=1:10
    close all;
fs=128
N=100
t=(1:N)/fs;
freq1=10;
freq2=9;
s=sin((2*pi*freq1)*t)+sin((2*pi*freq2)*t)
s=s+0.5*randn(size(s));
figure;subplot(211);stem(t,s)

padding=zeros(1,length(s)*0);
s=[s padding];
N=length(s)
S=fft(s);
freq=0:fs/N:fs-fs/N
subplot(212);stem(freq,log(abs(S).^2)+1)

padding=zeros(1,length(s)*5);
s=[s padding];
N=length(s)
S=fft(s);
freq=0:fs/N:fs-fs/N
hold all 
stem(freq,log(abs(S).^2))
hold off

Save(:,i)=S
end

Sm=mean(Save,2)
freq=0:fs/N:fs-fs/N
hold all 
stem(freq,log(abs(Sm).^2))
hold off

axis([0.5 30 min(log(abs(Sm).^2)) max(log(abs(Sm).^2))])