clear all; close all
Fs=200;
f1=0.7
f2=40
N=2
NoiseVar=0.1
Casual=1
w1=f1/(Fs/2); % LOW cutoff frequency
w2=f2/(Fs/2); % HIGHT cutoff frequency
%         [b,a]=ellip(N,3,20,[w1 w2]);
[b,a]=butter(N,[w1 w2]);

abs( pole((tf(b,a))))
ex_time=(-Fs:Fs)';%(2+1/Fs) seconds of signal


% ------------- PLOT1 ------------
%simulated low freq wave
ex_signal=zeros(size(ex_time));
ex_signal(ex_time==0)=10;%generate dirac of amplitude 10µV

w= window(@hamming,Fs/2) %generate very low freq (0.5 Hz)
ex_signal=conv(ex_signal,w,'same');
ex_signal=smooth(ex_signal,Fs/4);
ex_signal=ex_signal+randn(size(ex_signal))*NoiseVar;
subplot(521)
plot(ex_time,ex_signal)
if Casual
    ex_signal_filt = filter(b,a,ex_signal);
else
    ex_signal_filt = filtfilt(b,a,ex_signal);
end
hold all
plot(ex_time,ex_signal_filt)
hold off

ylabel('1 Hz')

subplot(522)
[Pxx]= fft(ex_signal-mean(ex_signal),Fs*4);
Pxx=sqrt(Pxx.*conj(Pxx));
Pxx=Pxx/length(Pxx);
freqs=0:Fs/length(Pxx):Fs-Fs/length(Pxx);
semilogx(freqs(1:length(Pxx)/2),Pxx(1:length(Pxx)/2))
axis([0 Fs 0.9*min(Pxx) 1.1*max(Pxx)])

[Pxx]= fft(ex_signal_filt-mean(ex_signal_filt),Fs*4);
Pxx=sqrt(Pxx.*conj(Pxx));
Pxx=Pxx/length(Pxx);
freqs=0:Fs/length(Pxx):Fs-Fs/length(Pxx);
hold all
semilogx(freqs(1:length(Pxx)/2),Pxx(1:length(Pxx)/2))
hold off
legend('pre',['post ' num2str(f1) '-' num2str(f2)])
% ------------- PLOT2 ------------
%simulated high freq wave
ex_signal=zeros(size(ex_time));
ex_signal(ex_time==0)=10;

w= window(@hamming,Fs/10) %generate low freq (5 Hz)
ex_signal=conv(ex_signal,w,'same');
ex_signal=smooth(ex_signal,Fs/20);
ex_signal=ex_signal+randn(size(ex_signal))*NoiseVar;

subplot(523)
plot(ex_time,ex_signal)
% axis([-Fs/2 Fs/2 -5 10])
if Casual
    ex_signal_filt = filter(b,a,ex_signal);
else
    ex_signal_filt = filtfilt(b,a,ex_signal);
end
hold all
plot(ex_time,ex_signal_filt)
hold off
ylabel('5 Hz')

subplot(524)
[Pxx]= fft(ex_signal-mean(ex_signal),Fs*4);
Pxx=sqrt(Pxx.*conj(Pxx));
Pxx=Pxx/length(Pxx);
freqs=0:Fs/length(Pxx):Fs-Fs/length(Pxx);
semilogx(freqs(1:length(Pxx)/2),Pxx(1:length(Pxx)/2))
axis([0 Fs 0.9*min(Pxx) 1.1*max(Pxx)])

[Pxx]= fft(ex_signal_filt-mean(ex_signal_filt),Fs*4);
Pxx=sqrt(Pxx.*conj(Pxx));
Pxx=Pxx/length(Pxx);
freqs=0:Fs/length(Pxx):Fs-Fs/length(Pxx);
hold all
semilogx(freqs(1:length(Pxx)/2),Pxx(1:length(Pxx)/2))
hold off
% ------------- PLOT3 ------------
%simulated high freq wave

ex_signal=zeros(size(ex_time));
ex_signal(ex_time==0)=10;

w= window(@hamming,Fs/20) %generate low freq (15 Hz)
ex_signal=conv(ex_signal,w,'same');
ex_signal=smooth(ex_signal,Fs/40);
ex_signal=ex_signal+randn(size(ex_signal))*NoiseVar;

subplot(525)
plot(ex_time,ex_signal)

if Casual
    ex_signal_filt = filter(b,a,ex_signal);
else
    ex_signal_filt = filtfilt(b,a,ex_signal);
end
hold all
plot(ex_time,ex_signal_filt)
hold off
ylabel('10 Hz')

subplot(526)
[Pxx]= fft(ex_signal-mean(ex_signal),Fs*4);
Pxx=sqrt(Pxx.*conj(Pxx));
Pxx=Pxx/length(Pxx);
freqs=0:Fs/length(Pxx):Fs-Fs/length(Pxx);
semilogx(freqs(1:length(Pxx)/2),Pxx(1:length(Pxx)/2))
axis([0 Fs 0.9*min(Pxx) 1.1*max(Pxx)])

[Pxx]= fft(ex_signal_filt-mean(ex_signal_filt),Fs*4);
Pxx=sqrt(Pxx.*conj(Pxx));
Pxx=Pxx/length(Pxx);
freqs=0:Fs/length(Pxx):Fs-Fs/length(Pxx);
hold all
semilogx(freqs(1:length(Pxx)/2),Pxx(1:length(Pxx)/2))
hold off

% ------------- PLOT4 ------------
%simulated high freq wave

ex_signal=zeros(size(ex_time));
ex_signal(ex_time==0)=10;

w= window(@hamming,Fs/40) %generate low freq (30 Hz)
ex_signal=conv(ex_signal,w,'same');
ex_signal=smooth(ex_signal,Fs/80);
ex_signal=ex_signal+randn(size(ex_signal))*NoiseVar;

subplot(527)
plot(ex_time,ex_signal)

if Casual
    ex_signal_filt = filter(b,a,ex_signal);
else
    ex_signal_filt = filtfilt(b,a,ex_signal);
end
hold all
plot(ex_time,ex_signal_filt)
hold off

ylabel('20Hz')

subplot(528)
[Pxx]= fft(ex_signal-mean(ex_signal),Fs*4);
Pxx=sqrt(Pxx.*conj(Pxx));
Pxx=Pxx/length(Pxx);
freqs=0:Fs/length(Pxx):Fs-Fs/length(Pxx);
semilogx(freqs(1:length(Pxx)/2),Pxx(1:length(Pxx)/2))
axis([0 Fs 0.9*min(Pxx) 1.1*max(Pxx)])

[Pxx]= fft(ex_signal_filt-mean(ex_signal_filt),Fs*4);
Pxx=sqrt(Pxx.*conj(Pxx));
Pxx=Pxx/length(Pxx);
freqs=0:Fs/length(Pxx):Fs-Fs/length(Pxx);
hold all
semilogx(freqs(1:length(Pxx)/2),Pxx(1:length(Pxx)/2))
hold off

% ------------- PLOT5 ------------
%simulated high freq wave

ex_signal=zeros(size(ex_time));
ex_signal(ex_time==0)=10;

w= window(@hamming,Fs/60) %generate low freq (60 Hz)
ex_signal=conv(ex_signal,w,'same');
ex_signal=smooth(ex_signal,Fs/120);
ex_signal=ex_signal+randn(size(ex_signal))*NoiseVar;

subplot(529)
plot(ex_time,ex_signal)

if Casual
    ex_signal_filt = filter(b,a,ex_signal);
else
    ex_signal_filt = filtfilt(b,a,ex_signal);
end
hold all
plot(ex_time,ex_signal_filt)
hold off

ylabel('30 Hz')
xlabel('time (s)')

subplot(5,2,10)
[Pxx]= fft(ex_signal-mean(ex_signal),Fs*4);
Pxx=sqrt(Pxx.*conj(Pxx));
Pxx=Pxx/length(Pxx);
freqs=0:Fs/length(Pxx):Fs-Fs/length(Pxx);
semilogx(freqs(1:length(Pxx)/2),Pxx(1:length(Pxx)/2))
axis([0 Fs 0.9*min(Pxx) 1.1*max(Pxx)])

[Pxx]= fft(ex_signal_filt-mean(ex_signal_filt),Fs*4);
Pxx=sqrt(Pxx.*conj(Pxx));
Pxx=Pxx/length(Pxx);
freqs=0:Fs/length(Pxx):Fs-Fs/length(Pxx);
hold all
semilogx(freqs(1:length(Pxx)/2),Pxx(1:length(Pxx)/2))
xlabel('freqs (Hz)')
hold off

% ------------- ANNOTATION and SAVE ------------
set(gcf,'color',[1 1 1])
print(gcf,['.\FiltResponse\lo' num2str(f1) '_hi' ...
    num2str(f2) '_ord' num2str(N) 'zeroDist.tif'],'-dtiff','-r450')