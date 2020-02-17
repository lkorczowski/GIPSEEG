close all
Fs=128;
timeS=0:1/Fs:1-1/Fs;
subjects=[04,06,08,09];
figure
compt=1;
for i=subjects
subplot(length(subjects),1,compt)
PLOTs(:,:,compt)=mean(P1all(7,:,i)',2);
plot(timeS,PLOTs(:,:,compt))
compt=compt+1;
end

figure
plot(-1+1/Fs:1/Fs:1-1/Fs,xcorr(PLOTs(:,:,1),PLOTs(:,:,4),'unbiased'))
%plot(-1+1/Fs:1/Fs:1-1/Fs,xcorr(PLOTs(:,:,1),5*randn(128,1)))