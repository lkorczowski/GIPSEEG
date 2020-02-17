%%
clear all
close all
data=load_EEG_data('E:\EEG Screening BI march-april 2015\','S10')
EEG=data.session{1}.online.s;
Tag=data.session{1}.online.Y;
Flash=data.session{1}.online.Flash;
Target=data.session{1}.online.Target;
Fs=data.session{1}.online.Fs;

IND=(find(Flash==1))
length(find(Target==1))
for i=1:length(IND)
epoch(:,i)=(IND(i):IND(i)+127);
epochEEG(:,:,i)=EEG(epoch(:,i),:);
end

MeanTarget=mean(epochEEG(:,:,Tag==1),3);


MeanNonTarget=mean(epochEEG(:,:,Tag==0),3);

  figure
  subplot(211)
  
plot(MeanTarget)
xlabel('Target')
axis([0 128 -10 10])
subplot(212)
plot(MeanNonTarget)
axis([0 128 -10 10])

xlabel('Non Target')
%find target
%take 1s after every target
%mean of all epoch

eegplot(MeanTarget')

figure
semilogy(abs(fft(MeanTarget(:,7))))