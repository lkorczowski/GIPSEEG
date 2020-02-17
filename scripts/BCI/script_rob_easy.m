%% load data
close all
clear all
tic
for PREP=0:1
    load(['example_rob.mat'])

DISP=0

if PREP
eeg=reshape(eegP300,[size(eegP300,1),size(eegP300,2)*size(eegP300,3)])';
Fs=500
if DISP, figure; 

Freqs=freqs_val(Fs,size(eeg,1)/2);
eegf=abs(fft(eeg)).^2;
subplot(211)
plot(Freqs,log(mean(eegf(1:end/2,:),2)))
end
[eeg Fs]=preprocessingEEG(eeg,Fs,[1 20 2 4]);

if DISP
eegf=abs(fft(eeg)).^2;
subplot(212)
plot(freqs_val(Fs,size(eeg,1)/2),log(mean(eegf(1:end/2,:),2)))
end
trigger=zeros(size(eeg,1),1);
trigger(1:Fs:end)=1;
eegP300=epoch_p300(eeg',trigger,Fs);

end

%% 
COMB{1}=1:22; %electrode combinaison
for ind=1:100
% compute everythings
X=eegP300;
Y=labels;
P=[0.5]; % pourcentage of test data (i.e 50% against 50% of training)

[Xtr Ytr Xte Yte RND]=Generate_Training_Test_Set(X(COMB{1},:,:),Y',P,[]);

[COVtr, P1] = covariances_p300(Xtr,Ytr);

COVte = covariances_p300(Xte,P1);

% Classification with minimum to mean distance (MDM)
[Yestimated Dist C COVte] = mdm(COVte,COVtr,Ytr);

Accuracy=sum(Yte'==Yestimated')/length(Yestimated);
[PerfX,PerfY,~,AUC(ind,PREP+1),OPTROCPT] = perfcurve(Yte,-diff(Dist'),1);
if DISP
ConfM=confusionmat(Yte,Yestimated)
disp(['Accuracy: ' num2str(Accuracy) '%'])
figure
plotroc(Yte,-diff(Dist'))
text(0.5,0.5,['AUC:' num2str(AUC) ' (' num2str(Accuracy) '%)' ])
end
end
disp('DONE')
end
toc

%%
mean(AUC)