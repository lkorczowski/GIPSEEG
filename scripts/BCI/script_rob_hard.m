%% load data
close all
clear all
load(['example_rob.mat'])

%% prepare to bruteforce best electrodes (useless when combinaison choosen)
COMB{1}=1:10
%check first electrodes
for indBrute=2:10
COMB{indBrute}=[COMB{indBrute-1} max(COMB{indBrute-1})+1];
end
%check last electrodes
COMB{11}=13:22
for indBrute=12:20
COMB{indBrute}=[COMB{indBrute-1} min(COMB{indBrute-1})-1];
end

%other combinaisons
COMB{21}=[1:5 18:22];
COMB{22}=[1:5 17:22];
COMB{23}=[1:5 16:22];
COMB{24}=[1:5 15:22];
COMB{25}=[1:5 14:22];
COMB{26}=[3:5 14:22];
COMB{27}=[3:5 13:22];
COMB{28}=[8 10 22 19];
COMB{29}=[1:22];
COMB{30}=[1:10 16:22];

%% compute everythings
for indElec=10 % did bruteforce, best combinaison: 10
for indR=26:50
    indR
DISPLAY=0;
X=eegP300;
Y=labels;
P=[0.25];

[Xtr Ytr Xte Yte RND]=Generate_Training_Test_Set(X(COMB{indElec},:,:),Y',P,[]);

if DISPLAY
figure
subplot(131)
plotEEG(mean(X(:,:,Y==0),3),5,128);
subplot(132)
plotEEG(mean(X(:,:,Y==1),3),5,128);
subplot(133)
plotEEG(mean(X(:,:,Y==1)-X(:,:,Y==0),3),5,128);
end

[COVtr, P1] = covariances_p300(Xtr,Ytr);
%eegplot([P1(:,:,1);P1(:,:,2)],'srate',128,'winlength',1)
%return
COVte = covariances_p300(Xte,P1);

%%% Classification par MDM riemannienne
[Yestimated Dist C COVte] = mdm(COVte,COVtr,Ytr);
Accuracy=sum(Yte'==Yestimated')/length(Yestimated);
                [PerfX,PerfY,~,AUC,OPTROCPT] = perfcurve(Yte,-diff(Dist'),1);
if DISPLAY
ConfM=confusionmat(Yte,Yestimated)
disp(['Accuracy: ' num2str(Accuracy) '%'])
figure
                plotroc(Yte,-diff(Dist'))
text(0.5,0.5,['AUC:' num2str(AUC) ' (' num2str(Accuracy) '%)' ])
end
outAUC(indR,1)=AUC; %WARNING indElec is selected only if doing bruteforce, if not put 1 instead
outPERF(indR,1)=Accuracy; %WARNING indElec is selected only if doing bruteforce, if not put 1 instead
end
end
%% plotting for bruteforce electrode selection (only if used)
subplot(411)
plot(mean(outAUC,1))
[tmp,ind]=max(mean(outAUC,1));hold on;plot(ind,tmp,'or')
ylabel('AUC')
subplot(412)
plot(std(outAUC))
[tmp,ind]=min(std(outAUC,1));hold on;plot(ind,tmp,'or')
ylabel('std AUC')
subplot(413)
plot(mean(outPERF,1))
[tmp,ind]=max(mean(outPERF,1));hold on;plot(ind,tmp,'or')
ylabel('%')
subplot(414)
plot(std(outPERF))
[tmp,ind]=min(std(outPERF,1));hold on;plot(ind,tmp,'or')
ylabel('std %')
%% show results perf

mean(outPERF(1:25,1))
std(outPERF(1:25,1))
mean(outPERF(26:50,1))
std(outPERF(25:50,1))

mean(outPERF(1:25,1))
std(outPERF(1:25,1))
mean(outPERF(26:50,1))
std(outPERF(25:50,1))