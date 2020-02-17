%% loading data and computing general P1
clear all;
close all
Users = {'06'}
Users = {'01','02','03','04','05','06','07','08','09','10',...
    '11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
% P300 moyen et structure de bruit
figure
for i = 1:length(Users)
    % load data
    data = load_erwan_data(Users{i},'D:\data\erwan');
    % Epoch signal
    
   
    Xtr=[];
    Y=[];
    for se=1:length(data.session)
        for ph=1:length(data.session{se}.phase)
            training = data.session{se}.phase{ph}.training;
            Fs = training.Fs;
            window = Fs; % 1s window
            Xtr =cat(3,Xtr, epoch_p300(training.s,training.Flash,window));
            Y=cat(1,Y,training.Y);
            
        end
    end
    
             
    P1(:,:,i) = mean(Xtr(:,:,Y==1),3);
    [temp I]=max(diag(P1(:,:,i)*P1(:,:,i)'));
    NoiseStructure(:,:,i) = cov(training.s(20*Fs:end,:));
    subplot(length(Users),1,i);plot(P1(I,:,i)); % show greatest P1 for each patient
end



%% Riemanian distance P1
C=[];
C(:,:,1)=eye(16);
C=cat(3,C,covariances(P1));
Dr=[];
for i=1:size(C,3)
    for j=1:size(C,3)
Dr(i,j)=distance_riemann(C(:,:,i),C(:,:,j));
    end
end

%% save distance
Riem_dist=Dr;
save('d:\data\summary.mat','Riem_dist','P1','NoiseStructure')

%% find 4 closest
p=(10e-1<Dr)
d=(Dr<9)
p==d
Close=(10e-1<Dr)&&(Dr<9);
%% Generation of synthetic data from real P300

clear all;
Users = {'15','16','18','21'};
Xtr=[];
% mean P300 and noise structure
for i = 1:length(Users)
    % load data
    data = load_erwan_data(Users{i});
    training = data.session{1}.phase{1}.training;
    test = data.session{1}.phase{1}.online;
    % Epoch signal
    Fs = training.Fs; 
    window = Fs; % 1s window
    Xtr(:,:,:,i) = epoch_p300(training.s,training.Flash,window);

    P1(:,:,i) = mean(Xtr(:,:,training.Y==1,i),3);
    P0(:,:,i) = mean(Xtr(:,:,training.Y==0,i),3);
    NoiseStructure(:,:,i) = cov(training.s(20*Fs:end,:));
    disp(['User ' Users{i} ' loaded'])
end
%% generate training signal
[b,a] = butter(5,[1 20]/(Fs/2));
NoiseLevel = 1;
trSignal = [];
for i = 1:length(Users)
   
    s = conv2(training.Target,P1(:,:,i)');
    s = s(1:length(training.Target),:);
   
    % generate noise
    noise = mvnrnd(zeros(size(P1,1),1),NoiseLevel*NoiseStructure(:,:,i),length(training.Target));
    noise = filtfilt(b,a,noise);
    
    s = s + noise;
   
    trSignal = cat(2,trSignal,s);
end

% generate test signal

[b,a] = butter(5,[1 20]/(Fs/2));
NoiseLevel = 1;
teSignal = [];
for i = 1:length(Users)
   
    s = conv2(test.Target,P1(:,:,i)');
    s = s(1:length(test.Target),:);
   
    % generate noise
    noise = mvnrnd(zeros(size(P1,1),1),NoiseLevel*NoiseStructure(:,:,i),length(test.Target));
    noise = filtfilt(b,a,noise);
    
    s = s + noise;
   
    teSignal = cat(2,teSignal,s);
end

%% test classification hyperscanning

% Epoch signal
Fs = training.Fs; 
window = Fs; % 1s window
Xtr = epoch_p300(trSignal,training.Flash,window);
Xte = epoch_p300(teSignal,test.Flash,window);

% Estimation des matrices de covariance spéciale P300
[COVtr, P1] = covariances_p300_hyper(Xtr,training.Y,4);
COVte = covariances_p300(Xte,P1);

%% Classification par MDM riemannienne
A=1
[Yte, distances] = mdm(COVte,COVtr,training.Y,{'',''}); 

% matrice de confusion
disp('matrice de confusion');
disp(confusionmat(test.Y,Yte));

% courbe roc
scores = -diff(distances')';
plotroc(test.Y',scores');

% Area under curve (higer is better, max=1);
auc = area_under_curve(scores',test.Y');
disp(['Area Under Curve : ' num2str(auc)]);

%% create C0 C1 for reference (problem non positivity)
clear all;
close all
load('d:\data\summary.mat','P1')
P=mean(P1,3);
clear P1
Users = {'15','16','18','21'};
%Users = {'01','02','03','04','05','06','07','08','09','10',...
 %   '11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
% P300 moyen et structure de bruit
figure
number=[];
for i = 1:length(Users)
    % load data
    data = load_erwan_data(Users{i});
    % Epoch signal
    
   
    Xtr=[];
    Xte=[];
    Y=[];
    Yt=[];
    for se=1:length(data.session)
        for ph=1:length(data.session{se}.phase)
            training = data.session{se}.phase{ph}.training;
            test = data.session{se}.phase{ph}.online;
            Fs = training.Fs;
            window = Fs; % 1s window
            Xtr =cat(3,Xtr, epoch_p300(training.s,training.Flash,window));
            Xte = cat(3,Xte,epoch_p300(test.s,test.Flash,window));
            Y=cat(1,Y,training.Y);
            Yt=cat(1,Yt,test.Y);
        end
    end
    
             
    P1(:,:,i) = mean(Xtr(:,:,Y==1),3);
    P0(:,:,i) = mean(Xtr(:,:,Y==0),3);
    [temp I]=max(diag(P1(:,:,i)*P1(:,:,i)'));
    NoiseStructure(:,:,i) = cov(training.s(20*Fs:end,:));
    subplot(length(Users),1,i);plot(P1(I,:,i)); % show greatest P1 for each patient
end
X1=unepoch_P300(P1);
X1=[P;X1];
X0=unepoch_P300(P0);
X0=[mean(P0,3);X0];

C1=cov(X1');
C0=cov(X0');

distance_riemann(C1,C0)

%% create C0 C1 for reference ERWAN
clear all;
close all
Users = {'15','16','18','21'};
%Users= {'15'};
%Users = {'01','02','03','04','05','06','07','08','09','10',...
 %   '11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
% P300 moyen et structure de bruit
figure
number=[];
for i = 1:length(Users)
    % load data
    data = load_erwan_data(Users{i});
    % Epoch signal
    
   
    Xtr=[];
    Xte=[];
    Y=[];
    Yt=[];
    for se=1:1%length(data.session)
        for ph=1:1%length(data.session{se}.phase)
            training = data.session{se}.phase{ph}.training;
            test = data.session{se}.phase{ph}.online;
            Fs = training.Fs;
            window = Fs; % 1s window
            Xtr =cat(3,Xtr, epoch_p300(training.s,training.Flash,window));
            Xte = cat(3,Xte,epoch_p300(test.s,test.Flash,window));
            Y=cat(1,Y,training.Y);
            Yt=cat(1,Yt,test.Y);
        end
    end
    % Estimation des matrices de covariance spéciale P300
[COVtr, P1(:,:,i)] = covariances_p300(Xtr,training.Y);
COVte = covariances_p300(Xte,P1(:,:,i));

% Classification par MDM riemannienne
[Yte, distances,C{i}] = mdm(COVte,COVtr,training.Y); 

% matrice de confusion
disp('matrice de confusion');
disp(confusionmat(Yt,Yte));

% courbe roc
scores = -diff(distances')';
%plotroc(test.Y',scores');

% Area under curve (higer is better, max=1);
auc = area_under_curve(scores',Yt');
disp(['Area Under Curve : ' num2str(auc)]);
             
    [temp I]=max(diag(P1(:,:,i)*P1(:,:,i)'));
    subplot(length(Users),1,i);plot(P1(I,:,i)); % show greatest P1 for each patient
end
%% create C0 C1 for reference without inter-stats LOUIS
clear all;
close all
Users = {'12','49'};
%Users= {'15'};
%Users = {'01','02','03','04','05','06','07','08','09','10',...
 %   '11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
% P300 moyen et structure de bruit
Directory='D:\data\Hyperscanning\BI-multiplayers\Training'
figure
number=[];
for i = 1:length(Users)
    % load data
    data = load_EEG_data(Directory,Users{i});
    % Epoch signal
    
   
    Xtr=[];
    Xte=[];
    Y=[];
    Yt=[];
    for se=1:1%length(data.session)
        for ph=1:1%length(data.session{se}.phase)
            training = data.session{se}.online;
            test = data.session{se}.online;
            Fs = training.Fs;
            window = Fs; % 1s window
            Xtr =cat(3,Xtr, epoch_p300(training.s,training.Flash,window));
            Xte = cat(3,Xte,epoch_p300(test.s,test.Flash,window));
            Y=cat(1,Y,training.Y);
            Yt=cat(1,Yt,test.Y);
        end
    end
    % Estimation des matrices de covariance spéciale P300
[COVtr, P1(:,:,i)] = covariances_p300(Xtr,training.Y);
COVte = covariances_p300(Xte,P1(:,:,i));

% Classification par MDM riemannienne
[Yte, distances,C{i}] = mdm(COVte,COVtr,training.Y); 

% matrice de confusion
disp('matrice de confusion');
disp(confusionmat(Yt,Yte));

% courbe roc
scores = -diff(distances')';
%plotroc(test.Y',scores');

% Area under curve (higer is better, max=1);
auc = area_under_curve(scores',Yt');
disp(['Area Under Curve : ' num2str(auc)]);
             
    %[temp I]=max(diag(P1(:,:,i)*P1(:,:,i)'));
    subplot(length(Users),1,i);plot(P1(:,:,i)'); % show greatest P1 for each patient
end
P1mean=zeros(size(P1(:,:,1)));
for i = 1:length(Users)
    P1mean=P1mean+P1(:,:,i);
end
P1mean=P1mean./length(Users);
figure;plot(P1mean')

covP1=cov(P1mean')

%% build C0 C1 without inter-stats
%REP4p=['D:\gipsa-lab-extensions\scenarios\share\MDM\P300\brain-invaders-mdm-adaptive-2p32\config\' cat(2,Users{:})];
%P1csv=['mdm-generic-P1'];
%P=csvread([REP4p P1csv '.csv'],C1);
nbp=2;
C0=zeros((nbp+1)*16);
C1=zeros((nbp+1)*16);
C1(1:16,1:16)=covP1;
C0(1:16,1:16)=covP1;
%C0(1:16,1:16)=C{1,1}{1,1}(1:16,1:16);
%C1(1:16,1:16)=C{1,1}{2,1}(1:16,1:16);
for i=1:nbp
    ind=(i)*16+1:(i+1)*16;
    C0(ind,ind)=C{1,i}{1,1}(17:end,17:end);
    C1(ind,ind)=C{1,i}{2,1}(17:end,17:end);
end

distance_riemann(C0,C1)
REP4p=['D:\gipsa-lab-extensions\scenarios\share\MDM\P300\brain-invaders-mdm-adaptive-2p32\config\' cat(2,Users{:}) '\'];
%REP4p='D:\data\'
if ~exist(REP4p)
    mkdir(REP4p)
    disp('Folder created')
end
C1csv=['mdm-INITIAL-C1-' int2str(nbp) 'p'];
C0csv=['mdm-INITIAL-C0-' int2str(nbp) 'p'];
P1csv=['mdm-INITIAL-P1-' int2str(nbp) 'p'];
csvwrite([REP4p C1csv '.csv'],C1);
csvwrite([REP4p C0csv '.csv'],C0);
csvwrite([REP4p P1csv '.csv'],P1mean);
%% create C0 C1 for reference with all-stats LOUIS
clear all;
close all
Users = {'01','19'};
%Users= {'15'};
%Users = {'01','02','03','04','05','06','07','08','09','10',...
 %   '11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
% P300 moyen et structure de bruit
Directory='D:\data\Hyperscanning\BI-multiplayers\Training'
%figure
number=[];
nusers=length(Users)
    Xtr=[];
    Xte=[];
    Ytr=[];
    Yte=[];
for i = 1:length(Users)
    % load data
    data = load_EEG_data(Directory,Users{i});
    % Epoch signal
    for se=1:1%length(data.session)
            training = data.session{se}.online;
            test = data.session{se}.online;
            Fs = training.Fs;
            window = Fs; % 1s window
            Xtr{i} =epoch_p300(training.s,training.Flash,window);
            Xte{i} = epoch_p300(test.s,test.Flash,window);
            Ytr{i}=training.Y;
            Yte{i}=test.Y;
    end
end

[Xout, Yout]=Multisubject_generator(Xtr,Ytr);
[Xtraining, Ytraining]=Randomize_trials(Xout, Yout);
disp('training done')
[Xout, Yout]=Multisubject_generator(Xte,Yte);
[Xtest, Ytest]=Randomize_trials(Xout, Yout);

[COVtr, P1] = covariances_p300_hyper(Xtraining,Ytraining,nusers);
COVte = covariances_p300(Xtest,P1);

cond(COVte(:,:,1))

%%% Classification par MDM riemannienne
[Yestimated, distances,C] = mdm(COVte,COVtr,Ytraining); 

% courbe roc
scores = -diff(distances')';
%plotroc(test.Y',scores');

% Area under curve (higer is better, max=1);
auc = area_under_curve(scores',Ytest');
disp(['Area Under Curve : ' num2str(auc)]);
disp(confusionmat(Ytest,Yestimated));

nbp=2;
C0=C{1};
C1=C{2};
REP4p=['D:\gipsa-lab-extensions\scenarios\share\MDM\P300\brain-invaders-mdm-adaptive-2p32\config\' cat(2,Users{:}) '\'];
%REP4p='D:\data\'
if ~exist(REP4p)
    mkdir(REP4p)
    disp('Folder created')
end
C1csv=['mdm-INITIAL-C1-' int2str(nbp) 'p'];
C0csv=['mdm-INITIAL-C0-' int2str(nbp) 'p'];
P1csv=['mdm-INITIAL-P1-' int2str(nbp) 'p'];
csvwrite([REP4p C1csv '.csv'],C1);
csvwrite([REP4p C0csv '.csv'],C0);
csvwrite([REP4p P1csv '.csv'],P1);
disp('SAVE DONE')
%% build N players synchroneous EEG (training and test)

clear all;
close all
%Users = {'15','16','18','21'}; % selected players
Users = {'01','02'}; % selected players

se=2; %selection session for users data
ph=1; %selection phase for users data
%Users= {'15'};
%Users = {'01','02','03','04','05','06','07','08','09','10',...
 %   '11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
% P300 moyen et structure de bruit
r=4; %decimation factor
Xtr={};
    Xte={};
    Y={};
    Yte={};
    %load data set
    
for i = 1:length(Users)
    % load data
    data = load_erwan_data(Users{i});
    % Epoch signal
    training = data.session{se}.phase{ph}.training;
    test = data.session{se}.phase{ph}.online;
    Fs = test.Fs;
    window = Fs; % 1s window
        Xtr{i} = epoch_p300(training.s,training.Flash,window);
    Xte{i} = epoch_p300(test.s,test.Flash,window);
    Yte{i}=test.Y;
    Ytr{i}=training.Y;
end
%% build N players synchroneous EEG (training and test) and save the matrix of performances
%% LOAD FILES WITH GDF
clear all;
close all
%Users = {'07','11','14','15','16','18','21'}; % selected players
Users = {'01','02','03','04'}; % selected players
se=1; %selection session for users data
ph=2; %selection phase for users data
%Users= {'15'};
%Users = {'01','02','03','04','05','06','07','08','09','10',...
 %   '11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
% P300 moyen et structure de bruit
Xtr={};
    Xte={};
    Ytr={};
    Yte={};
    %load data set
    
for i = 1:length(Users)
    % load data
    data = load_erwan_data(Users{i});
    % Epoch signal
    training = data.session{se}.phase{ph}.training;
    test = data.session{se}.phase{ph}.online;
    Fs = test.Fs;
    window = Fs; % 1s window
        Xtr{i} = epoch_p300(training.s,training.Flash,window);
    Xte{i} = epoch_p300(test.s,test.Flash,window);
    Yte{i}=test.Y;
    Ytr{i}=training.Y;
end

[Xout, Yout]=Multisubject_generator(Xtr,Ytr);
[Xtraining, Ytraining]=Randomize_trials(Xout, Yout);

[Xout, Yout]=Multisubject_generator(Xte,Yte);
[Xtest, Ytest]=Randomize_trials(Xout, Yout);
%% classification
REGU=0;
P300_ref_orientation='commonP1'
% Estimation des matrices de covariance spéciale P300
[COVtr, P1] = covariances_p300_hyper(Xtraining,Ytraining,length(Users),P300_ref_orientation);
COVte = covariances_p300_hyper(Xtest,P1,length(Users),P300_ref_orientation);
%condiWoutRegu=cond(COVte(:,:,1))
%condiWRegu=cond(COVte(:,:,1)+10e-5*eye(size(COVte,1)))
%%% Classification par MDM riemannienne
tic
if REGU
[Yestimated, distances,C] = mdm(COVte+repmat(10e-5*eye(size(COVte,1)),[1 1 size(COVte,3)]),COVtr+repmat(10e-5*eye(size(COVtr,1)),[1 1 size(COVtr,3)]),Ytraining); 
else
[Yestimated, distances,C] = mdm(COVte,COVtr,Ytraining); 
end
toc
%+repmat(eye(size(COVte,1)),[1 1 size(COVte,3)])
Perfclassif=length(find(Yestimated==Ytest))/length(Ytest);
disp(['Performance classification using closest reference :' num2str(Perfclassif) '%']);
% matrice de confusion
disp('matrice de confusion');
disp(confusionmat(Ytest,Yestimated));
% courbe roc 
scores = -diff(distances')';
 [TPR,FPR,TH] =roc(Ytest',scores');
%figure;plotroc(Ytest',scores');

% Area under curve (higer is better, max=1);
[X,Y,~,AUC,OPTROCPT] = perfcurve(Ytest',scores',1);
AUC
h=figure;plot(X,Y);
       xlabel('False positive rate'); ylabel('True positive rate')
       title('ROC for classification by logistic regression')
%auc = area_under_curve(scores',Ytest')
%disp(['Area Under Curve : ' num2str(auc)]);
%% LOAD FILES WITH MAT
clear all;
close all
AUsers = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
Users = {'11','15','18','21'}; % selected players

load('D:\data\ALLdata.mat')
RNDtr={};RNDte={};
%
for PCA=0:1;
switch PCA
    case 0
    Xtr=ALLdata.Xtr;
    Xte=ALLdata.Xte;
    case 1
        
      Xtr=epochPCA(ALLdata.Xtr);
      Xte=epochPCA(ALLdata.Xte);
      
end
    Ytr=ALLdata.Ytr;
    Yte=ALLdata.Yte;
 nusers=length(Users);
 allcombi=nchoosek(AUsers,nusers);
auc=[];
se=1; %selection session for users data 
ph=2; %selection phase for users data

 allcombi=nchoosek(AUsers,nusers);
%% generate multiplayers set
IND=find(ismember(AUsers,Users))
%%generate multiplayers set
[Xout, Yout,RNDtr]=Multisubject_generator(Xtr(IND),Ytr(IND),RNDtr);
[Xtraining, Ytraining]=Randomize_trials(Xout, Yout);
disp('training done')
[Xout, Yout, RNDte]=Multisubject_generator(Xte(IND),Yte(IND), RNDte);
[Xtest, Ytest]=Randomize_trials(Xout, Yout);
% end
%% classification
REGU=0;
P300_ref_orientation='commonP1'
% Estimation des matrices de covariance spéciale P300
[COVtr, P1] = covariances_p300_hyper(Xtraining,Ytraining,length(Users),P300_ref_orientation);
COVte = covariances_p300_hyper(Xtest,P1,length(Users),P300_ref_orientation);
%condiWoutRegu=cond(COVte(:,:,1))
%condiWRegu=cond(COVte(:,:,1)+10e-5*eye(size(COVte,1)))
%%% Classification par MDM riemannienne
tic
if REGU
[Yestimated, distances,C] = mdm(COVte+repmat(10e-5*eye(size(COVte,1)),[1 1 size(COVte,3)]),COVtr+repmat(10e-5*eye(size(COVtr,1)),[1 1 size(COVtr,3)]),Ytraining); 
else
[Yestimated, distances,C] = mdm(COVte,COVtr,Ytraining); 
end
toc
%+repmat(eye(size(COVte,1)),[1 1 size(COVte,3)])
Perfclassif=length(find(Yestimated==Ytest))/length(Ytest);
disp(['Performance classification using closest reference :' num2str(Perfclassif) '%']);
% matrice de confusion
disp('matrice de confusion');
disp(confusionmat(Ytest,Yestimated));
% courbe roc 
scores = -diff(distances')';
 [TPR,FPR,TH] =roc(Ytest',scores');
%figure;plotroc(Ytest',scores');

% Area under curve (higer is better, max=1);
[X,Y,~,AUC,OPTROCPT] = perfcurve(Ytest',scores',1);
AUC
figure;plot(X,Y);
       xlabel('False positive rate'); ylabel('True positive rate')
       title('ROC for classification by logistic regression')
%auc = area_under_curve(scores',Ytest')
%disp(['Area Under Curve : ' num2str(auc)]);
end

%% build N players synchroneous EEG (training and test) and save the matrix of performances (ROC area)
clear all;
close all
AUsers = {'01','02','03','04','05'}; % selected players
%AUsers = {'07','11','14'}; % selected players
%AUsers = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
 nusers=2;
 allcombi=nchoosek(AUsers,nusers);
auc=[];
se=1; %selection session for users data 
ph=2; %selection phase for users data

 allcombi=nchoosek(AUsers,nusers);
 load('D:\data\ALLdata.mat')
    Xtr=ALLdata.Xtr;
    Xte=ALLdata.Xte;
    Ytr=ALLdata.Ytr;
    Yte=ALLdata.Yte;
    
for Col=1:1%size(allcombi,1)
Users=allcombi(Col,:);
IND=find(ismember(AUsers,allcombi(Col,:)))
%%
RNDtr={};RNDte={};
%
for PCA=0;%:1
switch PCA
    case 0
    Xtr=ALLdata.Xtr;
    Xte=ALLdata.Xte;
    case 1
        
      Xtr=epochPCA(ALLdata.Xtr);
      Xte=epochPCA(ALLdata.Xte);
      
end
    Ytr=ALLdata.Ytr;
    Yte=ALLdata.Yte;
%%generate multiplayers set
[Xout, Yout,RNDtr]=Multisubject_generator(Xtr(IND),Ytr(IND),RNDtr);
[Xtraining, Ytraining]=Randomize_trials(Xout, Yout);
disp('training done')
[Xout, Yout, RNDte]=Multisubject_generator(Xte(IND),Yte(IND), RNDte);
[Xtest, Ytest]=Randomize_trials(Xout, Yout);
% end
%% classification
for CASE=0%:1
% Estimation des matrices de covariance spéciale P300
[COVtr, P1] = covariances_p300_hyper(Xtraining,Ytraining,nusers);
if CASE, COVtr=covariances_p300_clean_inter(COVtr,nusers); end
COVte = covariances_p300(Xtest,P1);
if CASE, COVte=covariances_p300_clean_inter(COVte,nusers); end

cond(COVte(:,:,1))

%%% Classification par MDM riemannienne
[Yestimated, distances,C] = mdm(COVte,COVtr,Ytraining); 

% matrice de confusion
disp('matrice de confusion');
disp(confusionmat(Ytest,Yestimated));

% courbe roc
scores = -diff(distances')';
[X,Y,~,AUC(Col,CASE+1,PCA+1),OPTROCPT] = perfcurve(Ytest',scores',1);
disp(['Area Under Curve : ' num2str(AUC(Col,CASE+1,PCA+1))]);
figure;f1=plotroc(Ytest',scores');
%legend('random',cat(2,allcombi{Col,:}))
legend('random',char(allcombi(Col,:)))
text(0.5,0.5,['AUC=' num2str(AUC(Col,CASE+1,PCA+1),4)],'FontSize',12)
%saveas(f1,['D:\data\figure\PCA' int2str(PCA) '\figure' cat(2,allcombi{Col,:}) 'inter' int2str(CASE) '.jpeg'])
%close
% Area under curve (higer is better, max=1);
%auc(Col) = area_under_curve(scores',Ytest')
%disp(['Area Under Curve : ' num2str(auc)]);

% h=figure;plot(X,Y);
%        xlabel('False positive rate'); ylabel('True positive rate')
%        title('ROC for classification by logistic regression')
%end
%end
end
end
end
results.AUC=AUC;
results.allcombi=allcombi;
%save('ALLresults.mat','results')


%% build N players synchroneous EEG (training and test) and save the matrix of performances (ROC area)
clear all;
close all
AUsers = {'01','02','03','04'}; % selected players
%AUsers = {'07','11','14'}; % selected players
%AUsers = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
 nusers=length(AUsers);
 allcombi=nchoosek(AUsers,nusers);
auc=[];
se=1; %selection session for users data 
ph=2; %selection phase for users data

 allcombi=nchoosek(AUsers,nusers);
 load('D:\data\ALLdata.mat')
    Xtr=ALLdata.Xtr;
    Xte=ALLdata.Xte;
    Ytr=ALLdata.Ytr;
    Yte=ALLdata.Yte;
    
for Col=1:size(allcombi,1)
Users=allcombi(Col,:);
IND=find(ismember(AUsers,allcombi(Col,:)))
%%
%%generate multiplayers set

[Xout, Yout]=Multisubject_generator(Xtr(IND),Ytr(IND),RNDtr1);
[Xtraining, Ytraining]=Randomize_trials(Xout, Yout);

[Xout, Yout]=Multisubject_generator(Xte(IND),Yte(IND), RNDte1);
[Xtest, Ytest]=Randomize_trials(Xout, Yout);
%% classification
for CASE=0:1
% Estimation des matrices de covariance spéciale P300
[COVtr, P1] = covariances_p300_hyper(Xtraining,Ytraining,nusers);
if CASE, COVtr=covariances_p300_clean_inter(COVtr,nusers); end
COVte = covariances_p300(Xtest,P1);
if CASE, COVte=covariances_p300_clean_inter(COVte,nusers); end

cond(COVte(:,:,1))

%%% Classification par MDM riemannienne
[Yestimated, distances,C] = mdm(COVte,COVtr,Ytraining); 

% matrice de confusion
disp('matrice de confusion');
disp(confusionmat(Ytest,Yestimated));

% courbe roc
scores = -diff(distances')';
[X,Y,~,AUC(Col),OPTROCPT] = perfcurve(Ytest',scores',1);
disp(['Area Under Curve : ' num2str(AUC(Col))]);
figure;f1=plotroc(Ytest',scores');
%legend('random',cat(2,allcombi{Col,:}))
legend('random',char(allcombi(Col,:)))
text(0.5,0.5,['AUC=' num2str(AUC(Col),4)],'FontSize',12)
%saveas(f1,['D:\data\figure\figure' cat(2,allcombi{Col,:}) '.jpeg'])
%close
% Area under curve (higer is better, max=1);
%auc(Col) = area_under_curve(scores',Ytest')
%disp(['Area Under Curve : ' num2str(auc)]);

% h=figure;plot(X,Y);
%        xlabel('False positive rate'); ylabel('True positive rate')
%        title('ROC for classification by logistic regression')
%end
%end
end
end
results.AUC=AUC;
results.allcombi=allcombi;
save('ALLresults.mat','results')