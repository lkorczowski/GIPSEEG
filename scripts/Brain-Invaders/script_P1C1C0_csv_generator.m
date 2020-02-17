%% create C0 C1 for reference with all-stats LOUIS
clear all;
close all
Users = {'27','34'};
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