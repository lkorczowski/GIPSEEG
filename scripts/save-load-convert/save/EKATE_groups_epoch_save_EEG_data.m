
clear all
Directory='D:\data\Hyperscanning\EKATE\Groups\';
load([Directory 'Groups.mat'])
    %SAVEPATH=[Directory 'Full_EEG.mat'];
AUsers=ALLgroups
 nusers=length(AUsers);

% P300 moyen et structure de bruit
Xtr={};
    Xte={};
    Yte={};
    data={};
    %load data set
    se=1; %selection session for users data 

      errors=[];  
for i = 1:length(AUsers)
    % load data
    %try
    data.subjects{i} = load_EEG_data(Directory,AUsers{i});
    for se=3
    % Epoch signal
    %selection session for users data 
    %    se=1; % player1
    %            se=2; % player2
    %se=3; % all data

    test = data.subjects{i}.session(se);
    Fs = test.Fs;
    window = Fs*1.5; % 1.5s window
    offset = -Fs*0.5; % -0.5s of offset
    tic
    Xte{i} = epoch_p300(test.Channels',test.Trigger,window,offset);
    toc
    Yte{i}=test.EpochClass;
%     EEGgroup(i)=test;
%     EEG1(i)=test;
%     EEG2(i)=test;
%     Flash{i,se}=test.Flash;
    end
    %catch e
     %   errors=[errors; AUsers{i}];
    %end
end

user=AUsers{i}
    ALLdata.Xte=Xte;
    ALLdata.Yte=Yte;
    save([Directory 'mat\epoch\ALL_epochs.mat'],'ALLdata') % save epochs (for immediate classification)
%     save([Directory 'Flash.mat'],'Flash') % save tags (for online BI behaviour analysis)
    %save([Directory 'Full_EEG.mat'],'EEGgroup','EEG1','EEG2') %save all data to mat