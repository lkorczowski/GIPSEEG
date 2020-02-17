%% load gdf TRAINING (1player) EEG data and save it as matlab file
clear all
for i=1:71
    if i<10
            AUsers{i}=['0' num2str(i)];
    else
    AUsers{i}=num2str(i);
    end
end
Directory='D:\data\Hyperscanning\BI-multiplayers\Training\';
SAVEPATH=[Directory 'ALLsolo.mat'];
if ~exist(SAVEPATH,'file')
SaveEEGdata(Directory,AUsers,SAVEPATH)
disp('saved')
else
    disp('file already exist')
end
%% load and check the sample rate of the data
clear all
Directory='D:\data\Hyperscanning\BI-multiplayers\Training\';
SAVEPATH=[Directory 'ALLsolo.mat'];
load(SAVEPATH)

%%get sample freq
for i=1:length(ALLdata.Xte)
S(i)=size(ALLdata.Xte{i},2);
end
S~=128

%% manual correction if needed + save
i=42 %subject
se=1 %session
if i<10
 data.subjects{i} = load_EEG_data(Directory,['0' num2str(i)]);
else
     data.subjects{i} = load_EEG_data(Directory,[num2str(i)]);
end
    % Epoch signal
    test = data.subjects{i}.session{se}.online;
    Fs = test.Fs;
    window = Fs; % 1s window
    ALLdata.Xte{i} = epoch_p300(test.s,test.Flash,window);
    ALLdata.Yte{i}=test.Y;
    
    save('D:\data\Hyperscanning\BI-multiplayers\Training\ALLsolo.mat','ALLdata')
    
%% check the data for NaN and size incoherence
for i=1:length(ALLdata.Xte)
NbSamples(i)=size(ALLdata.Xte{i},2); % check the epochs' sizes
SizeDifference(i)=-size(ALLdata.Xte{i},3)+size(ALLdata.Yte{i},1); % check the number of epochs
P1(:,:,i)=mean(ALLdata.Xte{i}(:,:,ALLdata.Yte{i}(1:size(ALLdata.Xte{i},3))==1),3);
NaNs(i)=length(find(isnan(ALLdata.Xte{i})));
end

Nb_errors=length(find(NbSamples~=128))+length(find(SizeDifference))
indNaN=find(NaNs)

%% correct size incoherence + save
for i=1:length(ALLdata.Xte)
ALLdata.Yte{i}=ALLdata.Yte{i}(1:size(ALLdata.Xte{i},3));
end
for i=1:length(ALLdata.Xte)
NbSamples(i)=size(ALLdata.Xte{i},2); % check the epochs' sizes
SizeDifference(i)=-size(ALLdata.Xte{i},3)+size(ALLdata.Yte{i},1); % check the number of epochs
P1(:,:,i)=mean(ALLdata.Xte{i}(:,:,ALLdata.Yte{i}(1:size(ALLdata.Xte{i},3))==1),3);
NaNs(i)=length(find(isnan(ALLdata.Xte{i})));
end


Nb_errors=length(find(NbSamples~=128))+length(find(SizeDifference))
indNaN=find(NaNs)

    save(SAVEPATH,'ALLdata')

%% load gdf GROUPS (2player) EEG data and save it as matlab file
clear all

Directory='D:\data\Hyperscanning\BI-multiplayers\Groups\';
SAVEPATH=[Directory 'ALLgroups.mat'];
load([Directory 'Groups.mat']) % load ALLgroups
if ~exist(SAVEPATH,'file')
SaveEEGepoch2p(Directory,ALLgroups,SAVEPATH)
disp('saved')
else
    disp('file already exist')
end