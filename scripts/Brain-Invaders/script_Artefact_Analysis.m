%% LOAD FILES WITH MAT
clear all;
close all
%%for 1 players from the training data
Fs=128;
window=1*Fs; %window
nb_players=71;
AUsers = Generate_Users_Numbers(1:nb_players);
%Users = {'11','15','18','21'}; % selected players
Directory= 'D:\data\Hyperscanning\BI-multiplayers\Groups\Rejected_parts\'
SAVEPATH=[Directory 'Artefacts.mat'];
load(SAVEPATH)
Directory= 'D:\data\Hyperscanning\BI-multiplayers\Groups\'
load([Directory 'ALLrawGroupsData.mat'])
load([Directory 'Groups.mat'])

clear Results

for i=1:length(Groups)
    Flash{i}=test{i}.Flash;
    Y{i}=test{i}.Flash;
    [X{i} isBad{i}]=epoch_p300(test{i}.s',Flash{i},window,0,Artefacts{i});
    Results{i,1}=ALLgroups{i};
    disp(['Group' 9 9 9 num2str(ALLgroups{i})])
        Results{i,2}=length(Flash{i})/Fs;
    disp(['Temps d enregistrement :' 9 num2str(length(Flash{i})/Fs) ' s'])
    Results{i,3}=length(isBad{i});
    disp(['Number of epoch corrupted:' 9 num2str(length(isBad{i}))])
     Results{i,4}=length(find((isBad{i})));
    disp(['Number of epoch corrupted:' 9 num2str(length(find((isBad{i}))))])
    Results{i,5}=length(find((isBad{i})))/length(isBad{i})*100;
    disp(['Pourcentage of epoch corrupted:' 9 num2str(length(find((isBad{i})))/length(isBad{i})*100) ' %'])
    
end

