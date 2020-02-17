% save EEG in TEXT
clear all
Directory='D:\data\Hyperscanning\BI-multiplayers\Groups\';
load([Directory 'ALLrawGroupsData.mat']) % load ALLgroups
load([Directory 'Groups.mat']) % load ALLgroups
SAVEPATH=[Directory 'texts_files\'];
AUsers=ALLgroups

%%
for i=1:length(ALLgroups)
A = [test{i}.s test{i}.Flash];

dlmwrite([SAVEPATH ALLgroups{i} '.txt'],A, ' ')
end

%%
Directory='D:\data\erwan\'
SAVEPATH=[Directory 'Weig Arit mean.txt'];
 m_P=textread(SAVEPATH) 
 SAVEPATH=[Directory 'Arit mean.txt'];
 m_Pw=textread(SAVEPATH)
 SAVEPATH=[Directory 'Weig Filt Arit mean.txt'];
 m_Xh=textread(SAVEPATH) 
 SAVEPATH=[Directory 'Weights.txt'];
 m_W=textread(SAVEPATH) 