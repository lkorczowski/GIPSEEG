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