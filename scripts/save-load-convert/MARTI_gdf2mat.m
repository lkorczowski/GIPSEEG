%%
clear all
close all
Users=Generate_Users_Numbers(1:50)
Directory='D:\data\Hyperscanning\MARTI\Screening\'
%% for saving every subjects and sessions into 1 matlab file
for indU=24:length(Users)
Screening(indU)=load_EEG_data(Directory,Users{indU})
end
%pure raw data !!!!
save([Directory 'RAW_screening_data.mat'],'Screening')
%% for saving every session into 1 matlab file for each subject
for indU=1:length(Users)
subject=load_EEG_data(Directory,Users{indU})
save([Directory 'RAW_screening_data_' Users{indU} '.mat'],'subject')
end
%pure raw data !!!!