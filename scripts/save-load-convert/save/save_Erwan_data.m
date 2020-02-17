clear data
Users=Generate_Users_Numbers(1:24)
BASEDIR='D:\data\erwanRAW'
user='24'
for i=1:length(Users)
data(i)=load_erwan_data(Users{i},BASEDIR,0);

end
%% save with format training
clear training online
for i=1:length(Users)
    i
    for j=1:length(data(i).session)
        Phase=find(cellfun(@(x) strcmp(x,'non-adaptive'), data(i).session{1}.phasetype));
        training(i,j)=data(i).session{j}.phase{Phase}.training;
        online(i,j)=data(i).session{j}.phase{Phase}.online;
    end
end

save('D:\data\erwan\FullErwanData_RAW.mat','training','online')
%% save with format daa

save('D:\data\erwan\FullErwanData_RAW.mat','data')
%save('ErwanData_S1_non-adaptive.mat','training','online')
%%
%% for comparison size with Marco files
k=1
BASEDIR='D:\data\Data-New-Marco-All\'
Names=Generate_Users_Numbers(1:24)
for i=1:24
    USERDIR = fullfile(BASEDIR,Names{i})

    listdir = dir(USERDIR);
Nsession = length(listdir)-2;
     phasetypefile=fullfile('D:\data\Data-New-Marco-All\',Names{i},listdir(k+2).name,'Order.txt')
    formatSpec = '%s';
    [Text]=fopen(phasetypefile)
    phasetype = textscan(Text,formatSpec,...
    'Delimiter', '\n', ...
    'CollectOutput', true);
    phasetype=phasetype{:};
    Phase=find(cellfun(@(x) strcmp(x,'non-adaptive'), phasetype))
    
sEEGtxt(i)=size(EEGtxt{i}.session{1}.phase{Phase}.training,2)
end
sEEG-sEEGtxt

for i=1:24
sEEG(i)=size(training(i).s,1);
end
sEEGtxt-sEEG

%%