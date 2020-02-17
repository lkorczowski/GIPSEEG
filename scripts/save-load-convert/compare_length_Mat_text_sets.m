for i=1:24
sEEG2(i)=size(online(i).s,1)
end
for i=1:24
sEEG2(i)=size(training(i).s,1)
end
for i=1:24
sEEGtxt(i)=size(EEGtxt{i}.session{1}.phase{1}.training.s,2)
end
sEEGtxt(i)=size(EEGtxt{i}.session{1}.phase{1}.training,2)
%%
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