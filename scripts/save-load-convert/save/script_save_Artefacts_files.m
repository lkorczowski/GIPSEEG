%% save rejection files
clear all
Directory= 'D:\data\Hyperscanning\BI-multiplayers\Groups\Rejected_parts\'
SAVEPATH=[Directory 'Artefacts.mat'];
        file = fullfile(Directory,'*.txt');
        FILES=dir(file);
        for i=1:length(FILES)
Artefacts{i}=load([Directory FILES(i).name]);
Groups{i}=FILES(i).name;
        end
        
        save(SAVEPATH, 'Artefacts', 'Groups')