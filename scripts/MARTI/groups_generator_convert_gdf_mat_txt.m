% This script cumulate the several small script to convert, preprocessing
% the data from the experiment MARTI. All the steps are not mandatory,
% please read carefull each line and put the input and output folders with
% particular care.
%
% see also: load_EEG_data, EEG_mat2txt, preprocessingEEG,
% Multisubject_MARTI_triggers, Multisubject_MARTI_convert_multi2solo

%% (1) Generate Groups Names and save it (usefull only one time but important for the loading)
outputfolder='.\MARTI_GroupsName.mat'
GroupsName={'G010732','G021314','G031136','G042248','G051718','G061647','G070841','G084049',...
    'G096870','G102138','G111045','G122531','G136169','G140346','G151534','G160643','G175164',...
    'G180204','G192350','G202829','G211242'};
save(outputfolder,'GroupsName')

%% (2) Load gdf files and save it to .mat
% WARNING the function load_EEG_data has been changed (new version, structures names changed)
outputfolder='D:\data\Hyperscanning\MARTI\Groups\MARTI_GroupsName.mat'
load(outputfolder) %load GroupsName
for indUser=1:length(GroupsName)
    Directory='D:\data\Hyperscanning\MARTI\Groups\'
    data = load_EEG_data(Directory,GroupsName{indUser})
    save([Directory GroupsName{indUser} '.mat'],'data') % one file per group. No filtering
end

%%  (3) Load .mat files and save it to .txt.  Steps  needed : (2)
% WARNING the functions needed:
% preprocessingEEG (last version)
% EEG_mat2txt (last version)
% WARNING : THIS FUNCTION DOESN'T RETRIEVE THE CONDITIONS (COOP/COMP
% SYNC/NOSYNC)
clear all

f1=1 %low cut-off freq
f2=30 %high cut-off freq
N=4 % filter order (the final order will be N^2)
decimationN=4 %

outputfolder='D:\data\Hyperscanning\MARTI\Groups\MARTI_GroupsName.mat'
load(outputfolder) %load GroupsName

for indUser=1:length(GroupsName)
    Directory='D:\data\Hyperscanning\MARTI\Groups\'
    load([Directory GroupsName{indUser} '.mat'],'data')
    
    for indSession=1:length(data.session)
        EEG=data.session(indSession);
        [EEG.Signal EEG.Fs EEG.Trigger]=preprocessingEEG(EEG.Signal,EEG.Fs,[f1 f2 N decimationN],EEG.Trigger);
        
        EEG.offset=0;
        FULLDIR=['D:\Mes Documents GIPSA\MATLAB\LK_TOOLBOX\data\MARTI\TXT\' GroupsName{indUser} '\Session' num2str(indSession) '\']
        EEG_mat2txt(FULLDIR,EEG)
    end
end
%eegplot(EEG.Signal')

%% (4) load .txt and save artefacts .mat (the artefact should be a .txt or .art file with 0 when sample is good, and something~=0 when sample is bad
outputfolder='D:\data\Hyperscanning\MARTI\Groups\MARTI_GroupsName.mat'
load(outputfolder) %load GroupsName

for indUser=1:length(GroupsName)
    Arte=[];
    for Session=1:4
        try
            Directory='D:\Mes Documents GIPSA\MATLAB\LK_TOOLBOX\data\MARTI\TXT\';% put here your input folder of the artefacts (.txt or .art)
            Arte=[Arte ; load([Directory GroupsName{indUser} '\Session' num2str(Session) '\EEG.art'])]; %put here the input name of the artefacts (.txt or .art)
        end
    end
    Artefacts{indUser}=Arte;
end

save('D:\data\Hyperscanning\MARTI\Groups\artefacts.mat','Artefacts') %save the artefact to mat
%% (5) concatenate all the sessions, filter and save. Steps  needed : (2)
%
%
% WARNING the functions needed:
% preprocessingEEG (last version)
% EEG_mat2txt (last version)
% Multisubject_MARTI_triggers (last version)
% Multisubject_MARTI_convert_multi2solo (IMPROPER, will be changed soon)
%
% output EEG structure format for each Player
%            Fs: scalar
%        StimPos: [nb_epochs*nb_sessions x1 double]
%       StimCode: [nb_epochs*nb_sessions x1 double]
%        Trigger: [nb_samples*nb_sessions x1 double]
%         Target: [nb_samples*nb_sessions x1 double]
%     EpochClass: [nb_epochs*nb_sessions x1 double]
%         Signal: [nb_samples*nb_sessions x nb_Channels]
%              h: [1x nb_sessions struct]
%           name: char
%     Conditions: [nb_epochs*nb_sessions x1 double]
%   Stimulations: [nb_samples*nb_sessions x1 double]

clear all
f1=1
f2=20
N=4
decimationN=4
load('D:\data\Hyperscanning\MARTI\Groups\MARTI_GroupsName.mat')
    Directory='D:\data\Hyperscanning\MARTI\Groups\'
pp=2;%nb players
for indUser=9%:length(GroupsName)
    load([Directory GroupsName{indUser} '.mat'],'data')
    EEG=struct;
    EEG.Trigger=[];
    EEG.Signal=[];
    EEG.StimPos=[];
    EEG.StimCode=[];
    EEG.EpochClass=[];
    
    for indSession=1:length(data.session)
        %concatenate the signals, all !!!!
        
        s=data.session(indSession).Channels;
        EEG.Fs=data.session(indSession).Fs;
                [EEG.EpochClass, EEG.Trigger]= Multisubject_MARTI_triggers(s(:,65),pp); % trigextr3(s(:,65),1,1);  %(3) Sweeps' indices   % don't change for TA2
                [EEG.Channels EEG.Fs EEG.Trigger]=preprocessingEEG(s(:,1:64),EEG.Fs,[f1 f2 N decimationN],EEG.Trigger);
                

        
    end
    %CONCATENATE THE SESSIONS HERE
    dataPlayers=Multisubject_MARTI_convert_multi2solo(data,pp) % WARNING NOT FINISHED (possibly improper function)
    
    % SAVE HERE FOR ONE FILE PER GROUP
        save([Directory GroupsName{indUser} '_allsessions.mat'],'dataPlayers')

end
    % SAVE HERE FOR ONE FILE FOR ALL

