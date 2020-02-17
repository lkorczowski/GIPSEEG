close all; clear all; clc;

cfg.path_input      = 'E:\CPS\BciZar\P300-Marco\Data-Old-Marco';
cfg.path_output     = 'E:\CPS\BciZar\P300-Marco\Data-New-Marco';
cfg.subjects        = 8:24; % 1:7; 8:24
cfg.sessions        = 1; % 1:8; 1
    
for subject = cfg.subjects
    for session = cfg.sessions
        disp(['Subject ' num2str(subject) ', session ' num2str(session)]);
        PHASE_FLAG = nan(2,1);
        for phase = 1:2
            inputFilePath = [cfg.path_input '\' num2str(subject, '%02d') ...
                '\Session' num2str(session) ...
                '\Phase' num2str(phase)];
            PHASE_TYPE = importdata([inputFilePath '\Results.csv-adaptive']);
            outputFilePath = [cfg.path_output '\' num2str(subject, '%02d') ...
                '\Session' num2str(session) ...
                '\Phase' num2str(phase-1)];
            if (strcmp(PHASE_TYPE,'False'))
                PHASE_FLAG(phase) = 0;
            elseif (strcmp(PHASE_TYPE,'True'))
                PHASE_FLAG(phase) = 1;
            else
                error('Unknown phase type');
            end
            % Training - read
            TRAIN_EEG = importdata([inputFilePath '\training.csv']);
            TRAIN_STM = importdata([inputFilePath '\training.csv-Stimulations']);
            % Online - read
            ONLINE_EEG = importdata([inputFilePath '\online.csv']);
            ONLINE_STM = importdata([inputFilePath '\online.csv-Stimulations']);
            % Training - write
            mkdir([outputFilePath '\Training\']);
            [EEG, STM] = parse(TRAIN_EEG, TRAIN_STM);
                writeMatrix([outputFilePath '\Training\EEG.txt'], EEG);
                writeStim([outputFilePath '\Training\STIM.txt'], STM);
                
            nTA = length(find(STM(:,2)==1)); nNT = length(find(STM(:,2)==0));
            writeMatrix([outputFilePath '\Training\Count.txt'], [nTA; nNT], '%0.0f ');
            % Online - write
            mkdir([outputFilePath '\Online\']);
            [EEG, STM] = parse(ONLINE_EEG, ONLINE_STM);
                writeMatrix([outputFilePath '\Online\EEG.txt'], EEG);
                writeStim([outputFilePath '\Online\STIM.txt'], STM);
            
            nTA = length(find(STM(:,2)==1)); nNT = length(find(STM(:,2)==0));
            writeMatrix([outputFilePath '\Online\Count.txt'], [nTA; nNT], '%0.0f ');
            if (phase == 2)
                % Write order
                text = cell(1,2);
                for i = 1:2
                    if (PHASE_FLAG(i))
                        text{i} = 'adaptive';
                    else
                        text{i} = 'non-adaptive';
                    end
                end
                outputFilePath = [cfg.path_output '\' num2str(subject, '%02d') ...
                    '\Session' num2str(session)];
                writeAscii([outputFilePath '\Order.txt'], text);
            end
        end
    end
end
