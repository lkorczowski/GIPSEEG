close all; clear all; clc;

cfg.path_input      = 'E:\CPS\BciZar\P300-Marco\Data-Marco';
cfg.path_output     = 'E:\CPS\BciZar\P300-Marco\Data-New-Marco';
cfg.subjects        = 1:2; % 1:7; 8:24
cfg.sessions        = 1:2; % 1:8; 1

for subject = cfg.subjects
    for session = cfg.sessions
        disp(['Subject ' num2str(subject) ', session ' num2str(session)]);
        for phase = 0:1
			% output
			order_output = importdata([cfg.path_output '\' num2str(subject, '%02d') ...
                '\Session' num2str(session) '\Order.txt']);
            outputFilePath = [cfg.path_output '\' num2str(subject, '%02d') ...
                '\Session' num2str(session) ...
                '\Phase' num2str(phase) '\Online'];
            count_output = dlmread([outputFilePath '\Count.txt']);
            % input
			if (strcmp(order_output{phase+1}, 'adaptive'))
				inputFilePath = [cfg.path_input '\' num2str(subject, '%02d') ...
					'\Session' num2str(session) ...
					'\Phase' num2str(1)];
			elseif (strcmp(order_output{phase+1}, 'non-adaptive'))
				inputFilePath = [cfg.path_input '\' num2str(subject, '%02d') ...
					'\Session' num2str(session) ...
					'\Phase' num2str(0)];
            else
                error('Order file...');
			end
            count_input = dlmread([inputFilePath '\Count.txt']);
			% compare...
            if ((count_input(1) ~= count_output(1)) || (count_input(2) ~= count_output(2)))
                error('Diferent number of targets');
            end
            EEG_output = dlmread([outputFilePath '\EEG.txt']);
            STIM_output = dlmread([outputFilePath '\STIM.txt']);
            idx_ta = find(STIM_output(:,2)==1);
            for i = 1:count_input(1)
                EEG_input = dlmread([inputFilePath '\@TA' num2str(i) '.txt']);
                range = STIM_output(idx_ta(count_input(1)-i+1),1);
                [d] = compare(EEG_input, EEG_output(:,range:(range+128-1)));
                if (d)
                    error('Difference...');
                end
            end
            idx_nt = find(STIM_output(:,2)==0);
            for i = 1:count_input(2)
                EEG_input = dlmread([inputFilePath '\@NT' num2str(i) '.txt']);
                range = STIM_output(idx_nt(count_input(2)-i+1),1);
                [d] = compare(EEG_input, EEG_output(:,range:(range+128-1)));
                if (d)
                    error('Difference...');
                end
            end
        end
    end
end
