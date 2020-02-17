% function load and manipulation ERWAN data


%% convert the matlab old structure to the new for ACSTP
directory='F:\data\erwan\mat\'
load([directory 'FullErwanData_RAW.mat'])
load('BImapping')

for indSubject=1:length(data)
    for indSession=1:length(data(indSubject).session);
        for indPhase=1:length(data(indSubject).session{indSession}.phasetype)
            
            EEG=convert_MATgdt_2_MAT3eg(data(indSubject).session{indSession}.phase{indPhase}.training);
            EEG.ElectrodesName=LABELS;
            EEG.SessionType=data(1).session{1}.phasetype{indPhase}
            save([directory 'ERWAN_SS' num2str(indSubject) '_s' num2str(indSession) '_training_' EEG.SessionType],'EEG')
            
               EEG=convert_MATgdt_2_MAT3eg(data(indSubject).session{indSession}.phase{indPhase}.online);
            EEG.ElectrodesName=LABELS;
            EEG.SessionType=data(1).session{1}.phasetype{indPhase}
            save([directory 'ERWAN_SS' num2str(indSubject) '_s' num2str(indSession) '_online_' EEG.SessionType],'EEG')
            
            
            
                end
    end
end