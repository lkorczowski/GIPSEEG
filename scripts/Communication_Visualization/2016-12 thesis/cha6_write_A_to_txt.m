    %% EXPORT A (for LORETA)
    FOLDER='D:\GoogleDrive\Présentations & Communications\2016-12 Thesis EEG-Hyperscanning\thesis\@LORETA\'
    TYPES={'AJD','BAJD','CAJD'};
    for tix=1:24
            for indT=1:length(TYPES)
    A=inv(results(tix,eix,six).(TYPES{indT}).B);
    dlmwrite([ FOLDER 's' num2str(tix) '_' TYPES{indT} '.txt'],A,'delimiter','\t')
            end
    end
    