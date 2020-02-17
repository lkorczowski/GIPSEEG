nbStims=[400 80]

NbStim2Take={round(nbStims*0.5) round(nbStims*0.25) round(nbStims*0.2)...
    round(nbStims*0.15) round(nbStims*0.1) round(nbStims*0.05)}

for indUser=13:24%:length(users) %select the subject
    for indSession=1%:8
        close all
        clear EEG
        nameSave='_variability_';
        Directory= ['D:\data\erwan\mat\']  ; %change path if needed
        for indBootstrap=1:length(NbStim2Take)
            clear FILTER R1 R2 INFOS
            outputFILE1=[Directory 'results\' 'ERWAN_' nameSave '_s' num2str(indUser) '_se' num2str(indSession) '_boot'  num2str(indBootstrap) '.mat'];
            R1=load(outputFILE1)
            outputFILE2=[Directory 'results\' '2ERWAN_' nameSave '_s' num2str(indUser) '_se' num2str(indSession) '_boot'  num2str(indBootstrap) '.mat'];
            R2=load(outputFILE2)
            for indR=1:length(R1.FILTER)
                R1.FILTER(indR).PzSNR=[];
            end
            FILTER=[R1.FILTER, R2.FILTER]
            INFOS=R1.INFOS
            INFOS.RNDseed= [R1.INFOS.RNDseed,R2.INFOS.RNDseed];
            INFOS.NbTrials=R1.INFOS.NbTrials+R2.INFOS.NbTrials;
            delete(outputFILE1,outputFILE2)
            disp(['old files ' '_s' num2str(indUser) '_boot'  num2str(indBootstrap) ' deleted'])
            save(outputFILE1,'FILTER','INFOS')
            disp(['file ' '_s' num2str(indUser) '_boot'  num2str(indBootstrap) ' saved'])
        end
        
        
    end
end