%%% (0) PREPARE Load the data player, filter, preprocessing
clear all
Directory='D:\data\Hyperscanning\MARTI\Groups\' %office
% Directory='E:\data\MARTI\' %home

load(['MARTI_GroupsName.mat'])
LABELS=getfield(load(['MARTI_ElectrodesName.mat']),'MARTImapping');
DISP=0;startDISP=1700;endDISP=2000;

for indGroup=1:length(GroupsName)
    disp(['group ' GroupsName{indGroup}])
    try
        load([Directory 'mat\dataPlayers\' GroupsName{indGroup} 'dataPlayers.mat'])
        for indPlayer=1:2
            EEG=dataPlayers(indPlayer);
            EEG.ElectrodesName=LABELS;
            
            
            if DISP, figure;
                Freqs=freqs_val(EEG.Fs,size(EEG.Channels(startDISP*EEG.Fs:endDISP*EEG.Fs,:),1)/2);
                eegf=abs(fft(EEG.Channels(startDISP*EEG.Fs:endDISP*EEG.Fs,:))).^2;
                subplot(211)
                plot(Freqs,log(mean(eegf(1:end/2,:),2)))
            end
            EEG=preprocessingEEG(EEG,[.1 40 4 4]);
            
            if DISP
                eegf=abs(fft(EEG.Channels(startDISP*EEG.Fs:endDISP*EEG.Fs,:))).^2;
                subplot(212)
                plot(freqs_val(EEG.Fs,size(EEG.Channels(startDISP*EEG.Fs:endDISP*EEG.Fs,:),1)/2),log(mean(eegf(1:end/2,:),2)))
                eegplot([EEG.Channels';EEG.Trigger'*100],'srate',128)

                ACSTPoptions.Epoch_size=128;
                ACSTPoptions.LatencyCorr_max=0
                ACSTPoptions.Mask_Electrodes=find(ismember(EEG.ElectrodesName,{'Pz','Cz','O1','O2','Oz'}));
                ACSTPoptions.Mask_Time=6:96;
                [~, ACSTPstruct]=ACSTP(EEG,ACSTPoptions);
            end
            
            %%
            prestim=1; % 3 seconds
            poststim=2;
            Overlap=0;
            EEGFT=EEG_LK2FT(EEG, prestim,poststim,Overlap);
            
            EEGFT.Conditions=EEG.Conditions;
            save([Directory '\mat\epoch\' GroupsName{indGroup} '_' num2str(indPlayer) '.mat'],'EEGFT')
            %             disp('file written')
        end
    catch e
        disp('error')
    end
end

%% When you saved Xte, Yte, ... concatenate that into ALLdata
clear all
Directory='D:\data\Hyperscanning\MARTI\Groups\' %office

load([Directory 'mat\epoch\ALL_epochs.mat'])

ALLdata.Xte=Xte; clear Xte;
ALLdata.Yte=Yte; clear Yte;
ALLdata.dimord='chan_time_rpt';
ALLdata.poststim=poststim;
ALLdata.prestim=prestim;
ALLdata.labels=LABELS;
ALLdata.conditions=Cond;
ALLdata.groups=GroupsName;
save([Directory 'mat\epoch\ALL_epochs.mat'], 'ALLdata')

%% When you concatenated into ALLdata, save ALLdata per condition
clear all
Directory='D:\data\Hyperscanning\MARTI\Groups\' %office

load([Directory 'mat\epoch\ALL_epochs.mat'])
load(['MARTI_GroupsName.mat'])
%check if condition are good

ALLdata.subject=GroupsName;
tmp1=ALLdata.Xte;
tmp2=ALLdata.Yte;
tmp3=ALLdata.conditions;
for indCond=unique(tmp3{1})'
    ALLdata.error=[];
    ALLdata=removefields(ALLdata,{'Xte','Yte','conditions'});
    for indG=1:length(GroupsName)
            try
            ALLdata.Xte{indG,1}=[];
            ALLdata.Xte{indG,1}=cat(1,tmp1{indG,1}(:,:,tmp3{indG,1}==indCond),...
                tmp1{indG,2}(:,:,tmp3{indG,2}==indCond));
            ALLdata.Yte(indG,:)={tmp2{indG,1}(tmp3{indG,1}==indCond),...
                tmp2{indG,2}(tmp3{indG,2}==indCond)};
            ALLdata.condition=indCond;
            catch e
                ALLdata.error=[ALLdata.error indG];
            end
    end
    save([Directory 'mat\epoch\ALL_epochs_' num2str(indCond) '.mat'], 'ALLdata')
end

