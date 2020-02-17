function EEG=EEG_cnt2mat(file)
% INPUTS :
% ------
% file is the file path of the .cnt and .dat, exemple:
% file='D:\data\ARNAUD\cba\cba1ff01'
% it will import both
% D:\data\ARNAUD\cba\cba1ff01.cnt
% D:\data\ARNAUD\cba\cba1ff01.dat

% OUTPUT :
% ------
% EEG is a structure with
%              Fs: scalar (sample rate in Hz)
%         Trigger: [nb samples x1 ] Trigger channel of '0' with '1' at the start
%                   of each sweep. There are [nb epochs] '1'.
%      EpochClass: [nb epochs x1] class of the sweeps (0 for Non-TARGET, 1
%                   for TARGET).
%       EpochType: [nb epochs x1] specific for the study 
%                   0 = target in the animal categorization task
%                   1 = target in the "easy" animal recognition task
%                   2 = target in the "hard" animal recognition task
%                   3 = target in the non-animal recognition task
%                   4 (not used)
%                   5 = distractor in the animal categorization task
%                   6 = distractor in the "easy" animal recognition task
%                   7 = distractor in the "hard" animal recognition task
%                   8 = distractor in the non-animal recognition task
%    EpochCorrect: [nb epochs x1] specific for the study 0 is incorect (subject
%                   responded on a distractor or failed to respond on a target)
%                   and 1 is correct (subject responded on a target or did not
%                   respond on a distractor).
%      EpochCount: [nb epochs x1] the number of the Epoch given in
%                  EEG.Trigger. If 0, it is not an Epoch.
%        Channels: [nb samples x nb channels] preprocessed EEG recordings
%   NoiseTrigger*: the equivalent of Trigger but for the sweep of the noise.
%                   By default, it takes the same.
%
% *** History: 05-May-2015
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
%%

cnt = loadcnt([file '.cnt']);%, 'dataformat','int32')
EEG.Channels=cnt.data';
EEG.ElectrodesLocation=cnt.electloc;
for indElec=1:length(EEG.ElectrodesLocation)
EEG.ElectrodesName{indElec}=EEG.ElectrodesLocation(indElec).lab;
end
for i=1:length(cnt.event)
    EEG.EpochCount(i)=cnt.event(i).stimtype;
    EEG.Trigger(i)=cnt.event(i).offset;
    EEG.Fs=cnt.header.rate;
end
EEG.Trigger=EEG.Trigger(find(EEG.EpochCount)); %remove 0
EEG.EpochCount=EEG.EpochCount(find(EEG.EpochCount)); %remove 0
fid = fopen([file '.dat'],'r');
    i = 1;
    countStim=1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
        if (i>20) && ~(strcmp(num2str(A{i}),'-1')) %get stimulation code and so
            tmp=str2num(A{i});
            TrialNumber(countStim)=tmp(1);
            TrialClass(countStim)=tmp(2);
            TrialType(countStim)=floor((tmp(3)-12000)*1e-4);
            TrialCorrect(countStim)=tmp(4);
            countStim=countStim+1;
        end
    end
    fclose(fid);
    EEG.EpochCorrect=TrialCorrect;
    EEG.EpochType=TrialType;
    EEG.EpochClass=TrialClass;
    
