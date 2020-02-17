function [EEG RND]=schuffle_continuousEEG(EEG,NbStim2Take,RND)
% schuffle_continuousEEG will remove from EEG an certain amount of
% stimulation code in order to keep only NbStim2Take. The result Stim2Take
% is a boolean vector that say if the Stimulation code is included or not.
%
% [EEG Stim2Take]=schuffle_continuousEEG(EEG,Stim2Take)
%    Generate a new set of data from EEG with the epoch selected in
%    Stim2Take.
%
% [EEG Stim2Take]=schuffle_continuousEEG(EEG,Stim2Take,Classes,NbStim2Take)
%    Generate a randomly selected set of data from EEG. The function will
%    keep only NbStim2Take for every Classes.

%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%
% EEG is a structure with
%         Trigger: [nb samples x1 ] Trigger channel of '0' with '1' at the start
%                   of each sweep. There are [nb epochs] '1'.
%      EpochClass: [nb epochs x1] class of the sweeps (0 for Non-TARGET, 1
%                   for TARGET).
%
% RND*:            [nb epochs x1] a randompermutation seed
%
% NbStim2Take:     [nb unique stim x1] the number of stimulation to keep
%                  for each class. WARNING: put thoses number in the
%                  ascending order of the StimCode for the classes
%                  contained in EEG.EpochClass.
%                  [scalar] In that case, it takes equaly the same amount
%                  of stimulation.
%
%%%%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%%%
% EEG, the modified continuous eeg structure
% RND, the random permutation seed




if nargin<3 || isempty(RND)
    RND=randperm(length(EEG.EpochClass));
end

if length(NbStim2Take)==1
    NbStim2Take=ones(1,length(unique(EEG.EpochClass)))*NbStim2Take;
end

RND2=sort(schuffle_indices_classwise(RND,EEG.EpochClass,NbStim2Take));


EEG.EpochClass=EEG.EpochClass(RND2); %take only need StimCode

tmp=EEG.Trigger;
if length(tmp)==length(EEG.EpochClass) % EEG.Trigger is StimPositions
    EEG.Trigger=tmp(RND2);
else % EEG.Trigger is a stimulation channel
    lentmp=length(tmp);
    tmp=find(tmp);
    tmp=tmp(RND2);
    EEG.Trigger=Position2Trigger(tmp,lentmp);
end