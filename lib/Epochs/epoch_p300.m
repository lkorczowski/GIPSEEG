function [X isBad] = epoch_p300(Signal,Target,window,offset,Artefacts)
% 1 %%%%%% FOR LK_TOOLBOX %%%%%%%%%
% [X isBad] = epoch_p300(Signal,Target,window,offset,Artefacts)
% compute epochs X from continuous EEG Signals with Target (flashes). The
% epochs will be of the length window. You can add an offset to all the
% epochs. Artefacts is a vector of the length of Target that notes the
% artefacts. Thus isBad returns a 1 when an epoch in X contains an
% artefact.
%
% INPUTS
% Signal: [nb_electrodes x nb_samples], 
% Target: [nb_samples x 1] The trigger channel with  zero everywhere but when the epochs starts.
%        There are [nb trials] non-zero elements
% window: [a scalar] the length in samples of the final epochs
% offset: [nb trials x1] a vector containing the possible lag of the epochs
%         from the trigger channel Target. Default: 0. It can be also a scalar, in
%         this case, the offset will be applied for all the epochs (in case you
%         want to study a second pre-stimulation, put -128 if the sample rate is
%         128 Hz).
% Artefacts: [nb_samples x 1] is a vector of zeros but where the signal is
%         supposed to contain artefacts.
% 
% OUTPUTS
% X: [nb_electrodes x window x nb_trials] the epoched EEG into a 3D matrix
% isBad: [nb_trials x 1] a vector of boolean. 1 if the correponding epoch
%        in X do contain Artefacts (see the corresponding INPUT). 0
%        otherwise.
%
% 2 %%%%% FOR FIELDTRIP %%%%%
%
% [X] = epoch_p300(Signal,trl)
% compute epochs X from continuous EEG Signals with trl.
%  The trial definition "trl" is an Nx3 matrix, N is the number of trials.
%   The first column contains the sample-indices of the begin of each trial 
%   relative to the begin of the raw data, the second column contains the 
%   sample-indices of the end of each trial, and the third column contains 
%   the offset of the trigger with respect to the trial. An offset of 0 
%   means that the first sample of the trial corresponds to the trigger. A 
%   positive offset indicates that the first sample is later than the trigger, 
%   a negative offset indicates that the trial begins before the trigger.
%
% Signal is [nb_electrodes x nb_samples]
%
% *** History: 2015-03-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
% From Alexandre Barachand (c) 2013, modified by Louis Korczowski 2015
%
% see also: unepoch_p300, epoch_EEG_struct

X=[];
isBad=[];
if nargin==2
    % convert to Fieldtrip structure
    for indTrl=1:size(Target,1)
    X{indTrl}=Signal(:,Target(indTrl,1)+Target(indTrl,3) :  Target(indTrl,2)+Target(indTrl,3));    
    end
else
    %extracts the epoch matrices with respect to the initial signal.
    if (nargin<5 || isempty(Artefacts))
        Artefacts=zeros(1,length(Target)); %no artifacts
    end
    if (nargin < 4 || isempty(offset))
        offset = 0; %no offset
    end
    ixTarget = find(Target);
    
    
    while (ixTarget(end)+window-1)>size(Signal,2)
       ixTarget(end) = []; 
    end
    
    NTarget = length(ixTarget);
    
    X = zeros(size(Signal,1),window,NTarget);
    isBad=zeros(1,NTarget);
    for i=1:NTarget
        if (ixTarget(i)+window-1+offset)<1 || (ixTarget(i)+offset)>size(Signal,2) % total out-of-bound epoch : error
            error('epoch_P300: out of bound epoch, check offset and epoch size.')
        elseif (ixTarget(i)+offset)<1 % partial out-of-bound epoch : warning + filling
            warning('epoch_P300: offset makes epochs out of bound (negative sample). Trunked epoch, filling with zeros (check isBad).')
            isBad(i)=1;
            X(:,2-(ixTarget(i)+offset):end,i) = Signal(:,1:ixTarget(i)+window-1+offset);
        elseif (ixTarget(i)+window-1+offset)>size(Signal,2) % partial out-of-bound epoch : warning + filling
            warning('epoch_P300: offset makes epochs out of bound (max reached). Trunked epoch, filling with zeros (check isBad).')
            isBad(i)=1;
            X(:,1:(Size(Signal,2)-(ixTarget(i)+window-1+offset)),i) = Signal(:,ixTarget(i)+offset:end);
        else
        X(:,:,i) = Signal(:,ixTarget(i)+offset:ixTarget(i)+window-1+offset);
        isBad(i)= length(find(Artefacts(ixTarget(i)+offset:ixTarget(i)+window-1+offset)>0));
        end

    end
end