function [Xte Yte isBad Flash]=epoch_EEG_struct(EEG,window,offset,Artefacts)
%  [Xte Yte isBad Flash]=epoch_EEG_struct(EEG,window,offset,Artefacts)
% see also: epoch_p300, unepoch_p300
if nargin<4 || isempty(Artefacts)
    
    Artefacts=cell(1,length(EEG));
end
if  nargin<3 || isempty(offset)
    offset=0;
end
if nargin<2 || isempty(window)
    window = EEG(1).Fs; % 1s window
end



for i = 1:length(EEG)
    % Epoch signal
    test = EEG(i);
    Fs = test.Fs;
    
    
    [Xte{i} isBad{i} ]= epoch_p300(test.s,test.Flash,window,offset,Artefacts{i});
    
    Yte{i}=test.Y;
    EEG(i)=test;
    Flash{i}=test.Flash;
    
end