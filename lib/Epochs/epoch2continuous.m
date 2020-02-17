function [s TriggerC]=epoch2continuous(X,Trigger,Format)
% Project back the data available into a epoched EEG 3D matrix X
% to a continuous EEG 2D matrix s with the help of a Trigger channel.
%
%
% [s]=epoch2continuous(X)
% Assume that the epochs are non-overlapping
%               X: [nb channels x nb samples x nb trials] epoched EEG
%               s: [nb channels x nb samples*nb trials] continuous EEG
%
% [s]=epoch2continuous(X,Trigger)
% Assume that the epochs are overlapping
%               X: [nb channels x nb samples epoch x nb trials] epoched EEG
%         Trigger: [nb samples continuous x 1] the trigger channel which is
%                  0 everywhere but at the begginning of each epoch contained in the
%                  third dimension of X. There are [nb trials] '1'.
%               s: [nb channels x length(Trigger)] continuous EEG
%
% [s TriggerC]=epoch2continuous(X,Trigger,Format)
% Assume that the epochs are overlapping
%          Format: string, 'compact' or 'full' (default). If 'full', the
%                  output will be the same length that Trigger. If  'compact', the
%                  continuous signal will be removed from all the gap between
%                  epochs (normaly zeros with the 'full')
%        TriggerC: the compacted EEG trigger

N=size(X,1);
T=size(X,2);
if(N>T),warning('There are more channels than samples in the epochs'); end
K=size(X,3);
s=nan(N,length(Trigger));

if nargin<3
    Format='full';
end

if nargin<2 %assuming non-overlapping epochs
    Trigger=zeros(T*K,1);
    Trigger(1:T:end)=1;
end
TriggerC=Trigger;
StimPos=find(Trigger)
if length(StimPos)~=K, error('The Number of epoch in X in incompatible with the number of element in Trigger');end

% the following loop could be optimize to avoid writting several time the
% same values. If the Trigger channel is wrong, s will have discontinuities
% in the EEG (a solution could be to smooth the signal thereafter)
for indK=1:K, s(:,StimPos(indK):StimPos(indK)+T-1)=X(:,:,indK);end

% for a compact version, we can remove all the gap (nan) between trials
if strcmp(Format,'compact')
    toRemove=isnan(s(1,:));
    s(:,toRemove)=[];
    TriggerC(toRemove)=[];
else
    %for non-compact version
    s(:,isnan(s(1,:)))=0;
end
