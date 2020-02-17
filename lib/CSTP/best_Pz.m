function [Pzind SNR]=best_Pz(Zbarz,Electrodes,Window,indClass)
% [Pzind SNR]=best_Pz(Zbarz,Electrodes,Window,indClass)
% find the best subspace dimension
% it is important that the Pz are sorted in descend order in Zbarz
% *** History: 2015-03-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
% see also: ACSTP
if nargin<4
    indClass=size(Zbarz{1},3); %take the "target" class that should be the last
end

for i=1:length(Zbarz)
    noise=Zbarz{i}(:,:,indClass);
    signal=noise(Electrodes,Window);
    SNR(i)=norm(signal,'fro')/norm(noise,'fro'); %%%%% IN PAPER: squared %%%%
end
        SNR=SNR-min(SNR);%normalize SNR

%% Find the maximum
% argmax_Pzind(SNR(Pzind))
if length(Zbarz)>3
    [maxs INDtmp]=findpeaks(SNR);
else
    INDtmp=[];
end
    [tmp Pzind]=max(SNR);
if ~isempty(INDtmp)
    IND=INDtmp(maxs>(tmp*0.66));
    if ~isempty(IND)
    Pzind=IND(1);
    end
end