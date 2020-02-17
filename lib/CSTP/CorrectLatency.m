function FlashCorrected=CorrectLatency(Flash,Latencies,RND)
% FlashCorrected=CorrectLatency(Flash,Latencies,RND)
%
% take the indice of the sweeps (Flash) and the Latencies (expressed in
% samples) and give the corrected indice of the sweeps (FlashCorrected)
% RND is the random seed if you had randomized the Flashs (useful for
% bootstrap analysis).
%
% *** History: 2015-03-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
% see also : LatencyEstimation


    if nargin<3
       RND=1:length(find(Flash)); %no random seed
    end
    if size(Flash,1) == 1, Flash = Flash'; end;
    FlashCorrected=zeros(size(Flash));
    indF=find(Flash);
    if length(indF)~=Latencies
        indF(RND)=indF(RND)+Latencies;
        indFcorrected=indF;
    else
    indFcorrected=indF+Latencies;
    end
    FlashCorrected(indFcorrected)=1;
    

end