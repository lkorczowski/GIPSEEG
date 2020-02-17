    function OutStim=DecimateTrigger(InStim,decimationfactor,SizeFinal)
    
% Decimate an single dimension trigger channel InStim
%
% OutStim=DecimateTrigger(InStim,l_decimationfactor,SizeFinal)
% 
%
%
% INPUTS :
% ------
%          
%          InStim: [nb samples x1 ] Trigger channel with a scaler for each
%                                   stimulation (0 otherwise)
%  decimationfactor: [scalar] the decimation factor
%      SizeFinal*: [scalar] the wanted outlength (default: length of InStim
%      rounded)
%
% *optional
%
% OUTPUT :
% ------
% OutStim with decimated trigger signal
%
%
% *** History: 16-Apr-2015
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
    if nargin<3 | isempty(SizeFinal)
        SizeFinal=length(InStim);
    end
    
        l_indF=find(InStim);
        l_indFdec=round(l_indF/decimationfactor);
        if (any(l_indFdec==0)|any(l_indFdec>SizeFinal))
            warning('DecimateTrigger: Trigger detected at the edge, the stimuli indexes can be biased')
            l_indFdec(l_indFdec==0)=1;
            l_indFdec(l_indFdec>SizeFinal)=SizeFinal;
        end
        OutStim=zeros(SizeFinal,1); %the stimulation channel is the size of the decimated data
        OutStim(l_indFdec)=InStim(l_indF); %put the stimcodes at the good place
    end