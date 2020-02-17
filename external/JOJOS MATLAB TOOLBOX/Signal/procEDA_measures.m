function [EDA_measures] = procEDA_measures(EDA, EVENTS, Fs, display_mode)
% function [EDA_measures] = procEDA_measures(EDA, EVENTS, Fs, display_mode)
%
% ************************************************************************
% Process common measures from Electrodermal Activity (EDA) signals.
% In particular, look at Event-Related Skin Conductance Responses (ER-SCR)
% using EVENTS vector.
% /!\ THESE ARE REALLY "HOME-MADE" MEASUREMENTS. USE IT WITH CAUTION.
%
% Input:
% ------
% EDA           --> input timeserie (nbChannels x nbSamples) (s)
% EVENTS        --> matrix with event time samples (nbChannels x nbEvents), used to measure SCRs
% Fs            --> (optional) EDA frequency rate (Hz)
% display_mode  --> (optional) if set to '1', various display 
%
% Output:
% -------
% EDA_measures  --> amplitude of ER-SCR, (nbChan x nbEvents)
%                   defined simply as max(EDA(0-3s0sec))-EDA(0), so the estimated 
%                   amplitudes are strictly positive and easy to compute
% ?..             
%
% History:
% --------
% Last version:  2013-08-01
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com
% based on (Boucsein et al., 2012) 

% assess inputs
if (nargin < 4) || isempty(display_mode)
    display_mode = 0;
end
if (nargin < 3) || isempty(Fs)
    Fs = 128;
end
if iscolumn(EVENTS)
   EVENTS = EVENTS'; 
end
if size(EDA,1) ~= size(EVENTS,1)
    error('[procEDA_measures] Not the same number of channels in EDA and EVENTS.')
end
nbChan = size(EDA,1);
nbEvents = size(EVENTS,2);

% estimation of ER-SCR amplitudes 
AMP = zeros(nbChan,nbEvents);
for event_ix = 1:nbEvents
    inital_val      = indexMat(EDA,1:nbChan,EVENTS(:,event_ix));
    run_begin_ix    = matFromVecs(EVENTS(:,event_ix),EVENTS(:,event_ix)+30*Fs);
    data_run_begin  = indexMat(EDA,1:nbChan,run_begin_ix);
    AMP(:,event_ix) = max(data_run_begin,[],2) - inital_val;
end

EDA_measures = AMP;




function [mat] = matFromVecs(v1,v2)
    mat = zeros(length(v1),v2(1)-v1(1));
    for i = 1:length(v1)
        mat(i,:) = v1(i):v2(i)-1;
    end
end


end



