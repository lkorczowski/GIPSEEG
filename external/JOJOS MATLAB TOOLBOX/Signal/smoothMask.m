function [SMOOTH_MASK] = smoothMask(MASK, SMOOTH, MIN)
% function [SMOOTH_MASK] = smoothMask(MASK, SMOOTH, MIN)
% ************************************************************************
%  Weight binary mask with smoothing windows to remove sharp
%  discontinuities. 
%
% INPUTS
% - MASK        --> vector of binary values or logical (1 x Ns)
% - SMOOTH      --> (optional)  specifies nb of samples to smooth at begin 
%                               and end of masked periods (with value = 1)
%                               default: 8 samples
% - MIN         --> (optional)  specifies nb of samples minimim for
%                               masked periods to be considered 
%                               (shorter periods are set to zero)
%                               default: 2 samples
%
% OUTPUTS
% - SMOOTH_MASK	 --> vector of double values (1 x Ns)
%
% HISTORY 
% Last version: 10/04/2013
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

if(nargin<2)
    SMOOTH = 8;
elseif(nargin<3)
    MIN = 2;
end

winSmooth = hann(2*SMOOTH+2)'; 
winSmooth_beg = winSmooth(2:SMOOTH+1);
winSmooth_end = winSmooth(SMOOTH+2:end-1);
SMOOTH_MASK	= zeros(1,length(MASK));    % weights for smoothing data at edged of masked periods
SMOOTH_MASK(MASK) = 1;
sig_beg = find(diff(MASK) == 1) +1;     % index of begin of masked signal periods
sig_end = find(diff(MASK) == -1);       % index of end of masked signal periods
if length(sig_end) < length(sig_beg)    % here MASK == 1 until end sample, so we have to add it manually
    sig_end = [sig_end length(MASK)];
end
for sig_ix = 1:length(sig_beg)
    sig_length = sig_end(sig_ix)-sig_beg(sig_ix)+1;
    if(sig_length >= 2*SMOOTH)
       SMOOTH_MASK(sig_beg(sig_ix):sig_beg(sig_ix)+SMOOTH-1) = winSmooth_beg;
       SMOOTH_MASK(sig_end(sig_ix)-SMOOTH+1:sig_end(sig_ix)) = winSmooth_end;
    elseif(sig_length > MIN)  % here masked period is too short so we have to crop smoothing 
       SMOOTH_MASK(sig_beg(sig_ix):sig_beg(sig_ix)+ceil(sig_length/2)-1) = winSmooth_beg(1:ceil(sig_length/2));
       SMOOTH_MASK(sig_end(sig_ix)-floor(sig_length/2)+1:sig_end(sig_ix)) = winSmooth_end(end-floor(sig_length/2)+1:end);
    else 
       SMOOTH_MASK(sig_beg(sig_ix):sig_end(sig_ix)) = 0;
    end
end
    
end