function [Phi_in, Phi_inout] = procFreqFit(Data, Faxis, FoI, Delta_f)
% This function computes the frequency fitting (to FoI) of signals with
% spectral content (Data)
%
% INPUTS
% - Data            --> dependence measure, format can be:
%                       - 2D: IC * freq
%                       - 3D: IC * IC * freq
% - Faxis           --> index of frequencies corresponding to Data
% - FoI             --> Frequencies of Interest (bands and single freq)
%                       eg.: for freq 5, 7 and 10Hz alone -> FoI = [5 7 10]
%                       ATTENTION:  these freqs should be present in Faxis
%                                   otherwise they are not considered!
% - Delta_f         --> (optionnal) specify size of frequency bands around
%                       FoIs. Eg.: Delta_f=.5 --> 1Hz constant bands centered on FoIs.
%                       Default value: 0 (consider only single freqs)
%                       
%
% OUTPUTS
% - Phi_in          --> dependence: (sum in-FoI)
% - Phi_inout       --> dependence: (sum in-FoI) / (total)
%               
% HISTORY 
% Last version:  11/11/2011
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

% Default freq band
if nargin < 4   Delta_f = 0; end

% Calculate indexes of FoI and non-FoI
InFreqIndex = zeros(size(Faxis));
for i = 1:length(FoI) % for each FoI
    InFreqIndex((Faxis >= FoI(i)-Delta_f) & (Faxis <= FoI(i)+Delta_f)) = 1;
end
InFreqIndex = logical(InFreqIndex);
if(~any(InFreqIndex))
    error('[procFreqAdjust] input arguments: no FoI in Faxis!')
end
OutFreqIndex = (ones(1,length(Faxis)) ~=InFreqIndex);


% 2D mode (IC*freq)
if(length(size(Data)) == 2)
    if(length(Faxis) ~= size(Data,2))
        error('[procFreqAdjust] input arguments: Faxis does not match Data!')
    end
    % Process frequency adjustement
    Phi_in      = sum(Data(:,InFreqIndex),2)';
	Phi_inout   = sum(Data(:,InFreqIndex),2)' ./ sum(Data,2)' ;
    
elseif(length(size(Data)) == 3)
    if(length(Faxis) ~= size(Data,3))
        error('[procFreqAdjust] input arguments: Faxis does not match Data!')
    end
    % Process frequency adjustement
    Phi_in      = sum(Data(:,:,InFreqIndex),3) ./ length(find(InFreqIndex));
    Phi_inout   = sum(Data(:,:,InFreqIndex),3) ./ sum(Data,3); %mean(Data(:,:,OutFreqIndex),3); %
    % same with geometrical mean
%     Phi_inout   = sum(Data(:,:,InFreqIndex),3) ./ (prod(Data(:,:,OutFreqIndex),3)^(1/length(OutFreqIndex)));     
else
    error('[procFreqFit] Data has a wrong format')
end






