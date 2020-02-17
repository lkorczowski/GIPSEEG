function [Trigger StimCode StimPos]= triggers_continuous2discrete(StimContinuous)
%% [Trigger StimCode StimPos]= triggers_continuous2discrete(StimContinuous)
%
% INPUTS
%  StimContinuous: [nb samples x1 ] the Stimulations codes in continuous
%                  channel, possibly as step signals when the stimulations
%                  codes are maintained during several samples.

%         Trigger: [nb samples x1 ] The discrete Trigger channel of '0' with '1' at the start
%                   of each sweep/flash/stimulation. There are [nb epochs] '1'.
%        StimCode: [nb epochs x1] class code of the sweeps
%         StimPos: [nb epochs x1] position of the sweeps (i.e.
%                  find(StimDiscrete))

Trigger=[0; diff(StimContinuous)>0];% put a zero at the first sample
StimPos=find(Trigger);
StimCode=StimContinuous(StimPos);

