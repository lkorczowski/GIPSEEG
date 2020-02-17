% pop_acstp - Compute the Adaptive Common Spatio-Temporal Pattern Filter (ACSTP) from
%             continuous or epoch EEG data.
%
% EEGOUT = pop_acstp(EEGIN);
%
% Input:
%    EEGIN - input EEGLAB dataset
%
% Output:
%    EEGOUT - updated EEGLAB structure
%
% Author: Louis Korczowski (EEGLAB interface by Arnaud Delorme)
%
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
% TO DO: 
% - ADD HISTORY
% - ADD OPTIONAL INPUTS
%
% see also : ACSTPshow, ACSTP

function [EEGOUT ACSTPstruct com] = pop_acstp(EEG)

if nargin < 1
    help pop_acstp;
    return;
end;
EEGOUT = EEG;
com    = [];

if nargin < 2

   cbevent = ['if ~isfield(EEG.event, ''type'')' ...
				   '   errordlg2(''No type field'');' ...
				   'else' ...
                   '   tmpevent = EEG.event;' ...
                   '   if isnumeric(EEG.event(1).type),' ...
				   '        [tmps,tmpstr] = pop_chansel(unique([ tmpevent.type ]));' ...
				   '   else,' ...
                   '        [tmps,tmpstr] = pop_chansel(unique({ tmpevent.type }));' ...
                   '   end;' ...
				   '   if ~isempty(tmps)' ...
				   '       set(findobj(''parent'', gcbf, ''tag'', tmpeventstr), ''string'', tmpstr);' ...
				   '   end;' ...
				   'end;' ...
				   'clear tmps tmpevent tmpeventstr tmpv tmpstr tmpfieldnames;' ];
               %%

   cbchan = ['if isempty(EEG.chanlocs)' ...
				   '   errordlg2(''No channel labels defined'');' ...
				   'else' ...
                   '   [tmps,tmpstr] = pop_chansel({ EEG.chanlocs.labels });' ...
				   '   if ~isempty(tmps)' ...
				   '       set(findobj(''parent'', gcbf, ''tag'', ''chans''), ''string'', tmpstr);' ...
				   '   end;' ...
				   'end;' ...
				   'clear tmps tmpevent tmpeventstr tmpv tmpstr tmpfieldnames;' ];
               
   geometry = { [2 1.5] 1   [2 1 0.5] [2 1 0.5] [2 1 0.5] [2 1 0.5] [2 1 0.5] 1 1 };
   geomvert = [ 1       0.5 1         1         1         1         1         1 1 ];
   uilist = { ...
       { 'style' 'text'       'string' 'Name for the new dataset' 'fontweight' 'bold' } ...
       { 'style' 'edit'       'string'  fastif(isempty(EEG.setname), '', [ EEG.setname ' epochs' ]) 'tag' 'name'} { } ...
       ...
       { 'style' 'text'       'string' 'Target time-locking event type(s) (default:all)' } ...
       { 'style' 'edit'       'string' '' 'tag' 'event1' } ...
       { 'style' 'pushbutton' 'string' '...' 'callback' [ 'tmpeventstr = ''event1'';' cbevent ] } ...
       ...
       { 'style' 'text'       'string' 'Non-target type(s) (default:all different from target)' } ...
       { 'style' 'edit'       'string' '' 'tag' 'event2' } ...
       { 'style' 'pushbutton' 'string' '...' 'callback' [ 'tmpeventstr = ''event2'';' cbevent ] } ...
       ...
       { 'style' 'text'       'string' 'Mask selection - electrodes ([]=all)' } ...
       { 'style' 'edit'       'string' '' 'tag' 'chans' } ...
       { 'style' 'pushbutton' 'string' '...' 'callback' cbchan } ...
       ...
       { 'style' 'text'       'string' 'Mask selection - time [start, end]' } ...
       { 'style' 'edit'       'string' ''  'tag' 'masktime'} ...
       { 'style' 'text'       'string' 'seconds' } ...
       ...
       { 'style' 'checkbox'   'string' 'Allow latency correction of at most' } ...  
       { 'style' 'edit'       'string' ''  'tag' 'latencycor'} ...
       { 'style' 'text'       'string' 'seconds' } ...
       ...
       { 'style' 'checkbox'   'string' 'Allow overlapping ERPs (continuous data only)'  'tag' 'erpoverlap'} ...
       ...
       { 'style' 'checkbox'   'string' 'Allow amplitude correction'  'tag' 'ampcor'} ...
       ...
       };
       %*LK: I can be wrong but the checkbox 'Allow latency correction of at most' doesn't have a tag

%        { 'style' 'text'       'string' 'Non-target time-locking event type(s)' } ...
%        { 'style' 'edit'       'string' '' 'tag' 'event2' } ...
%        { 'style' 'pushbutton' 'string' '...' 'callback' [ 'tmpeventstr = ''event2'';' cbevent ] } ...
%        ...
%        { 'style' 'text'       'string' 'Epoch limits [start, end]' } ...
%        { 'style' 'edit'       'string' '-1 2' 'tag' 'epoch'} ...
%        { 'style' 'text'       'string' 'seconds' } ...
%        ...
%   geometry = { geometry{:} [2 1 0.5] [2 1 0.5] };
%   geomvert = [ geomvert 1 1 ];
   
   [result userdat tmp res] = inputgui( 'geometry', geometry, 'uilist', uilist, 'helpcom', 'pophelp(''pop_epoch'')', 'title', 'Extract epochs with ACSTP filtering - pop_acstp()');
   if length(result) == 0 return; end;
end;

EEGCSTP.Fs         = EEG.srate;
%*LK: ligne suivante : attention ! Si EEG.data est un vector ligne [1 x EEG.nbchan*EEG.pnts*EEG.trials] , reshape
% ne va pas bien fonctionné les électrodes seront mélangées temporellement.
EEGCSTP.Channels   = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials)';


EEGCSTP.Channels(end+1,:) = 0;

% target trials
tmpEvent = EEG.event;
eventType = { tmpEvent.type };
if isempty(res.event1)
    indEpochs = 1:EEG.trials;
else
    res.event1 = parsetxt( res.event1 );
    indEpochs  = [];
    for iEvent = 1:length(res.event1)
        indEvent = strmatch(res.event1{iEvent}, eventType, 'exact');
        %*LK: si je comprends bien .epoch est l'index (en échantillon). Je
        %croyais que c'était .latency ?
        indEpochs = [ indEpochs EEG.event(indEvent).epoch ];
    end;
    if isempty(indEpochs), errordlg2('No epoch found with this target event type'); end;
end;
if isempty(res.event2)
    indNonEpochs = setdiff(1:EEG.trials, indEpochs);
else
    res.event2 = parsetxt( res.event2 );
    indNonEpochs  = [];
    for iEvent = 1:length(res.event2)
        indEvent = strmatch(res.event2{iEvent}, eventType, 'exact');
        indNonEpochs = [ indNonEpochs EEG.event(indEvent).epoch ];        
    end;
    if isempty(indNonEpochs), errordlg2('No epoch found with this non-target event type'); end;
end;

allEpochs = sort( [ indEpochs(:);indNonEpochs(:) ]);
if length(unique(allEpochs)) ~= length(allEpochs)
    errordlg2('Target and non-target selection overlap');
end;
EEGCSTP.Trigger = zeros(1, size(EEGCSTP.Channels,1));
% ACSTP options
%*LK: if allEpochs are the indices (in sample), I do not understand why we
%need to multiply by the length of the epoch EEG.pnts +1 ???
% wouldn't be EEGCSTP.Trigger(allEpochs)=1; instead ?
% old line : EEGCSTP.Trigger((allEpochs-1)*EEG.pnts+1) = 1;
EEGCSTP.Trigger(allEpochs) = 1;

EEGCSTP.EpochClass = zeros(1, length(allEpochs));
for iEpoch = 1:length(allEpochs)
    if any(allEpochs(iEpoch) == indEpochs)
        EEGCSTP.EpochClass(iEpoch) = 1;
    end;
end;

% ACSTP options
%*LK: latencycor doesn't exists (inside res). Line corrected
% if isempty(res.latencycor), res.latencycor = 0; else res.latencycor = str2num(latencycor); end;
%*LK: if latencycor==0 then, the latency correction is desactivated.
if isempty(res.latencycor), res.latencycor = 0; else res.latencycor = str2num(res.latencycor); end;
if isempty(res.chans)       chanInds = 1:EEG.nbchan; else chanInds = eeg_chaninds(EEG, parsetxt(res.chans)); end;
if isempty(res.masktime) 
    %*LK: if continuous EEG, EEG.pnts will be the size of the entire
    %recording.
    maskTime = 1:EEG.pnts;
else
    timeLim = str2num(res.masktime);
    if length(timeLim) ~= 2, error('Wrong mask time limits'); end;
    maskTime = eeg_lat2point(timeLim, 1, EEG.srate, [EEG.xmin EEG.xmax])
     %*LK: maskTime here has only 2 sample [beginning end] however, in ACSTP
    %we ask to the user to put exhautively the samples indexes such as
    % maskTime= 5:55 (samples between 5 and 55) because the user want maybe
    % discontinuities ( e.g. [5:30 40:55])
    % here is a proposition accordingly :
    maskTime = ceil(maskTime(1)):floor(maskTime(2));
end;
    
%*LK: if continuous EEG, EEG.pnts will be the size of the entire
    %recording. Same problem that line 171. I hardcoded 128 to test the other problems
    %maybe we would like to consider a box in the GUI of pop_acstp ?
ACSTPoptions.Epoch_size      = 128;warning('Epoch_size option hardcoded')%EEG.pnts; 
ACSTPoptions.LatencyCorr_max = ceil(res.latencycor*EEG.srate);
ACSTPoptions.Mask_Electrodes = chanInds;

ACSTPoptions.Mask_Time       = maskTime;
%*LK: added the option for the weights estimation should be 0 or 1
ACSTPoptions.Weights=res.ampcor;
%*LK: if the user what to set the Subspace dimension to evaluate.
% need to add the box 'subspace dimension
% (by default it will try to find the best subspace)'.
 ACSTPoptions.SubspaceDim= [16:-1:10];warning('SubspaceDim option hardcoded');
%
%*LK: the user, to go faster, could want to limit the latencies estimation
%to the target class (TA). (by default the latencies will be compute for
%every classes).
% we can use a 'callback' [ 'tmpeventstr = ''event1'';' cbevent ] for that
  ACSTPoptions.computeClassLat=[1];warning('computeClassLat option hardcoded');
ACSTPoptions.overlap=res.erpoverlap;

[ res ACSTPstruct ] = ACSTP(EEGCSTP, ACSTPoptions);
EEGOUT = EEG;
EEGOUT.data = res;

 EEGCSTP.ElectrodesName = { EEG.chanlocs.labels };
 figure
 ACSTPshow(EEGCSTP,ACSTPoptions,ACSTPstruct);


% write history TO DO
%com = sprintf('EEG = pop_acstp(EEG, %s)', vararg2str(


% TO DO
%*LK: corriger bug ACSTPshow 

    