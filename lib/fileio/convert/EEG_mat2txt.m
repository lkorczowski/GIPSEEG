function EEG_mat2txt(filepath, EEG)
% EEG_mat2txt(filepath, EEG)
%
% INPUTS
% ****************************
% filepath:     the full path and file for the EEG txt file. Each EEG
%               sample will be delimited by ' '. Each electrod signal will
%               be concatenated.
%
% EEG structure (gdf) such as :
% 'EEG.Channels':      the EEG signals dim [nb_samples x nb_electrods]
% 'EEG.Trigger':         binary vector for stimulation dim [1 x nb_samples] 0 for no Flash, 1
%                      for Flash. There are N 1s in 'EEG.Trigger'
% 'EEG.EpochClass':             class' tag vector of the stimulation dim [N x 1] 0 for "non-Target", 1 for "Target".
% 'EEG.Fs':            scalar for the sample rate (i.e 128,256, or 512Hz)
% 'EEG.offset':        a scalar or vector [N x 1] containing the offset(s) if
%                      necessary. By default: 0 for each stim.
%
% OUTPUTS
% ****************************
% in the filepath folder, 3 files will be generated:
% EEG.txt
% STIM.txt
% INFO.txt (read this file for details)
%
% exemple: script_MAT_to_TXT
% see also: writeStim, writeInfo
%
% *** 2015-03-09
% Created by L. KORCZOWSKI @GIPSA-Lab

if ~isfield(EEG,'offset')
    EEG.offset=0;
end
if length(EEG.offset)==1
    disp(['offset has been set to ' num2str(EEG.offset)])
    EEG.offset=ones(length(EEG.EpochClass),1)*EEG.offset;
end

tic
indFlash=find(EEG.Trigger);
if size(EEG.EpochClass,1)<size(EEG.EpochClass,2)
    EEG.EpochClass=EEG.EpochClass';
end
stim=[indFlash,EEG.offset,EEG.EpochClass];
FOLDER=fileparts(filepath);
if ~exist(FOLDER,'dir')
mkdir(FOLDER);
disp([FOLDER ' folder has been created'])
end
disp('Writting...')
dlmwrite(fullfile(filepath,'EEG.txt'),EEG.Channels(:)','delimiter',' ', 'precision','%.6f'); %write EEG signals
disp('EEG.txt written')
writeStim( fullfile(filepath,'STIM.txt'), stim ); %write stim signals
disp('STIM.txt written')
writeInfo( fullfile(filepath,'INFO.txt'), EEG );
disp('INFO.txt written')
disp(['all files have been written successfully in ' FOLDER])


