function [ ] = writeInfo( filename, EEG )
% this function generate a txt file 
% it writes in filename the informations about the EEG
%
% INPUTS
% ************************************
% filename:  the full path and file for the EEG txt file. Each EEG
%               sample will be delimited by ' '. Each electrod signal is 
%
% EEG structure (gdf) such as :
% 'EEG.Channels':      the EEG signals dim [nb_samples x nb_electrods]
% 'EEG.Trigger':  binary vector for stimulation dim [1 x nb_samples] 0 for no Flash, 1
%               for Flash. There are N 1s in 'EEG.Flash'
% 'EEG.EpochClass':      class' tag vector of the stimulation dim [N x 1] 0 for "non-Target", 1 for "Target".
% 'EEG.Fs':     scalar for the sample rate (i.e 128,256, or 512Hz)
% 'EEG.offset': a scalar of vector [N x 1] containing an offset if
%               necessary. By default: 0 for each stim.


    FID = fopen( filename, 'w' );
    fprintf(FID,'EEG file (EEG.txt) and Stimulation file written (STIM.txt) \r\n');
    fprintf(FID,[datestr(now) '\r\n\nContent:\n']);
    fprintf(FID, ['Nb Samples:\t' num2str(size(EEG.Channels,1)) ' (approx. %d seconds)\n'],size(EEG.Channels,1)/EEG.Fs);
    fprintf(FID, ['Nb Electrodes:\t' num2str(size(EEG.Channels,2)) '\n']);
    fprintf(FID, ['Sample Rate:\t' num2str(EEG.Fs) ' Hz\n']);
    fprintf(FID, ['Nb Epochs:\t' num2str(length(EEG.EpochClass)) '\n\n']);
    fprintf(FID, ['EEG.txt contains the EEG samples separated by '' ''. \nEach electrode is consecutive to each other\n']);
    fprintf(FID, ['STIM.txt contains the stimulations position, stimulation offset, and stimulation tag (NT/TA)\n\n']);
    fprintf(FID, ['related matlab functions: EEG_mat2txt.m, writeStim.m, writeInfo.m']);
    fclose(FID);
end

