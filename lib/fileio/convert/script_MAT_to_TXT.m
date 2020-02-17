 clear all
load( ['ss1_ERWAN_RAW.mat']) %your .mat file containing EEG and signals data
EEG=data.session{1}.phase{1}.online
% EEG structure such as :
%
% EEG structure (gdf) such as :
% 'EEG.Channels':      the EEG signals dim [nb_samples x nb_electrods]
% 'EEG.Trigger':         binary vector for stimulation dim [1 x nb_samples] 0 for no Flash, 1
%                      for Flash. There are N 1s in 'EEG.Trigger'
% 'EEG.EpochClass':             class' tag vector of the stimulation dim [N x 1] 0 for "non-Target", 1 for "Target".
% 'EEG.Fs':            scalar for the sample rate (i.e 128,256, or 512Hz)
% 'EEG.offset':        a scalar or vector [N x 1] containing the offset(s) if
%                      necessary. By default: 0 for each stim.

% functions needed : writeStim, preprocessingEEG (optional)
% Louis Korczowski, GIPSA-lab, 2015

f1=0.25 ;%lowcut freq
f2=30 ;%highcut freq
N=4;%filter order
decimationN=4; %decimation factor

tic
FULLDIR = '.\txt\'; % PUT output directory here

fs=EEG.Fs;
E=EEG.s;
Flash=EEG.Flash;

% if needed, you can use preprocessing with bandpass filter with lowcut
% frequency 'f1', highcut frequency 'f2', Order 'N' and decimation factor
% 'decimationN'
% Too avoid the preprocessing, just save the output to the corresponding
% input (E->Ef, fs->fsf, flashf->Flash)
[EEG.s EEG.Fs EEG.Flash]=preprocessingEEG(EEG.s,EEG.Fs,[f1 f2 N decimationN],EEG.Flash);
EEG.offset=0;

EEG_mat2txt(FULLDIR,EEG)

