clear all
load( [Directory 'FullErwanData_RAW.mat']) %your .mat file containing EEG and signals data
% EEG structure such as :
% 'EEG.s' the EEG signals dim [nb_electrodes x nb_samples]
% 'EEG.Flash' binary vector for stimulation dim [1 x nb_samples] 0 for no Flash, 1
% for Flash. There are N 1s in 'EEG.Flash'
% 'EEG.Y' : class' tag vector of the stimulation dim [N x 1] 0 for "non-Target", 1 for "Target".
% 'EEG.fs' : scalar for the sample rate (i.e 512Hz)

% functions needed : writeStim, preprocessingEEG (optional)
% Louis Korczowski, GIPSA-lab, 2015

f1=1 ;%lowcut freq
F2=20 ;%highcut freq
N=4;%filter order
decimationN=4; %decimation factor

tic
FULLDIR = ''; % PUT output directory here
filenameEEG=['\EEG.txt'];
filenameSTIM=['\STIMnew.txt'];

fs=EEG.Fs;
E=EEG.s;
Flash=EEG.Flash;

% if needed, you can use preprocessing with bandpass filter with lowcut
% frequency 'f1', highcut frequency 'f2', Order 'N' and decimation factor
% 'decimationN'
% Too avoid the preprocessing, just save the output to the corresponding
% input (E->Ef, fs->fsf, flashf->Flash)
[Ef fsf flashf]=preprocessingEEG(E,fs,[f1 f2 N decimationN],Flash);


indFlash=find(flashf);
stim=[indFlash,EEG.Y];

toc
dlmwrite([FULLDIR filenameEEG],EEG(:)','delimiter',' ', 'precision','%.6f') %write EEG signals
writeStim( [FULLDIR filenameSTIM], stim ) %write stim signals
toc
disp('File written')

