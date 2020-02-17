% loadONLY(FILEpath) - SCRIPT TO COPY TO HAVE EXPECTED RESULT
% load all variables from FILEpath.mat only if non-existing in workspace
FILEpath='file_path';

%compare workspace and file
VARs=who;VARs2=who('-file',FILEpath);
%load non-existing variable from file
for indv=1:length(VARs2);if ~exist(VARs2{indv},'var');load(FILEpath,VARs2{indv});end;end