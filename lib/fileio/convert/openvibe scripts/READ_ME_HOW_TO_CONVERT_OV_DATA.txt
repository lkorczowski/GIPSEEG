-you need to have openvibe installed
-you need to have ruby installed (optional, however it's important to do the process automatically) and openvibe-gem ruby

############### MANUAL ###############################
-open gdf_converter.xlm in openvibe designer
-modify Generic stream reader in order to open the .ov file
-modify GDF written to save it where you want

############## IN MATLAB #############################
-for the gdf loading function you need the biosig toolbox the folder t200 access file is needed
-use load_EEG_data.m for automatic convertion from gdf to mat. For manual use [s h] = sload(file);
-convert h to the session structure with
% EEG is a structure with
%             Fs: scalar (sample rate in Hz)
%       StimCode*: [nb epochs x1 double] (useless for CSTP)
%     StimCodePos*: [nb epochs x1 double] (useless for CSTP) 
%           Flash: [nb samples x1 double] position of the sweeps
%          Target*: [nb samples x1 double] position of the sweeps of class
%          TARGET
%               Y: [nb epochs x1 double] class of the sweeps (0 for
%               Non-TARGET, 1 for TARGET)
%               s: [nb samples x nb channels double] raw EEG recording
%               h*: [1x1 struct] software and hardware information (useless for CSTP)
%   *optionnal



-before analysing the data, use preprocessingEEG
-check the data with script_ACSTP (modify carefully the EEG structure and do not preprocess if you did it already)

