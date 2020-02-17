function sLoreta_loc2(A, LeadField_file, powSource_file)
% function sLoreta_loc2(A, LeadField_file, powSource_file)

% Source localization from BSS mixing matrix by sLoreta inverse problem 
% resolution.
% ************************************************************************
%
% Inputs:
% ------
%   A     	 -->  Mixing matrix (columns are scalp maps)
%   transMat_file   -->  Transformation matrix matrix filename 
%                        (must include path). File should be created using:
%                           1)  LORETA software -> lead field matrix 'a'
%                           2)  Marco's tool -> binary version with no 
%                               header of transformation matrix
%   powSource_file  -->  (optionnal) Filename used for writing and subsequent 
%                        reading of source power (must include path). 
%                       
%
% References
% -------
% Congedo (2006) 
% Pascual-Marqui (2002)
% Pascual-Marqui (2007)
%
% History:
% --------
% created by Guillaume Lio @ GIPSA-Lab, November 2009
% modified by Jonas Chatel-Goldman @ GIPSA-Lab,February 2013


    if(nargin < 3)
        path           = fileparts(LeadField_file);
        powSource_file = [path '\sLORETA_POW.lor'];
    end

    %  Transformation matrix file opening...
    s           = size(A);    	% Size of the BSS mixing matrix
    Nb_voxels   = 2394;         % Number of voxels in the head model for inverse problem resolution         
    Nb_sources  = s(2);         % Number of sources
    Nb_electrodes = s(1);       % Number of electrodes
        % this file has a header with format:
        % int32 (line number) / int32 (column number) / float32 / float32
        % i.e. 4 + 4 + 4 + 4 = 16 bytes
    [fid,msg] = fopen(LeadField_file,'r'); 
    fseek(fid,16,'bof'); % skip header
    T = fread(fid, [Nb_electrodes, Nb_voxels * 3], 'float32');   
    fclose(fid);
	T = T';     % T is N_voxels*3 by N_chan matrix
 
    % inverse solution sLoreta Calculation                     
    disp('Inverse solution calculation...');
    SourcesPOW = zeros(Nb_sources, Nb_voxels);
    TA = T * A;     % format: (Nb_voxels * 3, Nb_sources)
    for m = 1 : Nb_sources
        for i = 1 : (Nb_voxels * 3)
            SourcesPOW(m, (fix((i-1)/3) +1)) = SourcesPOW(m, (fix((i-1)/3) +1)) + TA(i, m).^2;
        end
    end
    clear 'TA' 'T' ;   

    % Loreta key viewer
    disp('Source localization visualization...');
    dlmwrite(powSource_file, SourcesPOW, 'delimiter', '\t', 'precision', 6);
    winopen(powSource_file);
end
