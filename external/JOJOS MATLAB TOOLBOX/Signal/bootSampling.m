function [BOOTMAT] = bootSampling(TOTAL_SAMPLE, BOOT_SAMPLE, N_BOOT)
%function [BOOTMAT] = bootSampling(TOTAL_SAMPLE, BOOT_SAMPLE, N_BOOT)

% ************************************************************************
% This function computes a matrix of indexes for bootstrap sampling.
%
% Input:
% ------
% TOTAL_SAMPLE 	--> scalar, length of data on which to perform bootstrap sampling
% BOOT_SAMPLE	--> scalar, nb of samples to draw at each bootstrap sampling
% N_BOOT        --> scalar, nb of random sampling to perform      
%
% Output:
% -------
% BOOTMAT   	--> matrix with N_BOOT random sampling in interval [1:TOTAL_SAMPLE]
%                   format: (N_BOOT,BOOT_SAMPLE). Each row has unique values 
%                   only, but same row may be sampled multiple times
%
% History:
% --------
% First version: 2013-06-20
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

BOOTMAT = zeros(N_BOOT,BOOT_SAMPLE);
for boot_ix = 1:N_BOOT
    uniqueVals = 0;
    while(uniqueVals==0)
        BOOTMAT(boot_ix,:) = randi(TOTAL_SAMPLE,1,BOOT_SAMPLE);
        uniqueVals = (length(unique(BOOTMAT(boot_ix,:))) == BOOT_SAMPLE);
    end 
end

end

