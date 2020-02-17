function SCHUFFLE_INDEX = uniqueShuffle2(N_SHUFFLE, NB_WIN)
% function  SCHUFFLE_INDEX = uniqueShuffle2(N_SHUFFLE, NB_WIN)
%
% ************************************************************************
% This function return a random permutation matrix of size 
% (N_SHUFFLE,NB_WIN) with only unique rows filled with random indexes 
% in range [1->NB_WIN].
% Main purpose is to compute statistical permutation tests between two 
% vectors of identical length NB_WIN.
%
% 
%
% Input:
% ------
% N_SHUFFLE     --> number of shuffling requested
% NB_WIN        --> scalar with length of vectors to shuffle 
%
% Output:
% -------
% SCHUFFLE_INDEX    --> matrix of shuffling indexes with decimal values (N_SHUFFLE,NB_WIN)
%
% History:
% --------
% *** 2012-01-14
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


if(N_SHUFFLE > (2^NB_WIN))
    error(['Requested number of permutations exceed max possible perm with ' int2str(NB_WIN) ' time windows']);
% handle the particular case where all possible different permutations must be found
elseif (N_SHUFFLE == (2^NB_WIN))
    SCHUFFLE_INDEX = dec2bin(0:(2^NB_WIN)-1);   
else 
    SCHUFFLE_INDEX      = zeros(N_SHUFFLE,NB_WIN);
    SCHUFFLE_INDEX(1,:) = randperm(NB_WIN); 
    for p_ix = 2:N_SHUFFLE
    %         disp(['progress: ' int2str(100*p_ix/NB_WIN) '%'])
        allUnique = false;
        while(~allUnique) 
            % create random vector with unique values
            SCHUFFLE_INDEX(p_ix,:) = randperm(NB_WIN); 

            % reiterate until this vector is unique within SCHUFFLE_INDEX
            sameShuff_prev = false(1,NB_WIN);
            sameShuff_prev(1:p_ix-1) = true;
            for win_ix = 1:NB_WIN
                sameShuff_curr = (SCHUFFLE_INDEX(sameShuff_prev,win_ix) == SCHUFFLE_INDEX(p_ix,win_ix));
                if ~any(sameShuff_curr)
                    allUnique = true;
                    break;
                end
                sameShuff_prev = sameShuff_curr;
            end
        end
    end
end


