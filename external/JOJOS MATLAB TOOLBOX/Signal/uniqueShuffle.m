function SCHUFFLE_INDEX = uniqueShuffle(N_SHUFFLE, VEC_LENGTH, SYMETRICAL)
% function  SCHUFFLE_INDEX = uniqueShuffle(N_SHUFFLE, VEC_LENGTH, SYMETRICAL)
%
% ************************************************************************
% This function return a random permutation matrix of size 
% (N_SHUFFLE,VEC_LENGTH) with only unique permutations. 
% Main purpose is to compute statistical permutation tests between two 
% vectors of identical length VEC_LENGTH.
%
% This problem can be resumed as finding a combination of N_SHUFFLE unique
% integer values drawn from uniform distribution on the open interval: 
% [0 (2^VEC_LENGTH)-1]
% 
%
% Input:
% ------
% N_SHUFFLE     --> number of shuffling requested
% VEC_LENGTH    --> scalar with length of vectors to shuffle (max 53 due to 
%                   randi function)
% SYMETRICAL    --> (optional) if set to 0, do not consider mirrored 
%                   permutations (e.g. 101 and 010), i.e. draw unique 
%                   values on half-interval [0 (2^VEC_LENGTH-1)-1]
%
% Output:
% -------
% SCHUFFLE_INDEX --> binary matrix of shuffling indexes with logical values (N_SHUFFLE,VEC_LENGTH)
%
% History:
% --------
% *** 2012-12-03
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


% default mode
if(nargin < 3)||(isempty(SYMETRICAL)) SYMETRICAL = 1;  end

% 1) calculate max number of permutations possible
if SYMETRICAL == 1 
    Max_shuffling = 2^(VEC_LENGTH); 
else
    Max_shuffling = 2^(VEC_LENGTH-1); 
    VEC_LENGTH = VEC_LENGTH-1;
end
allUnique = false;
if(N_SHUFFLE>Max_shuffling)
    disp(['[uniqueShuffle] --- WARNING --- ']);
    disp(['Requested number of permutations (' int2str(N_SHUFFLE) ...
            ') is >= to max possible permutations (' int2str(Max_shuffling) ').']);
    disp(['--> Forcing permutation number to ' int2str(Max_shuffling)]);
    N_SHUFFLE = Max_shuffling;
    permutDec = [0:N_SHUFFLE-1];  % handle the particular case where all possible different permutations must be found
    allUnique = true;
    
elseif (N_SHUFFLE==Max_shuffling)
	N_SHUFFLE = Max_shuffling;
    permutDec = [0:N_SHUFFLE-1];  % handle the particular case where all possible different permutations must be found
    allUnique = true;
    
% 2) Find a combination of N_SHUFFLE unique integer values
else
    permutDec = zeros(1,N_SHUFFLE)-1; % full of -1
    for shuf_ix = 1:N_SHUFFLE
        while(~allUnique)
            permutDec(shuf_ix) = randi(Max_shuffling-1)-1; 
            allUnique = length(unique(permutDec(permutDec>=0)))==shuf_ix;
        end    
        allUnique = false;
    end
end


% 3) Convert decimal values to logical permutation indexes
permutBin = dec2bin(permutDec);     % this is a char matrix
SCHUFFLE_INDEX = true(N_SHUFFLE,VEC_LENGTH); % this is a logical one
for rows_ix = 1:N_SHUFFLE
    for col_ix = 1:VEC_LENGTH
        SCHUFFLE_INDEX(rows_ix,col_ix) = logical(str2double(permutBin(rows_ix,col_ix)));
    end
end


% TRASH
% if(N_SHUFFLE <= Max_shuffling/2)
%     while(~allUnique)
%         permutDec = randi(Max_shuffling,[1 N_SHUFFLE])-1; 
%         allUnique = length(unique(permutDec))==N_SHUFFLE;
%     end
% else
%     while(~allUnique)
%         permutDec = randi(Max_shuffling,[1 Max_shuffling-N_SHUFFLE])-1; 
%         allUnique = length(unique(permutDec))==(Max_shuffling-N_SHUFFLE);
%     end
%     permutDec = setdiff(0:Max_shuffling-1, permutDec);
% end
