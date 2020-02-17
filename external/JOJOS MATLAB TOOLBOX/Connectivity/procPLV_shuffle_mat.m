function [PLV_obs PLV_perm] = procPLV_shuffle_mat(DATA1, DATA2, WINDOW_SIZE, OVERLAP, PERMS, PUT_NaN)
% [PLV_obs PLV_perm] = procPLV_shuffle_mat(DATA1, DATA2, WINDOW_SIZE, OVERLAP, PERMS, PUT_NaN)
%
% This function is identical to 'procPLV_mat' except that: 
% 1) it computes PLV either on corresponding time windows (observed set) or 
% on uniquely shuffled windows NBPERM times (permutation set).
% 2) PLV values are no longer interpolated.
%
% Inputs:
% - DATA1           --> 2D matrix with the data (P channels, N samples)
% - DATA2           --> 2D matrix with the data (P channels, N samples)
% - WINDOW_SIZE     --> (optionnal) size of the window in number of samples 
% - OVERLAP         --> (optionnal) percentage of the overlapping window
%                                   (from 0 to 0.99)
% - PERMS           --> (optionnal) can be either a scalar with number of 
%                                   unique permutations for shuffling of 
%                                   time windows; or a matrix of size 
%                                   (NBPERM,nbWin) with unique rows filled 
%                                   with random indexes in range [1->nbWin]
% - PUT_NaN         --> (optionnal) if set to 0, put zeros instead of NaNs
%                                   when PLV value is not available
%                                   (default is 1)
%
% Outputs:
% - PLV_obs         --> 2D matrix of PLV values computed on corresponding 
%                       time windows (P, nbWin)
% - PLV_perm        --> 3D matrix of PLV values computed on shuffled time 
%                       windows (P, nbWin, NBPERM)
%
% History
% First version: Jonas Chatel-Goldman @ GIPSA-Lab, 10/01/2013

% default NaN or 0 vals
if(nargin < 6)  PUT_NaN = 1;  end
% default number of permutations
if(nargin < 5)  NBPERM = 1000;  end
% default overlapping
if(nargin < 4)  OVERLAP = 0;  end
% default windows size
if(nargin < 3)  WINDOW_SIZE = 128;  end


% assess inputs
s_data1 = size(DATA1);
s_data2 = size(DATA1);
nbChan  = s_data1(1);  % nb channels
if find(s_data1 ~= s_data2)
   error('[function procPLV] Input matrices must have the same size (P channels by N samples)'); 
end
if s_data1(2) < WINDOW_SIZE
    error('[function procPLV] Input data is shorter than specified window length!');
end
if(exist('PERMS','var') && length(PERMS)==1)
    NBPERM = PERMS;
end

% Calculate the number of overlapped windows 
if(OVERLAP == 0)
    nbWin = floor(s_data1(2) / WINDOW_SIZE);
else
    nbFullWin = floor(s_data1(2)/WINDOW_SIZE); 
    nbWin = 1 + floor((nbFullWin-1)/(1-OVERLAP) + (s_data1(2)-((nbFullWin)*WINDOW_SIZE))/((1-OVERLAP)*WINDOW_SIZE));
end

% Check input variable 'PERMS'
if(exist('PERMS','var') && length(PERMS)>1)
   if(size(PERMS,2)~=nbWin)
       error('Input index matrix ''PERMS'' must have as many columns as time windows ''nbWin''');
   end
   doComputeShuffle = 0;
   NBPERM = size(PERMS,1);
   SHUFFLE_INDEX = PERMS;
else
   doComputeShuffle = 1;
end

% Allocation of PLV matrices
PLV_obs     = zeros(nbChan,nbWin);
PLV_perm    = zeros(nbChan,nbWin,NBPERM);

% Compute NBPERM unique shuffling indexes for time windows permutations
if doComputeShuffle
    SHUFFLE_INDEX = compute_shuffle_ix();
end

% Phase extraction using Hilbert transform
DATA_h1  = hilbert(DATA1');
DATA_h2  = hilbert(DATA2');

% perc10w = floor(s_data1(2)/5); % discard 5% of hilbert transform (spurious borders)
% [PHASE] = unwrap(angle(DATA_h(perc10w:end-perc10w,:)));
[PHASE1] = unwrap(angle(DATA_h1),[],2)';
[PHASE2] = unwrap(angle(DATA_h2),[],2)';


%% Loop on all time windows
for window_ix = 1 : nbWin
    if(mod(window_ix,100) == 0)
        fprintf('%s',' ');
        disp(['Processing PLV on ' int2str(NBPERM) ' permutations. Progress: ' int2str(100*window_ix/nbWin) '%...']);
    end
    % time markers to select the data
    t1_ref  = floor((window_ix-1) * (1-OVERLAP) * WINDOW_SIZE) +1;   % begin sample of reference window
    t2_ref  = t1_ref + WINDOW_SIZE -1;                                   % end sample of reference window
    
    % select current window   
    cPHASE1 = PHASE1(:,t1_ref:t2_ref);  % format: (nbChan, WINDOW_SIZE)
    cPHASE2 = PHASE2(:,t1_ref:t2_ref);  % corresponding window

    % compute observed PLV value (corresponding windows)
    RP = cPHASE1-cPHASE2;       % relative phase
    PLV_obs(:,window_ix)=abs(sum(exp(1i*RP),2)/WINDOW_SIZE);    
    
    % compute permuted PLV value (shuffled windows)
    for perm_ix = 1:NBPERM
        permWin_ix = SHUFFLE_INDEX(perm_ix,window_ix);
        t1_shuf = floor((permWin_ix-1) * (1-OVERLAP) * WINDOW_SIZE) +1;  % begin sample of shuffled window
        t2_shuf = t1_shuf + WINDOW_SIZE -1;                             % end sample of shuffled window    
        cPHASE2 = PHASE2(:,t1_shuf:t2_shuf);
        RP = cPHASE1-cPHASE2;       % relative phase
        PLV_perm(:,window_ix,perm_ix)=abs(sum(exp(1i*RP),2)/WINDOW_SIZE); 
    end
end   
return



function [Shuf_index] = compute_shuffle_ix
    if(NBPERM > (2^nbWin))
        error(['Requested number of permutations exceed max possible perm with ' int2str(nbWin) ' time windows']);
    end
    Shuf_index      = zeros(NBPERM,nbWin);
    Shuf_index(1,:) = randperm(nbWin); 
    for p_ix = 2:NBPERM
%         disp(['progress: ' int2str(100*p_ix/NBPERM) '%'])
        allUnique = false;
        while(~allUnique) 
            % create random vector with unique values
            Shuf_index(p_ix,:) = randperm(nbWin); 
            
            % reiterate until this vector is unique within Shuf_index
            sameShuff_prev = false(1,nbWin);
            sameShuff_prev(1:p_ix-1) = true;
            for win_ix = 1:nbWin
                sameShuff_curr = (Shuf_index(sameShuff_prev,win_ix) == Shuf_index(p_ix,win_ix));
                if ~any(sameShuff_curr)
                    allUnique = true;
                    break;
                end
                sameShuff_prev = sameShuff_curr;
            end
        end
    end
end

end




% TRASH
%     t_plv_1 = t1 + WINDOW_SIZE/2;  % PLV values are centered in the middle of sliding window
%     t_plv_2 = t_plv_1 + (1-OVERLAP)*WINDOW_SIZE;
%     PLV(t_plv_1:t_plv_2)=abs(sum(exp(1i*RP))/length(RP));
% PLV(1:floor(WINDOW_SIZE/2):end) = PLV_interp(1:s_data1(2)-floor(WINDOW_SIZE/2)+1);
% 
% if nbWin < 9
%     PLV_interp = Hold(PLV,ceil(s_data(2)/nbWin));  
% function y=Hold(x,r) 
% % SigEch        ->  line data vector
% % r             ->  nb sample to hold
%     y = ones(r,1)*x; 
%     y = reshape(y,1,size(y,1)*size(y,2));
% return

% function [Shuf_index] = compute_shuffle_ix
%     Max_shuffling   = 2^(nbWin); 
%     Shuf_index      = zeros(NBPERM,nbWin);
%     Shuf_index(1,:) = randperm(nbWin); 
%     for p_ix = 2:NBPERM
%         disp(['progress: ' int2str(100*p_ix/NBPERM) '%'])
%         allUnique = false;
%         while(~allUnique) 
%             Shuf_index(p_ix,:) = randperm(nbWin); 
%             for p_ix2 = 1:p_ix-1
%                 allUnique = (length(find(Shuf_index(p_ix,:)==Shuf_index(p_ix2,:))) ~= nbWin); % false only when all values are identical 
%                 if(~allUnique)  
%                     break;  
%                 end
%             end
%         end
%     end
% end

