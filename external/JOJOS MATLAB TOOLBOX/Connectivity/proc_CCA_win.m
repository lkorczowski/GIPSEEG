function [A,B,R,stats,t,masked] = proc_CCA_win(DATA, SIZE_SETS, WINDOW_SIZE, WINDOW_TYPE, OVERLAP, MASK)
% [A,B,r,stats,t,masked] = proc_CCA_win(DATA, SIZE_SETS, WINDOW_SIZE, WINDOW_TYPE, OVERLAP, MASK)
%
% This function performs Canonical Correlation Analysis (CCA) on two 
% datasets sliced in a number of time windows. One CCA is computed for each 
% of these time window. 
%
% *** Inputs *** 
% - DATA            --> matrix with concatenated datasets X and Y [N samples x d1+d2 channels]
% - SIZE_SETS       --> vector with respective sizes of X and Y dataset: (d1,d2)
% - WINDOW_SIZE     --> (optionnal) size of the window in number of samples
% - WINDOW_TYPE     --> (optionnal) string with the type of the window
% - OVERLAP         --> (optionnal) percentage of the overlapping window (from 0 to 0.99)
% - MASK            --> (optionnal) binary mask [1 x DATA_LENGTH] 
%                                   -> epochs with at least one sample set
%                                      to 0 are not considered
%
%  *** Outputs *** 
% - A,B         --> sample canonical coefficients: [d1 x d x w] and [d2 x d x w]
%                   where d = min(rank(X),rank(Y)) and w = number of windows
% - R           --> vectors of sample canonical correlations: [w x d] 
% - stats       --> statistical info (see "canoncorr" function documentation)
% - t           --> time vector containing indexes of each window middle
%                   sample: [1 x w] 
% - masked      --> binary vector with same size of 't' containing 0 when
%                   period is masked, 1 otherwise
%
%  *** History *** 
% First version: 12/07/2012
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com



% default mask
if(nargin < 6)||(isempty(MASK))                 MASK  = true(1,size(DATA,1));  end
% default overlapping
if(nargin < 5)||(isempty(OVERLAP))              OVERLAP = .5;  end
% default windows type
if(nargin < 4)||(isempty(WINDOW_TYPE))          WINDOW_TYPE = 'rectangular';  end
% default windows size
if(nargin < 3)||(isempty(WINDOW_SIZE))          WINDOW_SIZE = 128;  end
% check size sets
if(nargin < 2)||(isempty(SIZE_SETS))            error('Size of each dataset must be specified');  end

UIdisplay = 0;


%% Adjust window from the inputs 
s_data = size(DATA);

switch lower(WINDOW_TYPE)
    case 'hamming'
        win = hamming(WINDOW_SIZE);
    case 'hann'
        win = hann(WINDOW_SIZE);
    case 'hanning'
        win = hanning(WINDOW_SIZE);
    case 'blackman'
        win = blackman(WINDOW_SIZE);
    case 'barthannwin'
        win = barthannwin(WINDOW_SIZE);
    case 'blackmanharris'
        win = blackmanharris(WINDOW_SIZE);
    case 'bohmanwin'
        win = bohmanwin(WINDOW_SIZE);
    case 'chebwin'
        win = chebwin(WINDOW_SIZE);
    case 'gausswin'
        win = gausswin(WINDOW_SIZE);
    case 'kaiser'
        win = kaiser(WINDOW_SIZE);
    case 'nuttallwin'
        win = nuttallwin(WINDOW_SIZE);
    case 'parzenwin'
        win = parzenwin(WINDOW_SIZE);
    case 'tukeywin'
        win = tukeywin(WINDOW_SIZE);
    case 'rectangular'
        win = ones(WINDOW_SIZE,1);
    case 'flattopwin'
        win = flattopwin(WINDOW_SIZE);
    case 'welch'
         %win = welchwin(WINDOW_SIZE, 0);
        Nd2=(WINDOW_SIZE-1)/2;
        for i=1 : WINDOW_SIZE
            win(i)=1-((i-Nd2)/Nd2)^2;
        end
        win = win';
    otherwise
        disp('WARNING - Unknown window type!');
        disp('Switching to default window type : Hamming Window');
        win = hamming(WINDOW_SIZE);
end

if s_data(1) < WINDOW_SIZE
    error('Input data too short for specified window length! You may have to transpose input data or change windows size...');
end


%% Calculate the number of overlapped windows and pre-allocate memory
if(OVERLAP == 0)
    number_windows = floor(s_data(1) / WINDOW_SIZE);
else
    nbFullWin = floor(s_data(1)/WINDOW_SIZE); 
    number_windows = 1 + (nbFullWin-1)/(1-OVERLAP) + floor((s_data(1)-((nbFullWin)*WINDOW_SIZE))/((1-OVERLAP)*WINDOW_SIZE));
end

% pre-allocation of memory 
d = min(SIZE_SETS);
A = zeros(d,SIZE_SETS(1),number_windows);
B = zeros(d,SIZE_SETS(2),number_windows);
R = zeros(number_windows,d); 
t = zeros(1,number_windows);
masked = zeros(1,number_windows);


%% Loop on all time windows
for window_ix = 1 : number_windows
    % UI display
    if(UIdisplay) && ~mod(window_ix, 50)  
        disp(['Processing CCA on time windows (' num2str(100*window_ix/number_windows,2) '%)']); 
    end
    
    % time markers to select the data, and time vector update
    t1 = floor((window_ix-1) * (1-OVERLAP) * WINDOW_SIZE) +1;   % marker of the beginning of the time window
    t2 = t1 + WINDOW_SIZE -1;                           % marker of the end of the time window
    t(window_ix) = fix(t1+WINDOW_SIZE/2);
    
    % if at least one mask sample is set to zero, do not consider this epoch
    if MASK(t1:t2) 
        % Selection of current window and apodization   
        cdata = DATA(t1:t2, :) .* (win*ones(1,s_data(2)));        
        % CCA calculation
        [A(:,:,window_ix),B(:,:,window_ix),R(window_ix,:),temp1,temp2,stats{window_ix}] = ...
                canoncorr(cdata(:,1:SIZE_SETS(1)),cdata(:,SIZE_SETS(2)+1:end)); 
        % update output 'masked' vector
        masked(window_ix) = 1; 
    end
end

masked = logical(masked);       % enables fast subsequent indexing
if ~(masked)
   error('[proc_CCA_win] WARNING : nothing computed because all windows are masked (at least in part)') 
end





