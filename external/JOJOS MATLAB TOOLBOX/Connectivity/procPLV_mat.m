function [PLV] = procPLV_mat(DATA1, DATA2, WINDOW_SIZE, OVERLAP, PUT_NaN)
% [PLV] = procPLV_mat(DATA1, DATA2, WINDOW_SIZE, OVERLAP, PUT_NaN)
%
% This function is identical to 'procPLV' except that it computes PLV on
% multiple pairs of vectors at the same time (paired vectors are
% corresponding rows of DATA1 and DATA2)
%
% Inputs:
% - DATA1           --> 2D matrix with the data (P channels by N samples)
% - DATA2           --> 2D matrix with the data (P channels by N samples)
% - WINDOW_SIZE     --> (optionnal) size of the window in number of samples 
% - OVERLAP         --> (optionnal) percentage of the overlapping window
%                                   (from 0 to 0.99)
% - PUT_NaN         --> (optionnal) if set to 0, put zeros instead of NaNs
%                                   when PLV value is not available
%                                   (default is 1)
%
% Outputs:
% - PLV             --> matrix of PLV values (P channels by N samples)

%
% History
% Last version:  20/11/2012
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

% default overlapping
if(nargin < 5)  PUT_NaN = 1;  end
% default overlapping
if(nargin < 4)  OVERLAP = .5;  end
% default windows size
if(nargin < 3)  WINDOW_SIZE = 128;  end


s_data1 = size(DATA1);
s_data2 = size(DATA2);
nbChan  = s_data1(1);  % nb channels
if find(s_data1 ~= s_data2)
   error('[function procPLV] Input matrices must have the same size (P channels by N samples)'); 
end
if s_data1(2) < WINDOW_SIZE
    error('[function procPLV] Input data is shorter than specified window length!');
end

% Calculate the number of overlapped windows 
if(OVERLAP == 0)
    nbWin = floor(s_data1(2) / WINDOW_SIZE);
else
    nbFullWin = floor(s_data1(2)/WINDOW_SIZE); 
    nbWin = 1 + (nbFullWin-1)/(1-OVERLAP) + floor((s_data1(2)-((nbFullWin)*WINDOW_SIZE))/((1-OVERLAP)*WINDOW_SIZE));
end

% creation of PLV vector
PLV = zeros(nbChan,nbWin);

% Phase extraction using Hilbert transform
DATA_h1  = hilbert(DATA1');
DATA_h2  = hilbert(DATA2');
clear DATA1 DATA2

% perc10w = floor(s_data1(2)/5); % discard 5% of hilbert transform (spurious borders)
% [PHASE] = unwrap(angle(DATA_h(perc10w:end-perc10w,:)));
[PHASE1] = unwrap(angle(DATA_h1),[],2)';
[PHASE2] = unwrap(angle(DATA_h2),[],2)';
clear DATA_h1 DATA_h2

%% Loop on all time windows
for window_ix = 1 : nbWin
    
    % time markers to select the data
    t1 = floor((window_ix-1) * (1-OVERLAP) * WINDOW_SIZE) +1;   % marker of the beginning of the time window
    t2 = t1 + WINDOW_SIZE -1;                                   % marker of the end of the time window
    
    % select current window   
    cPHASE1 = PHASE1(:,t1:t2);  % format: (nbChan, WINDOW_SIZE)
    cPHASE2 = PHASE2(:,t1:t2);

    % compute relative phase and PLV value
    RP = cPHASE1-cPHASE2;       % relative phase
    PLV(:,window_ix)=abs(sum(exp(1i*RP),2)/WINDOW_SIZE);
end
clear PHASE1 PHASE2 RP cPHASE1 cPHASE2

%% --- Interpolation of PLV values 
% Interpolation using Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)
% This method is not time consuming and more importantly it has no overshoots
axis_interp	= 1:s_data1(2);
axis_rho 	= 1:s_data1(2)/nbWin:s_data1(2);
PLV         = pchip(axis_rho,PLV,axis_interp);

% Pad first half window with NaN or zeros, and crop last half window
HalfWin = floor(WINDOW_SIZE/2);
if PUT_NaN == 1         % put NaN values where measures are not available
	PLV = [ NaN(nbChan,HalfWin) PLV(:,1:end-HalfWin) ]; 
elseif PUT_NaN ==0      % put zero values where measures are not available
	PLV = [ zeros(nbChan,HalfWin) PLV(:,1:end-HalfWin) ]; 
end



return








% TRASH
% 
% % PLV values are centered in the middle of the sliding window
% PLV_interp  = resample(PLV',round(s_data1(2)/nbWin),1)';
% PLV_interp(PLV_interp>1) = 1 ; % threshold to '1', since interpolation sometimes makes spurious overrun
% HalfWin     = floor(WINDOW_SIZE/2);
% if PUT_NaN == 1
% 	PLV = NaN(nbChan,s_data1(2)); % put NaN values where PLV is not available
% elseif PUT_NaN ==0
% 	PLV = zeros(nbChan,s_data1(2)); % put zeros values where PLV is not available
% end
% PLV(:,HalfWin:length(PLV_interp)+HalfWin-1) = PLV_interp;
% PLV = PLV(:,1:s_data1(2));    
% 
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
