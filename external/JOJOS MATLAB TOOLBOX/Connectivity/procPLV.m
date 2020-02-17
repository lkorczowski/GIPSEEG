function [PLV] = procPLV(DATA, WINDOW_SIZE, OVERLAP, PUT_NaN)
% [PLV] = procPLV(DATA, WINDOW_SIZE, OVERLAP)
%
% This function computes Phase Locking Value (PLV) between two time series 
% on sliding windows using Hilbert transform. It outputs a vector of PLV 
% values interpolated so as to obtain as many time samples as input data.
% Input data must first be band-pass filtered around freq of interest.
%
% Caution: interpolation at the end might not be optimal, and border
% effects due to hilbert transform are not discared here. 
%
% Inputs:
% - DATA            --> 2D matrix with the data (2 channels by N samples)
% - WINDOW_SIZE     --> (optionnal) size of the window in number of samples 
% - OVERLAP         --> (optionnal) percentage of the overlapping window
%                                   (from 0 to 0.99)
% - PUT_NaN         --> (optionnal) if set to 0, padd with zeros instead of 
%                                   NaNs when PLV value is not available
%                                   (default is 1)
%
% Outputs:
% - PLV             --> vector of PLV values (N samples)

%
% History
% Last version:  20/11/2012
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


% default padding
if(nargin < 4)  PUT_NaN = 1;  end
% default overlapping
if(nargin < 3)  OVERLAP = .5;  end
% default windows size
if(nargin < 2)  WINDOW_SIZE = 128;  end


s_data = size(DATA);
if(s_data(1)~=2)
   error('[function procPLV] Input data should be in format (2 channels by N samples)'); 
end
if s_data(2) < WINDOW_SIZE
    error('[function procPLV] Input data is shorter than specified window length!');
end

% Calculate the number of overlapped windows 
if(OVERLAP == 0)
    number_windows = floor(s_data(2) / WINDOW_SIZE);
else
    nbFullWin = floor(s_data(2)/WINDOW_SIZE); 
    number_windows = 1 + (nbFullWin-1)/(1-OVERLAP) + floor((s_data(2)-((nbFullWin)*WINDOW_SIZE))/((1-OVERLAP)*WINDOW_SIZE));
end

% creation of PLV vector
% PLV = zeros(1,s_data(2));
PLV = zeros(1,number_windows);

% Phase extraction using Hilbert transform
DATA_h  = hilbert(DATA');
% perc10w = floor(s_data(2)/5); % discard 5% of hilbert transform (spurious borders)
% [PHASE] = unwrap(angle(DATA_h(perc10w:end-perc10w,:)));
[PHASE] = unwrap(angle(DATA_h'));


%% Loop on all time windows
for window_ix = 1 : number_windows
    
    % time markers to select the data
    t1 = floor((window_ix-1) * (1-OVERLAP) * WINDOW_SIZE) +1;   % marker of the beginning of the time window
    t2 = t1 + WINDOW_SIZE -1;                                   % marker of the end of the time window
    
    % select current window   
    cPHASE = PHASE(:,t1:t2);

    % compute relative phase and PLV value
    RP = cPHASE(1,:)-cPHASE(2,:);       % relative phase
    PLV(window_ix)=abs(sum(exp(1i*RP))/length(RP));
end


%% --- Interpolation of PLV values 
% Interpolation using Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)
% This method is not time consuming and more importantly it has no overshoots
axis_interp	= 1:s_data(2);
axis_rho 	= 1:s_data(2)/number_windows:s_data(2);
PLV         = pchip(axis_rho,PLV,axis_interp);

% Pad first half window with NaN or zeros, and crop last half window
HalfWin = floor(WINDOW_SIZE/2);
if PUT_NaN == 1         % put NaN values where measures are not available
	PLV = [ NaN(1,HalfWin) PLV(:,1:end-HalfWin) ]; 
elseif PUT_NaN ==0      % put zero values where measures are not available
	PLV = [ zeros(1,HalfWin) PLV(:,1:end-HalfWin) ]; 
end


return



% TRASH
% 
% if number_windows < 9
%     PLV_interp = Hold(PLV,ceil(s_data(2)/number_windows));  
% function y=Hold(x,r) 
% % SigEch        ->  line data vector
% % r             ->  nb sample to hold
%     y = ones(r,1)*x; 
%     y = reshape(y,1,size(y,1)*size(y,2));
% return


% % Test
% x   = [0:10]; 
% y   = sin(x); 
% xi  = [0:.25:10]; 
% yi  = interp1(x,y,xi,'spline'); 
% figure,
% plot(x,y,'o',xi,yi)
