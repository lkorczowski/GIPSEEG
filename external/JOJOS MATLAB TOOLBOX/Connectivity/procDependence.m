function [RHO_ORDINARY, RHO_INSTANT, RHO_LAGGED] = procDependence(DATA, S_DATA, LINEARITY, WINDOW_SIZE, OVERLAP, FREQ_RANGE, WINDOW_TYPE, PUT_NaN, UI_DISP)
% This function computes the squared linear dependence (coherence) and non
% linear PUT_NaN dependence (phase synchronisation) of two univariate/multivariate time series.
%
% 
% INPUTS
% - DATA            --> multivariate time serie, with concatenated datasets 
%                       (M channels by N time samples)
% - S_DATA         	--> (optional) size of each dataset (S must be a column vector)
%                       if not set each row is considered as a single dataset
% - LINEARITY       --> (optional) specifies if DFT is normalized (non-linear 
%                       case, value 0) or not (linear case, value 1)
% - WINDOW_SIZE    	--> (optional) size of the window in number of samples
% - OVERLAP         --> (optional) percentage of the overlapping window (from 0 to 0.99)
%                       all samples are considered uniformely only if OVERLAP has value 1 - 2^(-n)
% - FREQ_RANGE    	--> (optional) average output over specified frequency range (only one freq then)
%                       format must be [sample_rate min_freq_in_Hz max_freq_in_Hz]
% - WINDOW_TYPE   	--> (optional) string with the type of the window
% - PUT_NaN         --> (optionnal) if set to 0, padd with zeros instead of 
%                                   NaNs when dependance value is not available
%                                   (default is to pad with zeros)
% - UI_DISP         --> (optional) display waiting bar during execution (0 or 1)
%
% OUTPUTS
% - RHO_ORDINARY    --> square ordinary dependence      (NbFreq * N time samples)
% - RHO_INSTANT     --> square instantaneous coherence  (NbFreq * N time samples)
% - RHO_LAGGED      --> square lagged coherence         (NbFreq * N time samples)
%   
% REFERENCES
% RD Pascual-Marqui: Instantaneous and lagged measurements of linear and nonlinear dependence between groups of 
% multivariate time series: frequency decomposition. arXiv:0711.1455 [stat.ME], 2007-November-09, http://arxiv.org/abs/0711.1455 
%
% HISTORY 
% Last version: 14/02/2013
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com



%% --- Assessing inputs and various init
% default UI display
if(nargin < 9)                          	UI_DISP  = 0;           end
% default padding
if(nargin < 8) || (isempty(PUT_NaN))        PUT_NaN  = 1;            end
% default windows type for cospectra processing
if(nargin < 7) || (isempty(WINDOW_TYPE))    WINDOW_TYPE  = 'hanning';  end
 % default frequency range 
if(nargin < 6) || (isempty(FREQ_RANGE)) 	FREQ_RANGE   = [];      end
% default overlapping for cospectra processing
if(nargin < 5) || (isempty(OVERLAP))        OVERLAP     = .75;      end
% default windows size for cospectra processing
if(nargin < 4) || (isempty(WINDOW_SIZE))    WINDOW_SIZE  = 128; 	end
% default measure is linear
if(nargin < 3) || (isempty(LINEARITY))      LINEARITY   = 1;    	end
% default data size is to consider each row independently (univariate measure)
if(nargin < 2) || (isempty(S_DATA))         S_DATA   = ones(size(DATA,1),1); end

if size(S_DATA,2) > 1    S_DATA = S_DATA';     end  
cix         = [cumsum(S_DATA)-S_DATA+1, cumsum(S_DATA)]; % indexes of each dataset
nbChan      = size(DATA,1);
nbSamples   = size(DATA,2);
nbDatasets  = length(S_DATA);

% Adjust the window size to the next power of 2 (for FFT calculation)
WINDOW_SIZE = 2 ^ nextpow2(WINDOW_SIZE);

% Calculate the number of frequencies and time windows
ZP_factor = 64;     % factor for fft zero-padding
nbFreqs = (WINDOW_SIZE / 2) * ZP_factor;
nbWin = floor(nbSamples / WINDOW_SIZE);
if(UI_DISP) 
    progress = 0; 
end

   
% Memory pre-allocation 
RHO_ORDINARY    = zeros(1, nbFreqs);  
RHO_INSTANT     = zeros(1, nbFreqs);
RHO_LAGGED      = zeros(1, nbFreqs);
S               = zeros(nbChan, nbChan, nbFreqs);

%% --- Create window from the inputs 
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
    case 'rect'
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


%% --- Loop on all time windows
for win_ix = 1 : nbWin
    % UI display
    if(UI_DISP)   
        progress2 = 10*fix(10*win_ix/nbWin);
        if progress2 ~= progress || win_ix==1
            progress = progress2;
            fprintf('%s','')
            disp(['Processing dependance, progression ' int2str(progress) '%']); 
        end
    end
    
    % various initializations
    t1 = floor((win_ix-1) * (1-OVERLAP) * WINDOW_SIZE) +1;   % marker of the beginning of the time window
    t2 = t1 + WINDOW_SIZE -1;                           % marker of the end of the time window
    

    % select current window and apodize it   
    cdata = DATA(:,t1:t2)' .* (win*ones(1,nbChan));
%         fdata = fft(cdata) ./ nbWin_averaged; % FFT calculation (complex data)
    fdata = fft(cdata,ZP_factor*WINDOW_SIZE) ./ nbWin; % FFT calculation with zero-padding
%     fdata = reshape(fdata(1:ZP_factor*WINDOW_SIZE/2,:)',[nbChan ZP_factor WINDOW_SIZE/2]);
%     fdata = squeeze(mean(fdata,2))';

    % cross-spectra calculation after normalization for non-linear case
    for f_ix = 1 : nbFreqs   
        if LINEARITY == 0
            for d_ix = 1 : nbDatasets   % for each dataset
                fd_sel = fdata(f_ix, cix(d_ix,1):cix(d_ix,2));
                fdata(f_ix,cix(d_ix,1):cix(d_ix,2)) = fd_sel*((fd_sel*fd_sel')^(-1/2)); % factoring out amplitude
            end
        end
        % cross-spectra calculation and sommation
        S(:,:,f_ix) = S(:,:,f_ix) + (fdata(f_ix,:)' * fdata(f_ix,:)) ;  
    end
end

% -- compute measures of dependence for all the frequencies 
for f_ix = 1 : nbFreqs   

    S_tot = S(:,:,f_ix); % index win_ix is for averaging over time

    % construct S_part, eg. for case 2/2: S_part = (S11(f),0 ; 0,S22(f))
    S_part = S_tot(cix(1,1):cix(1,2),cix(1,1):cix(1,2));
    for d_ix = 2 : nbDatasets        
        S_part  = blkdiag(S_part,S_tot(cix(d_ix,1):cix(d_ix,2),cix(d_ix,1):cix(d_ix,2)));   
    end
    % process determinants
    det_S_tot       = det(S_tot);    % ??? is it normal to have some determinants with complex values? 
    det_S_part      = det(S_part);
    det_Sreal_tot   = det(real(S_tot));
    det_Sreal_part  = det(real(S_part));
    % (test regularization)
%         det_S_tot       = det(S_tot + eye(nbChan));    % ??? is it normal to have some determinants with complex values? 
%         det_S_part      = det(S_part + eye(nbChan));
%         det_Sreal_tot   = det(real(S_tot + eye(nbChan)));
%         det_Sreal_part  = det(real(S_part + eye(nbChan)));
    % compute dependence measures
    RHO_ORDINARY(f_ix)  = real( 1 - (det_S_tot/det_S_part) );                                  % Ordinary dependence (eq.50)
    RHO_INSTANT(f_ix)   = real( 1 - (det_Sreal_tot/det_Sreal_part) );                          % Instantaneous dependence (eq.51)
    RHO_LAGGED(f_ix)    = real( 1 - (det_S_tot/det_S_part)/(det_Sreal_tot/det_Sreal_part) );   % Lagged dependence (eq. 52)
end


end






