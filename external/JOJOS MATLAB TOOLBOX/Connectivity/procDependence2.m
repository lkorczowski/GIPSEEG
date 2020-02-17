function [RHO_ORDINARY, RHO_INSTANT, RHO_LAGGED] = procDependence2(DATA, S_DATA, LINEARITY, AVERAGING, WINDOW_SIZE, OVERLAP, FREQ_RANGE, WINDOW_TYPE, PUT_NaN, UI_DISP)
% This function computes the squared linear dependence (coherence) and non
% linear PUT_NaN
% dependence (phase synchronisation) of two univariate/multivariate time series.
%
% As compared to function 'procDependence', here output is not averaged
% over time windows. It is interpolated to obtain time a varying course 
% with as many time samples as input data.
% 
%
% INPUTS
% - DATA            --> multivariate time serie, with concatenated datasets 
%                       (M channels by N time samples)
% - S_DATA         	--> (optional) size of each dataset (S must be a column vector)
%                       if not set each row is considered as a single dataset
% - LINEARITY       --> (optional) specifies if DFT is normalized (non-linear 
%                       case, value 0) or not (linear case, value 1)
% - AVERAGING       --> (optional) vector with 2 positive integers [av_t av_f],
%                       specifies number of time windows 'av_t' and/or frequency bins 'av_f' 
%                       on which cospectra are averaged for a better estimation, i.e.
%                       output frequency resolution is: (av_f * sample_rate) / WINDOW_SIZE
%                       output time resolution is:      WINDOW_SIZE / (av_t * sample_rate)
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
% Last version:  14/02/2013
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


%% --- Assessing inputs and various init
% default UI display
if(nargin < 10)                          	UI_DISP  = 0;           end
% default padding
if(nargin < 9) || (isempty(PUT_NaN))        PUT_NaN  = 1;            end
% default windows type for cospectra processing
if(nargin < 8) || (isempty(WINDOW_TYPE))    WINDOW_TYPE  = 'hanning';  end
 % default frequency range 
if(nargin < 7) || (isempty(FREQ_RANGE)) 	FREQ_RANGE   = [];      end
% default overlapping for cospectra processing
if(nargin < 6) || (isempty(OVERLAP))        OVERLAP     = .75;      end
% default windows size for cospectra processing
if(nargin < 5) || (isempty(WINDOW_SIZE))    WINDOW_SIZE  = 128; 	end
% default averaging is on 4 time windows
if(nargin < 4) || (isempty(AVERAGING))      AVERAGING   = [4 1];    	end
% default measure is linear
if(nargin < 3) || (isempty(LINEARITY))      LINEARITY   = 1;    	end
% default data size is to consider each row independently (univariate measure)
if(nargin < 2) || (isempty(S_DATA))         S_DATA   = ones(size(DATA,1),1); end

if size(S_DATA,2) > 1    S_DATA = S_DATA';     end  
if any(mod(AVERAGING,1)) || any(AVERAGING<=0) || (length(AVERAGING)~=2)
    error('Input AVERAGING must be a 2-element positive integer vector.'); 
end
cix         = [cumsum(S_DATA)-S_DATA+1, cumsum(S_DATA)]; % indexes of each dataset
nbChan      = size(DATA,1);
nbSamples   = size(DATA,2);
nbDatasets  = length(S_DATA);

% Adjust the window size to the next power of 2 (for FFT calculation)
WINDOW_SIZE = 2 ^ nextpow2(WINDOW_SIZE);

% Calculate the number of frequencies
nbFreqs = (WINDOW_SIZE / 2);
nbFreqs_averaged = floor(nbFreqs/AVERAGING(2)); % number of freq is decreased when averaging over frequency
ZP_factor = 64;     % factor for fft zero-padding

% Calculate the number of main time windows 
winSize_averaged = WINDOW_SIZE*AVERAGING(1); % window size is increased when averaging over time
nbWin_averaged = floor(nbSamples / winSize_averaged);
if nbSamples < winSize_averaged
    error('Cospectra computation: data segment too short for specified window length');
end

if(UI_DISP) progress = 0; end

    
% Memory pre-allocation 
RHO_ORDINARY    = zeros(nbWin_averaged, nbFreqs_averaged);  
RHO_INSTANT     = zeros(nbWin_averaged, nbFreqs_averaged);
RHO_LAGGED      = zeros(nbWin_averaged, nbFreqs_averaged);


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
for main_win_ix = 1 : nbWin_averaged
    % UI display
    if(UI_DISP)   
        progress2 = 10*fix(10*main_win_ix/nbWin_averaged);
        if progress2 ~= progress || main_win_ix==1
            progress = progress2;
            disp(['Processing dependance, progression ' int2str(progress) '%']); 
        end
    end
    
    % various initializations
%     S = zeros(nbChan, nbChan, nbFreqs);     % cross-spectal matrix must be reset here
    S = zeros(nbChan, nbChan,WINDOW_SIZE);
    t_average_beg = (main_win_ix-1)*winSize_averaged + 1; 	% begin sample of main time window    
    t_average_end = t_average_beg + winSize_averaged - 1;   % end sample of main time window 
    t1 = t_average_beg;                                     % begin sample of current time window
    t2 = t1 + WINDOW_SIZE -1;                               % end sample of current time window
    win_ix = 1; 
    
    % estimate cospectra while current window is at least in part comprised within main window
    while (t1 < t_average_end) && (t2 <= nbSamples)
        % select current window and apodize it   
        cdata = DATA(:,t1:t2)' .* (win*ones(1,nbChan));
%         fdata = fft(cdata) ./ nbWin_averaged; % FFT calculation (complex data)
        fdata = fft(cdata,ZP_factor*WINDOW_SIZE) ./ nbWin_averaged; % FFT calculation with zero-padding
        fdata = reshape(fdata(1:ZP_factor*WINDOW_SIZE/2,:)',[nbChan ZP_factor WINDOW_SIZE/2]);
        fdata = squeeze(mean(fdata,2))';
        
        % cross-spectra calculation after normalization for non-linear case
        for f_ix = 1 : size(fdata,1)%nbFreqs   
            if LINEARITY == 0
                for d_ix = 1 : nbDatasets   % for each dataset
                    fd_sel = fdata(f_ix, cix(d_ix,1):cix(d_ix,2));
                    fdata(f_ix,cix(d_ix,1):cix(d_ix,2)) = fd_sel*((fd_sel*fd_sel')^(-1/2)); % factoring out amplitude
                end
            end
            % cross-spectra calculation and sommation
            S(:,:,f_ix) = S(:,:,f_ix) + (fdata(f_ix,:)' * fdata(f_ix,:)) ;  
        end

        % calculate indexes of next time window (overlapping is performed WITHIN main window)
        t1 = t1 + (1-OVERLAP)*WINDOW_SIZE;	% begin sample of current time window
        t2 = t1 + WINDOW_SIZE -1;          	% end sample of current time window
        win_ix = win_ix+1;
    end
    
    
    % -- compute measures of dependence for all the frequencies 
    for f_ix = 1 : nbFreqs_averaged   
        
        % average cross-spectra as required
        freq_beg = (f_ix-1)*AVERAGING(2) + 1;       % first index of current frequency bin 
        freq_end = freq_beg + AVERAGING(2) - 1;     % last index of current frequency bin 
     	S_tot = mean( S(:,:,freq_beg:freq_end)./win_ix, 3); % index win_ix is for averaging over time
        
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
        RHO_ORDINARY(main_win_ix,f_ix)  = real( 1 - (det_S_tot/det_S_part) );                                  % Ordinary dependence (eq.50)
        RHO_INSTANT(main_win_ix,f_ix)   = real( 1 - (det_Sreal_tot/det_Sreal_part) );                          % Instantaneous dependence (eq.51)
        RHO_LAGGED(main_win_ix,f_ix)    = real( 1 - (det_S_tot/det_S_part)/(det_Sreal_tot/det_Sreal_part) );   % Lagged dependence (eq. 52)
    end
end

%% --- Average output over specified frequency range
if ~isempty(FREQ_RANGE)
    freqRes = (.5*FREQ_RANGE(1))/nbFreqs_averaged;
    f_bins  = .5*freqRes : freqRes : .5*(FREQ_RANGE(1)-freqRes);  % these are the central freq of bins considered here
    f_min   = FREQ_RANGE(2) - .5*freqRes; 
    f_max   = FREQ_RANGE(3) + .5*freqRes;
    f_binsToAverage = f_bins>f_min & f_bins<f_max;
	RHO_ORDINARY    = mean(RHO_ORDINARY(:,f_binsToAverage),2);
    RHO_INSTANT     = mean(RHO_INSTANT(:,f_binsToAverage),2);
    RHO_LAGGED      = mean(RHO_LAGGED(:,f_binsToAverage),2);
    nbFreqs_averaged = 1; % now outputs only one averaged frequency 
end


%% --- Interpolation of dependence values 
% Interpolation using Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)
% This method is very low time consuming and more importantly it makes no overshoots
axis_interp	= 1:nbSamples;
axis_rho 	= 1:nbSamples/nbWin_averaged:nbSamples;
RHO_ORDINARY    = pchip(axis_rho,RHO_ORDINARY',axis_interp);
RHO_INSTANT     = pchip(axis_rho,RHO_INSTANT',axis_interp);
RHO_LAGGED      = pchip(axis_rho,RHO_LAGGED',axis_interp);

% Pad first half window with NaN or zeros, and crop last half window
HalfWin         = floor(winSize_averaged/2);
if PUT_NaN == 1         % put NaN values where measures are not available
	RHO_ORDINARY    = [ NaN(nbFreqs_averaged,HalfWin) RHO_ORDINARY(:,1:end-HalfWin) ]; 
    RHO_INSTANT     = [ NaN(nbFreqs_averaged,HalfWin) RHO_INSTANT(:,1:end-HalfWin) ]; 
    RHO_LAGGED      = [ NaN(nbFreqs_averaged,HalfWin) RHO_LAGGED(:,1:end-HalfWin) ]; 
elseif PUT_NaN ==0      % put zero values where measures are not available
	RHO_ORDINARY    = [ zeros(nbFreqs_averaged,HalfWin) RHO_ORDINARY(:,1:end-HalfWin) ]; 
    RHO_INSTANT     = [ zeros(nbFreqs_averaged,HalfWin) RHO_INSTANT(:,1:end-HalfWin) ]; 
    RHO_LAGGED      = [ zeros(nbFreqs_averaged,HalfWin) RHO_LAGGED(:,1:end-HalfWin) ]; 
end

    
end






