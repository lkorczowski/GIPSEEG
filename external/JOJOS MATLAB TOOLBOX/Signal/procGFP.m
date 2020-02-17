function [GFP] = procGFP(DATA)
% function  [GFP] = procGFP(DATA)
%
% ************************************************************************
% Calculate Global Field Power (GFP), Lehmann & Skrandies (1984, p. 235)
% (merely the standard deviation of the data)
%
% Input:
% ------
% DATA	--> input timeseries, Nchannel x Ntime
%
% Output:
% -------
% GFP   --> vector of Global Field Power (reference-free)
%
% History:
% --------
% Last version:  2013-02-08
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com



    % simply process std on centered data 
    GFP = std(DATA,1);

    % then average it: moving average with null phase shift
	GFP = filtfilt([1/3 1/3 1/3],1,GFP); 


end





%% test function
% nbSample = 1000;
% t_axis = 1:nbSample;
% x1 = sin(t_axis);
% x2 = .8*sin(t_axis+0.25);
% x3 = .2*sin(1.5*t_axis+0.4);
% x1(50:100) = 2*randn(1,100-50+1);
% x2([50:150,700:710]) = [2*randn(1,150-50+1) , 3*randn(1,710-700+1)];
% x3(200:400) = 3*randn(1,400-200+1);
% DATA = [x1 ; x2 ; x3];
% GFP = procGFP(DATA);
% figure,
% PlotEEG_maison([DATA ; GFP]),    






