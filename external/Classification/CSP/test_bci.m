function y = test_bci(X,S,T,w,b)
% Prediction = test_bci(Raw-Block, Spatial-Flt, Temporal-Flt, Weights, Bias)
% A linear online detector for oscillatory processes.
%
% In:
%   X : incoming raw sample data [Samples x Channels]
%   S,T: spatio-temporal filter [Channels x Filters], [Samples x Samples]
%   w,b: linear classifier [Filters x 1], [1 x 1]
%
% Out:
%   y : the processed output
%	
% Notes:	
%	y can be post-processed by sign(y) for pure classification, or 1./(1+exp(-y)) for logistic regression
%
 
global B; % B is the buffer
if any(size(B) ~= [length(T),length(S)]) 
    B = zeros(length(T),length(S)); 
end
B = [B;X]; 
B = B(end-length(T)+1:end,:);
y = log(var(T*(B*S)))*w - b;