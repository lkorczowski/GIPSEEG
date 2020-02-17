function Conv=ConvergenceLatencies(Lat1,Lat2)
% Conv=ConvergenceLatencies(Lat1,Lat2)
%   Converge criterion between 2 iterations for the Latency estimation algorithm
% *** History: 2015-03-19
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
% see also: LatencyEstimation
 Conv=mean(abs(Lat1-Lat2));
end