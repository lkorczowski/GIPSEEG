function [Y] = zMean(X, dim)
% function [Y] = zMean(X, dim)
% Average or mean value of array after Fisher transform (z-score).
% Fisher transform is a variance-stabilizing transformation necessary when 
% suming correlation coefficients.
%
% INPUTS
% - X       --> data vector or matrix 
% - dim     --> returns the mean values for elements along the dimension of
%               X specified by scalar dim.
%
% OUTPUTS
% - Y       --> array of mean values of X after Fisher transform. 
%
% HISTORY 
% First version: Jonas Chatel-Goldman @ GIPSA-Lab, 24/06/2013


if (nargin<2) || (isempty(dim))
   dim = 1; 
end
if isrow(X)&&(dim==1)
    dim = 2;
end
Y =  tanh(mean(atanh(X),dim));
