function [Y] = zStd(X, dim)
% function [Y] = zStd(X, dim)
% Standard deviation of array after Fisher transform (z-score).
% Fisher transform is a variance-stabilizing transformation necessary when 
% suming correlation coefficients.
%
% INPUTS
% - X       --> data vector or matrix 
% - dim     --> returns the mean values for elements along the dimension of
%               X specified by scalar dim.
%
% OUTPUTS
% - Y       --> array of std values of X after Fisher transform. 
%
% HISTORY 
% First version: Jonas Chatel-Goldman @ GIPSA-Lab, 06/08/2013


if (nargin<2) || (isempty(dim))
   dim = 1; 
end
if isrow(X)
    X = X';
end
Y =  tanh(std(atanh(X),[],dim));