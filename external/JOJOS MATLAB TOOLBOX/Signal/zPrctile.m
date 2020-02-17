function [Y] = zPrctile(X, p, dim)
% function [Y] = zPrctile(X, p, dim)
% percentiles of the values in X after Fisher transform (z-score).
% Fisher transform is a variance-stabilizing transformation necessary when 
% suming correlation coefficients.
%
% INPUTS
% - X       --> data vector or matrix 
% - p       --> scalar or a vector of percent values (from 0 to 100)
% - dim     --> (optional) calculates percentiles along dimension dim. 
%               The dim'th dimension of Y has length length(p;
%
% OUTPUTS
% - Y       --> array of percentile values of X after Fisher transform. 
%
% HISTORY 
% First version: Jonas Chatel-Goldman @ GIPSA-Lab, 24/06/2013


if (nargin<3) || (isempty(dim))
   dim = 1; 
end

Y = tanh(prctile(atanh(X),p,dim));

