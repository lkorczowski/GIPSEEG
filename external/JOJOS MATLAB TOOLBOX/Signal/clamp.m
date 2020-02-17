function [y] = clamp(x, range)
% function [x] = clamp(x)
% ************************************************************************
% Restricts a value to be within a specified range defined by values min
% and max.
%
% INPUTS
% - x       --> scalar, vector or matrix 
% - range   --> [min max] specifies restricting value range
%
% OUTPUTS
% - y       --> "clamped" output
%
% HISTORY 
% First version: 07/02/2013
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

    y = x;
    y(x<range(1)) = range(1);
    y(x>range(2)) = range(2);

end