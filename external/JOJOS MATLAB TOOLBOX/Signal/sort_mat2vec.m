function [A_sort,A_ix] = sort_mat2vec(A,mode)
% This function sorts a 2D matrix into a vector, giving the associated 2D indexes
%
% INPUTS
% - Data            --> 2D matrix      
% - mode            --> 'ascend' ascending order (default), 'descend' descending order
%
% OUTPUTS
% - A_sort          --> vector with sorted A values
% - A_ix            --> sort indexes (1st row corresponds to A lines)
%   
% try an example with this:     A = [4 2 3 ; 1 5 6 ; 8 9 7]
%
% HISTORY 
% Last version: 23/11/2011
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

if(nargin<2) mode = 'ascend'; end
    


B = reshape(A',1,[]); 
[A_sort temp_order] = sort(B,mode);

A_ix        = zeros(2,length(temp_order)); 
A_ix(1,:)   = ceil(temp_order/size(A,2)); 
A_ix(2,:)   = mod(temp_order,size(A,2)); 
A_ix(2,(find(A_ix(2,:)==0))) = size(A,2); 