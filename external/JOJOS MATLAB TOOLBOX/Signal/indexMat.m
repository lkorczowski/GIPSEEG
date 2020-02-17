function [indexed_data] = indexMat(data, row_ix, col_ix)
% function [indexed_data] = indexMat(data, row_ix, col_ix)
% ************************************************************************
% Performs element-wise indexing of a 2D matrix from row and columns
% indexes. 
% --> Very useful: allows for indexing different rows with different columns!
%
% INPUTS
% - data   	--> 2D matrix (r1 x c1)
% - row_ix  --> vector with indexes of rows (r2 x 1)        
% - col_ix  --> matrix with indexes of columns for each row (r2 x c2)     
%
% OUTPUTS
% - indexed_data --> 2D indexed matrix (r2 x c2)
%
% EXAMPLE: 
% data2  = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15] 
% row   = [1 3 2]'
% col   = ([1 2 3 4; 2 3 4 5; 3 4 5 1])
% indexed_data = indexMat(data, row, col)
%
% HISTORY 
% *** 2013-09-09
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

    % assess inputs 
    if(~isrow(row_ix))
        row_ix = row_ix';
    end

    % various init
    r1 = size(data,1);
    c1 = size(data,2);
    r2 = length(row_ix);
    c2 = size(col_ix,2);

    % perform indexing
    row_ix          = repmat(row_ix,c2,1);
    col_ix          = col_ix';
    elements_ix     = sub2ind([r1 c1],row_ix(:),col_ix(:)); 
    indexed_data 	= reshape(data(elements_ix),c2,r2)';

end