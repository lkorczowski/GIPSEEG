function Cells=mat3D2cell(Matrix3D)
% convert vector cell to 3D matrix

% see also: cell2mat3D
Cells=cell(1,size(Matrix3D,3));
for indCell=1:size(Matrix3D,3)
    Cells{indCell}=Matrix3D(:,:,indCell);
end

