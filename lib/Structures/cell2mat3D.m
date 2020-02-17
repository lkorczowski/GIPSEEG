function Matrix3D=cell2mat3D(Cells)
% convert vector cell to 3D matrix

% see also: mat3D2cell
Matrix3D=[];
for indCell=1:length(Cells)
    Matrix3D(:,:,indCell)=Cells{indCell};
end

