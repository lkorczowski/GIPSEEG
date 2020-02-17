function [out header]=gp_csv2cell(CSV,hasheader)
if nargin<2
    hasheader=0;
end
if hasheader
header=CSV(1,:);
CSV=CSV(2:end,:);
end

for Col=1:size(CSV,2)
    for Row=1:size(CSV,1)
        try
        out{Row,Col}=str2num(CSV{Row,Col});
        if isempty(out{Row,Col})
            out{Row,Col} = strsplit(CSV{Row,Col},', ');
        end
        catch
            out{Row,Col}=CSV{Row,Col};
        end
    end
end
