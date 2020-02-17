function Index=ParametersIndex(P,name,ParamValue)

RawIndex=cellfun(@(x) isequal(x,ParamValue),{P.(name)},'UniformOutput',false);Index=[RawIndex{:}];
end