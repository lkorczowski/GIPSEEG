function struct2var(X)
for n = fieldnames(X)'
        name = n{1};
        value = X.(name);
        assignin('caller',name,value);
end
end