function len=structlength(X)
Ind=1;
for n = fieldnames(X)'
        name = n{1};
        len(Ind) = length(X.(name));
        Ind=Ind+1;
end
end